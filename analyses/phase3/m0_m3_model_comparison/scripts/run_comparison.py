from __future__ import annotations

import csv
import hashlib
import math
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares

from MEA.common.analysis_io import write_json_file
from MEA.common.config import REPO_ROOT
from MEA.epcsaft_ionic.speciation_feasibility import (
    ActivitySpeciationResult,
    ActivityState,
    solve_activity_speciation,
)
from MEA.smith_missen.ideal_speciation import SPECIES_9, solve_ideal_speciation

ROOT = REPO_ROOT
ANALYSIS = ROOT / "analyses/phase3/m0_m3_model_comparison"
RESULTS = ANALYSIS / "results"
BASE_BUNDLE = (
    ROOT
    / "data/reference/epcsaft_bundles"
    / "mea-co2-h2o-nine-species-regression-input/1"
)
PREREGISTRATION = (
    ROOT
    / "analyses/phase3/ionic_epcsaft_regression"
    / "ionic_volumetric_fit_preregistration.json"
)
VLE_OBSERVATIONS = (
    ROOT
    / "data/reference/MEA/observations/vapor_liquid_equilibrium"
    / "Canonical_VLE_Observations.csv"
)
PCO2_METROLOGY = ROOT / "data/reference/MEA/manifests/pco2_metrology_manifest.csv"
TEMPERATURE_K = 313.15
PRESSURE_PA = 7326.7
MEA_MASS_FRACTION = 0.30
TRACER_LOADING = 0.466
OBSERVED_PCO2_PA = 574.0
OBSERVED_CARBAMATE_X = 0.0502
LOADING_GRID = (0.10, 0.20, 0.30, 0.40, TRACER_LOADING, 0.50, 0.60)
MODEL_ORDER = ("M0", "M1", "M2", "M3")
COMPONENT_IDS = (
    "carbon-dioxide",
    "monoethanolamine",
    "water",
    "protonated-monoethanolamine",
    "carbamate-anion",
    "bicarbonate-anion",
    "carbonate-anion",
    "hydronium-cation",
    "hydroxide-anion",
)
ION_SIGMA_START = (3.48508556586, 3.53543525721)
ION_SIGMA_CONFIRMATION = (2.53508556586, 4.48543525721)
ION_SIGMA_BOUNDS = ((2.0, 2.0), (5.8, 5.8))
ION_SIGMA_SCALE = (1.9, 1.9)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_identity(path: Path) -> dict[str, str]:
    return {
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=path, text=True
        ).strip(),
        "tree": subprocess.check_output(
            ["git", "rev-parse", "HEAD^{tree}"], cwd=path, text=True
        ).strip(),
    }


def _installed_provider_identity() -> dict[str, str]:
    distribution = metadata.distribution("epcsaft")
    record = distribution.read_text("RECORD")
    header = distribution.locate_file("epcsaft/include/epcsaft/native_sdk_v1.h")
    if record is None or not header.is_file():
        raise ValueError("installed Provider RECORD or public header is unavailable")
    return {
        "installed_record_sha256": hashlib.sha256(record.encode("utf-8")).hexdigest(),
        "installed_header_sha256": _sha256(header),
        "execution_environment": (
            "task-owned isolated wheel installation built from the recorded commit"
        ),
    }


def _rewrite_single(path: Path, overrides: dict[tuple[str, str], float]) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = tuple(rows[0])
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["component_id"], row["family"])
        if key in overrides:
            row["value"] = format(overrides[key], ".15g")
            row["source_id"] = "mea-model-comparison"
            row["locator"] = f"M0--M3 controlled override: {key[0]}::{key[1]}"
            seen.add(key)
    if seen != set(overrides):
        raise ValueError(f"single-parameter overrides were not found: {set(overrides) - seen}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _append_csv(path: Path, rows: tuple[dict[str, object], ...]) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        existing = list(csv.DictReader(handle))
        fields = tuple(existing[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows((*existing, *rows))


def _pabsch_association(bundle: Path) -> None:
    domain = "mea-tracer-313-15-k-fit-range"
    source = "pabsch-2020-induced-association"
    common = {"source_id": source, "domain_id": domain}
    _append_csv(
        bundle / "sites.csv",
        (
            {
                "record_id": "carbon-dioxide-site-a",
                "component_id": "carbon-dioxide",
                "site_id": "a",
                "site_class": "a",
                "multiplicity": 1,
                "locator": "Pabsch2020 Table 2: CO2 association scheme N=1/1",
                **common,
            },
            {
                "record_id": "carbon-dioxide-site-b",
                "component_id": "carbon-dioxide",
                "site_id": "b",
                "site_class": "b",
                "multiplicity": 1,
                "locator": "Pabsch2020 Table 2: CO2 association scheme N=1/1",
                **common,
            },
        ),
    )
    water_energy = 2425.7
    cross_energy = 0.5 * water_energy
    water_kappa = 0.04509
    co2_kappa = 0.0450
    water_sigma = 2.7927
    co2_sigma = 2.7852
    cross_kappa = math.sqrt(water_kappa * co2_kappa) * (
        math.sqrt(water_sigma * co2_sigma) / (0.5 * (water_sigma + co2_sigma))
    ) ** 3
    association_rows: list[dict[str, object]] = []
    for co2_site, water_site in (("a", "b"), ("b", "a")):
        prefix = f"carbon-dioxide-{co2_site}-water-{water_site}"
        endpoints = {
            "component_id_a": "carbon-dioxide",
            "site_id_a": co2_site,
            "component_id_b": "water",
            "site_id_b": water_site,
            **common,
        }
        association_rows.extend(
            (
                {
                    "record_id": prefix + "-association-energy-over-k",
                    "family": "association_energy_over_k",
                    "value": format(cross_energy, ".15g"),
                    "unit": "kelvin",
                    "locator": (
                        "Pabsch2020 Table 2 and Eq. 16: induced CO2--water "
                        "association; arithmetic-mean cross energy"
                    ),
                    **endpoints,
                },
                {
                    "record_id": prefix + "-association-volume",
                    "family": "association_volume",
                    "value": format(cross_kappa, ".15g"),
                    "unit": "",
                    "locator": (
                        "Pabsch2020 Table 2 kappa values with the existing "
                        "Wolbach--Sandler cross-volume rule"
                    ),
                    **endpoints,
                },
            )
        )
    _append_csv(bundle / "association.csv", tuple(association_rows))


def _prepare_bundle(
    destination: Path,
    model_id: str,
    *,
    sigma_values: tuple[float, float] | None = None,
) -> None:
    shutil.copytree(BASE_BUNDLE, destination)
    manifest = destination / "bundle.toml"
    text = manifest.read_text(encoding="utf-8")
    text = text.replace(
        'bundle_id = "mea-co2-h2o-nine-species-regression-input"',
        f'bundle_id = "mea-model-comparison-{model_id.lower()}"',
    )
    text += (
        "\n[[sources]]\n"
        'source_id = "mea-model-comparison"\n'
        'citation = "MEA-Thermodynamics controlled M0--M3 experiment"\n'
        'use_basis = "Nonpromoting parameter override for a bounded model comparison."\n'
    )
    if model_id in ("M2", "M3"):
        text += (
            "\n[[sources]]\n"
            'source_id = "pabsch-2020-induced-association"\n'
            'citation = "Pabsch, Held, and Sadowski (2020), Ind. Eng. Chem. Res. 59, 16790--16801"\n'
            'doi = "10.1021/acs.iecr.0c01888"\n'
            'use_basis = "Fixed CO2--water induced-association topology and pure association inputs."\n'
        )
    manifest.write_text(text, encoding="utf-8")

    overrides: dict[tuple[str, str], float] = {}
    if model_id == "M1":
        overrides = {
            ("protonated-monoethanolamine", "segment_diameter"): 3.0435,
            ("protonated-monoethanolamine", "dispersion_energy_over_k"): 277.174,
        }
    if sigma_values is not None:
        overrides.update(
            {
                ("protonated-monoethanolamine", "segment_diameter"): sigma_values[0],
                ("carbamate-anion", "segment_diameter"): sigma_values[1],
            }
        )
    if overrides:
        _rewrite_single(destination / "single.csv", overrides)
    if model_id in ("M2", "M3"):
        _pabsch_association(destination)


@dataclass
class Prediction:
    model_id: str
    loading: float
    temperature_k: float
    pressure_pa: float
    parameter_fingerprint: str
    success: bool
    max_abs_reaction_balance_residual: float
    pco2_pa: float
    density_kg_m3: float
    mole_fractions: dict[str, float]


class ProviderActivityEvaluator:
    def __init__(
        self,
        model: Any,
        unit_registry: Any,
        *,
        temperature_k: float,
        pressure_pa: float,
    ) -> None:
        self.model = model
        self.u = unit_registry
        self.temperature_k = temperature_k
        self.pressure_pa = pressure_pa

    def _liquid_state(self, x: np.ndarray) -> Any:
        def pressure_residual(density: float) -> float:
            state = self.model.state(
                T=self.temperature_k * self.u.kelvin,
                rho=density * self.u.mole / self.u.meter**3,
                x=tuple(float(value) for value in x),
            )
            return float(state.pressure.to("pascal").magnitude) - self.pressure_pa

        density = brentq(pressure_residual, 20_000.0, 70_000.0, xtol=1.0e-7)
        return self.model.state(
            T=self.temperature_k * self.u.kelvin,
            rho=density * self.u.mole / self.u.meter**3,
            x=tuple(float(value) for value in x),
        )

    def evaluate(
        self,
        temperature_K: float,
        pressure_Pa: float,
        mole_fractions: np.ndarray,
    ) -> ActivityState:
        if temperature_K != self.temperature_k or pressure_Pa != self.pressure_pa:
            raise ValueError("comparison evaluator state does not match the requested state")
        x = np.asarray(mole_fractions, dtype=float)
        finite = self._liquid_state(x)
        reference_x = x.copy()
        reference_x[3:] = 0.0
        reference_x /= float(np.sum(reference_x))
        reference = self._liquid_state(reference_x)
        if finite.fugacity is None or reference.fugacity is None:
            raise ValueError("Provider did not return fugacity coefficients")
        log_gamma = np.asarray(finite.fugacity.ln_coefficient) - np.asarray(
            reference.fugacity.ln_coefficient
        )
        return ActivityState(
            log_activities=np.log(np.clip(x, 1.0e-300, None)) + log_gamma,
            convention="mole_fraction_activity",
            diagnostics={
                "mapping": "ln(a_i)=ln(x_i)+ln(phi_i)-ln(phi_i,neutral-pool)",
                "density_mol_m3": float(finite.molar_density.magnitude),
            },
        )


def _molar_masses(bundle: Path) -> np.ndarray:
    with (bundle / "single.csv").open(newline="", encoding="utf-8") as handle:
        values = {
            row["component_id"]: float(row["value"])
            for row in csv.DictReader(handle)
            if row["family"] == "molar_mass"
        }
    return np.asarray([values[component] for component in COMPONENT_IDS])


def _predict(
    bundle: Path,
    model_id: str,
    loading: float,
    *,
    temperature_k: float = TEMPERATURE_K,
    pressure_pa: float = PRESSURE_PA,
) -> Prediction:
    import epcsaft

    model = epcsaft.Mixture(
        epcsaft.Parameters.from_bundle(bundle, components=COMPONENT_IDS)
    )
    evaluator = ProviderActivityEvaluator(
        model,
        epcsaft.unit_registry,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
    )
    initial = solve_ideal_speciation(
        loading, MEA_MASS_FRACTION, temperature_k
    ).mole_fractions
    solved: ActivitySpeciationResult = solve_activity_speciation(
        loading=loading,
        mea_weight_fraction=MEA_MASS_FRACTION,
        temperature_K=temperature_k,
        pressure_Pa=pressure_pa,
        evaluator=evaluator,
        initial_mole_fractions=initial,
        max_nfev=500,
    )
    state = evaluator._liquid_state(solved.mole_fractions)
    if state.fugacity is None:
        raise ValueError("Provider did not return component fugacities")
    pco2 = float(state.fugacity.value[0].to("pascal").magnitude)
    density = float(state.molar_density.magnitude) * float(
        solved.mole_fractions @ _molar_masses(bundle)
    )
    return Prediction(
        model_id=model_id,
        loading=loading,
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        parameter_fingerprint=model.parameter_fingerprint,
        success=solved.success,
        max_abs_reaction_balance_residual=solved.max_abs_residual,
        pco2_pa=pco2,
        density_kg_m3=density,
        mole_fractions=dict(zip(SPECIES_9, solved.mole_fractions.tolist(), strict=True)),
    )


def _pco2_observations() -> list[dict[str, object]]:
    with VLE_OBSERVATIONS.open(newline="", encoding="utf-8") as handle:
        canonical = list(csv.DictReader(handle))
    with PCO2_METROLOGY.open(newline="", encoding="utf-8") as handle:
        metrology = {
            row["observation_id"]: row for row in csv.DictReader(handle)
        }

    rows: list[dict[str, object]] = []
    for row in canonical:
        required = (
            "temperature_canonical_C",
            "MEA_weight_fraction",
            "CO2_loading",
            "CO2_pressure",
            "total_pressure",
        )
        if row["source_key"] != "Hilliard2008" or not all(row[key] for key in required):
            continue
        if not math.isclose(float(row["temperature_canonical_C"]), 40.0):
            continue
        if not math.isclose(float(row["MEA_weight_fraction"]), MEA_MASS_FRACTION):
            continue
        role = metrology[row["observation_id"]]
        pressure_pa = 1000.0 * float(row["total_pressure"])
        observed_pco2_pa = 1000.0 * float(row["CO2_pressure"])
        if role["target_eligible"] != "yes":
            raise ValueError(f"ineligible pCO2 row selected: {row['observation_id']}")
        if role["measurement_origin"] != "calibration_derived_partial_pressure":
            raise ValueError(f"pCO2 metrology role drift: {row['observation_id']}")
        if role["pressure_specification"] != "row_reported_total_pressure":
            raise ValueError(f"pressure role drift: {row['observation_id']}")
        if not math.isclose(float(role["state_pressure_pa"]), pressure_pa):
            raise ValueError(f"state-pressure drift: {row['observation_id']}")
        if not math.isclose(1000.0 * float(role["observed_pco2_kpa"]), observed_pco2_pa):
            raise ValueError(f"observed-pCO2 drift: {row['observation_id']}")
        rows.append(
            {
                "observation_id": row["observation_id"],
                "active_row_id": row["active_row_id"],
                "source_row": row["source_row"],
                "loading": float(row["CO2_loading"]),
                "temperature_k": TEMPERATURE_K,
                "pressure_pa": pressure_pa,
                "observed_pco2_pa": observed_pco2_pa,
                "measurement_origin": role["measurement_origin"],
                "source_locator": role["source_locator"],
            }
        )
    if len(rows) != 24:
        raise ValueError(f"expected 24 state-complete Hilliard pCO2 rows, found {len(rows)}")
    return rows


def _fit_m3(work: Path) -> tuple[np.ndarray, dict[str, Any]]:
    cache: dict[tuple[float, float], Prediction] = {}

    def prediction(values: np.ndarray) -> Prediction:
        key = tuple(float(value) for value in values)
        if key not in cache:
            bundle = work / "m3-evaluation"
            if bundle.exists():
                shutil.rmtree(bundle)
            _prepare_bundle(bundle, "M3", sigma_values=key)
            cache[key] = _predict(bundle, "M3", TRACER_LOADING)
        return cache[key]

    def residual(values: np.ndarray) -> np.ndarray:
        row = prediction(values)
        if not row.success:
            raise RuntimeError("reactive state failed during M3 fitting")
        return np.asarray(
            [
                math.log10(row.pco2_pa / OBSERVED_PCO2_PA),
                math.log10(
                    row.mole_fractions["MEACOO-"] / OBSERVED_CARBAMATE_X
                ),
            ]
        )

    fits = []
    for start in (ION_SIGMA_START, ION_SIGMA_CONFIRMATION):
        fits.append(
            least_squares(
                residual,
                np.asarray(start),
                bounds=ION_SIGMA_BOUNDS,
                x_scale=np.asarray(ION_SIGMA_SCALE),
                diff_step=2.0e-4,
                max_nfev=30,
                ftol=1.0e-9,
                xtol=1.0e-9,
                gtol=1.0e-9,
            )
        )
    primary, confirmation = fits
    singular_values = np.linalg.svd(primary.jac, compute_uv=False)
    condition = float(singular_values[0] / singular_values[-1])
    agreement_affine = float(
        np.max(np.abs((primary.x - confirmation.x) / np.asarray(ION_SIGMA_SCALE)))
    )
    diagnostics = {
        "parameters": ["MEAH+::sigma", "MEACOO-::sigma"],
        "units": ["angstrom", "angstrom"],
        "starts": [list(ION_SIGMA_START), list(ION_SIGMA_CONFIRMATION)],
        "bounds": [list(ION_SIGMA_BOUNDS[0]), list(ION_SIGMA_BOUNDS[1])],
        "primary": {
            "success": bool(primary.success),
            "message": str(primary.message),
            "values": primary.x.tolist(),
            "residuals": primary.fun.tolist(),
            "cost": float(primary.cost),
            "nfev": int(primary.nfev),
        },
        "confirmation": {
            "success": bool(confirmation.success),
            "message": str(confirmation.message),
            "values": confirmation.x.tolist(),
            "residuals": confirmation.fun.tolist(),
            "cost": float(confirmation.cost),
            "nfev": int(confirmation.nfev),
        },
        "jacobian_singular_values": singular_values.tolist(),
        "jacobian_rank": int(np.linalg.matrix_rank(primary.jac)),
        "jacobian_condition_number": condition,
        "confirmation_affine_max_delta": agreement_affine,
        "confirmation_passed": bool(
            primary.success
            and confirmation.success
            and agreement_affine <= 1.0e-6
            and np.max(np.abs(primary.fun - confirmation.fun)) <= 1.0e-8
        ),
        "promotion_allowed": False,
        "evaluation_count": len(cache),
    }
    return primary.x, _adjudicate_fit(diagnostics)


def _adjudicate_fit(diagnostics: dict[str, Any]) -> dict[str, Any]:
    values = np.asarray(diagnostics["primary"]["values"], dtype=float)
    lower = np.asarray(ION_SIGMA_BOUNDS[0], dtype=float)
    upper = np.asarray(ION_SIGMA_BOUNDS[1], dtype=float)
    scaled_bound_distance = np.minimum(values - lower, upper - values) / np.asarray(
        ION_SIGMA_SCALE
    )
    active_bound_passed = bool(np.min(scaled_bound_distance) > 1.0e-7)
    diagnostics["scaled_distance_to_nearest_bound"] = scaled_bound_distance.tolist()
    diagnostics["active_bound_passed"] = active_bound_passed
    diagnostics["tracer_accepted"] = bool(
        diagnostics["confirmation_passed"]
        and diagnostics["jacobian_rank"] == 2
        and diagnostics["jacobian_condition_number"] <= 1.0e6
        and active_bound_passed
    )
    diagnostics["acceptance_reason"] = (
        "accepted" if diagnostics["tracer_accepted"] else "rejected_active_parameter_bound"
    )
    return diagnostics


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    predictions: list[Prediction] = []
    pco2_predictions: list[dict[str, object]] = []
    pco2_observations = _pco2_observations()
    with tempfile.TemporaryDirectory(prefix="mea-m0m3-") as temporary:
        work = Path(temporary)
        bundles: dict[str, Path] = {}
        for model_id in ("M0", "M1", "M2"):
            bundle = work / model_id.lower()
            _prepare_bundle(bundle, model_id)
            bundles[model_id] = bundle
        m3_values, fit = _fit_m3(work)
        bundle = work / "m3-final"
        _prepare_bundle(bundle, "M3", sigma_values=tuple(m3_values))
        bundles["M3"] = bundle
        for model_id in MODEL_ORDER:
            for loading in LOADING_GRID:
                predictions.append(_predict(bundles[model_id], model_id, loading))
            for observation in pco2_observations:
                row = _predict(
                    bundles[model_id],
                    model_id,
                    float(observation["loading"]),
                    temperature_k=float(observation["temperature_k"]),
                    pressure_pa=float(observation["pressure_pa"]),
                )
                pco2_predictions.append(
                    {
                        "model_id": model_id,
                        **observation,
                        "predicted_pco2_pa": row.pco2_pa,
                        "success": row.success,
                        "max_abs_reaction_balance_residual": (
                            row.max_abs_reaction_balance_residual
                        ),
                        "parameter_fingerprint": row.parameter_fingerprint,
                    }
                )

    if not all(row.success for row in predictions) or not all(
        bool(row["success"]) for row in pco2_predictions
    ):
        raise RuntimeError("one or more comparison states failed")
    tracer = {
        row.model_id: row
        for row in predictions
        if math.isclose(row.loading, TRACER_LOADING)
    }
    summary_rows = []
    for model_id in MODEL_ORDER:
        row = tracer[model_id]
        summary_rows.append(
            {
                "model_id": model_id,
                "pco2_observed_pa": OBSERVED_PCO2_PA,
                "pco2_predicted_pa": row.pco2_pa,
                "pco2_log10_residual": math.log10(row.pco2_pa / OBSERVED_PCO2_PA),
                "meacoo_observed_mole_fraction": OBSERVED_CARBAMATE_X,
                "meacoo_predicted_mole_fraction": row.mole_fractions["MEACOO-"],
                "meacoo_log10_residual": math.log10(
                    row.mole_fractions["MEACOO-"] / OBSERVED_CARBAMATE_X
                ),
                "density_predicted_kg_m3": row.density_kg_m3,
                "max_abs_reaction_balance_residual": row.max_abs_reaction_balance_residual,
                "parameter_fingerprint": row.parameter_fingerprint,
                "model_status": (
                    "bound_limited_diagnostic"
                    if model_id == "M3" and not fit["tracer_accepted"]
                    else "evaluated"
                ),
            }
        )
    _write_csv(RESULTS / "model_summary.csv", summary_rows)
    species_rows = []
    loading_rows = []
    for row in predictions:
        loading_rows.append(
            {
                "model_id": row.model_id,
                "loading_mol_co2_per_mol_mea": row.loading,
                "temperature_k": row.temperature_k,
                "pressure_pa": row.pressure_pa,
                "pco2_pa": row.pco2_pa,
                "density_kg_m3": row.density_kg_m3,
                "max_abs_reaction_balance_residual": row.max_abs_reaction_balance_residual,
            }
        )
        for species, value in row.mole_fractions.items():
            species_rows.append(
                {
                    "model_id": row.model_id,
                    "loading_mol_co2_per_mol_mea": row.loading,
                    "temperature_k": row.temperature_k,
                    "pressure_pa": row.pressure_pa,
                    "species": species,
                    "mole_fraction": value,
                }
            )
    _write_csv(RESULTS / "loading_predictions.csv", loading_rows)
    _write_csv(RESULTS / "species_predictions.csv", species_rows)
    _write_csv(RESULTS / "pco2_loading_comparison.csv", pco2_predictions)

    import epcsaft

    receipt = {
        "schema_version": 1,
        "experiment": "m0_m3_reactive_model_comparison",
        "status": (
            "COMPLETED_NONPROMOTING_EXPERIMENT"
            if fit["tracer_accepted"]
            else "COMPLETED_REJECTED_BOUND_LIMITED_M3_EXPERIMENT"
        ),
        "manuscript_modified": False,
        "scientific_claim_boundary": (
            "M3 is fitted to two training observations at one state; no predictive, "
            "identifiability, validation, or parameter-promotion claim is admitted."
        ),
        "conditions": {
            "temperature_k": TEMPERATURE_K,
            "pressure_pa": PRESSURE_PA,
            "mea_mass_fraction_unloaded": MEA_MASS_FRACTION,
            "loading_grid": list(LOADING_GRID),
        },
        "models": {
            "M0": "current authority-neutral regression-input bundle",
            "M1": "M0 with neutral-MEA sigma and epsilon transferred to MEAH+",
            "M2": "M0 with fixed Pabsch2020 CO2--water induced association",
            "M3": "M2 with two fitted ion segment diameters",
        },
        "m3_fit": fit,
        "density_context": {
            "status": "PLOTTED_NOT_FITTED",
            "reason": "Amundsen source rows do not report pressure in the admitted table.",
        },
        "pco2_loading_comparison": {
            "source": "Hilliard2008",
            "observation_count": len(pco2_observations),
            "loading_range": [
                pco2_observations[0]["loading"],
                pco2_observations[-1]["loading"],
            ],
            "state_policy": "row_reported_total_pressure",
            "measurement_origin": "calibration_derived_partial_pressure",
            "model_count": len(MODEL_ORDER),
        },
        "provider": {
            "version": getattr(epcsaft, "__version__", "0.2.0.dev0"),
            "module": "epcsaft/__init__.py",
            **_git_identity(
                ROOT.parent / "ePC-SAFT-project/ePC-SAFT-eos"
            ),
            **_installed_provider_identity(),
        },
        "source_hashes": {
            str(BASE_BUNDLE.relative_to(ROOT) / name): _sha256(BASE_BUNDLE / name)
            for name in (
                "bundle.toml",
                "single.csv",
                "sites.csv",
                "association.csv",
            )
        }
        | {
            str(PREREGISTRATION.relative_to(ROOT)): _sha256(PREREGISTRATION),
            str(VLE_OBSERVATIONS.relative_to(ROOT)): _sha256(VLE_OBSERVATIONS),
            str(PCO2_METROLOGY.relative_to(ROOT)): _sha256(PCO2_METROLOGY),
        },
        "outputs": {
            name: _sha256(RESULTS / name)
            for name in (
                "model_summary.csv",
                "species_predictions.csv",
                "loading_predictions.csv",
                "pco2_loading_comparison.csv",
            )
        },
    }
    write_json_file(RESULTS / "comparison_receipt.json", receipt)
    print(RESULTS / "comparison_receipt.json")


if __name__ == "__main__":
    main()
