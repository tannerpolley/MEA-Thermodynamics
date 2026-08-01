from __future__ import annotations

import json
import math
import tempfile
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np

import run_comparison as comparison
import run_pressure_sensitivity_fit as pressure_fit
from MEA.epcsaft_ionic.speciation_feasibility import (
    ActivitySpeciationResult,
    ActivityState,
    solve_activity_speciation,
)
from MEA.smith_missen.ideal_speciation import solve_ideal_speciation

RESULTS = comparison.RESULTS
PROVIDER_SOURCE = comparison.ROOT.parent / "ePC-SAFT-project/ePC-SAFT-eos"
REFERENCE_FLOOR = 1.0e-12
REFERENCE_ION_MULTIPLIERS = np.asarray((2.5, 1.0, 1.0, 1.0, 2.5, 1.0))
MODEL_ORDER = ("M0", "M1", "M2", "M3", "M4A", "M4B", "M5Q", "M5")
GROSS_CO2 = {
    "segment_count": 1.5131,
    "segment_diameter": 3.1869,
    "dispersion_energy_over_k": 163.33,
    "quadrupole_moment": 4.4,
}
NEUTRAL_DIPOLES = {
    "water": 1.8546,
    "monoethanolamine": 2.27,
}


MODELS = {
    "M0": ("base", "fixed regression input"),
    "M1": ("neutral-parent transfer", "controlled diagnostic"),
    "M2": ("induced association", "fixed literature topology"),
    "M3": ("two-row ionic fit", "bound-limited diagnostic"),
    "M4A": ("Hilliard pressure fit", "nonpromoting fit"),
    "M4B": ("Hilliard + Jou pressure fit", "nonpromoting fit"),
    "M5Q": ("M2 + CO2 QQ", "component-consistent CO2 diagnostic"),
    "M5": ("M5Q + H2O/MEA DD + DQ", "unrefitted full-polar diagnostic"),
}

_rows = pressure_fit._rows
_write_csv = pressure_fit._write_csv
_metric = pressure_fit._metrics


def _append_source(manifest: Path, source_id: str, citation: str, use_basis: str) -> None:
    text = manifest.read_text(encoding="utf-8")
    text += (
        "\n[[sources]]\n"
        f'source_id = "{source_id}"\n'
        f'citation = "{citation}"\n'
        f'use_basis = "{use_basis}"\n'
    )
    manifest.write_text(text, encoding="utf-8")


def _prepare_polar_bundle(destination: Path, *, full: bool) -> None:
    comparison._prepare_bundle(destination, "M2")
    manifest = destination / "bundle.toml"
    for source in (
        (
            "gross-vrabec-polar-equations",
            "Gross (2005), AIChE J. 51, 2556; Gross and Vrabec (2006), AIChE J. 52, 1194; Vrabec and Gross (2009), J. Phys. Chem. B 113, 10935",
            "Fixed Gross--Vrabec DD/QQ/DQ35 equations implemented by the recorded Provider commit.",
        ),
        (
            "gross-2005-co2-qq",
            "Gross (2005), doi:10.1002/aic.10502, Table A1 first CO2 row",
            "CO2 pure parameters fitted with QQ active and measured Q=4.4 D angstrom.",
        ),
    ):
        _append_source(manifest, *source)
    comparison._rewrite_single(
        destination / "single.csv",
        {
            ("carbon-dioxide", "segment_count"): GROSS_CO2["segment_count"],
            ("carbon-dioxide", "segment_diameter"): GROSS_CO2["segment_diameter"],
            ("carbon-dioxide", "dispersion_energy_over_k"): GROSS_CO2[
                "dispersion_energy_over_k"
            ],
        },
    )
    gross_families = set(GROSS_CO2) - {"quadrupole_moment"}
    single_rows = _rows(destination / "single.csv")
    for row in single_rows:
        if row["component_id"] == "carbon-dioxide" and row["family"] in gross_families:
            row["source_id"] = "gross-2005-co2-qq"
            row["locator"] = "Gross (2005), Table A1, first carbon-dioxide row"
    _write_csv(destination / "single.csv", single_rows)
    cross_kappa = math.sqrt(0.04509 * 0.0450) * (
        math.sqrt(2.7927 * GROSS_CO2["segment_diameter"])
        / (0.5 * (2.7927 + GROSS_CO2["segment_diameter"]))
    ) ** 3
    association_rows = _rows(destination / "association.csv")
    for row in association_rows:
        if (
            row["source_id"] == "pabsch-2020-induced-association"
            and row["family"] == "association_volume"
        ):
            row["value"] = format(cross_kappa, ".15g")
            row["locator"] += "; recomputed with Gross (2005) CO2 segment diameter"
    _write_csv(destination / "association.csv", association_rows)
    polar_rows = [
        {
            "record_id": "carbon-dioxide-quadrupole-moment",
            "component_id": "carbon-dioxide",
            "family": "quadrupole_moment",
            "value": GROSS_CO2["quadrupole_moment"],
            "unit": "debye * angstrom",
            "source_id": "gross-2005-co2-qq",
            "locator": "Table A1, first carbon-dioxide row",
            "domain_id": "mea-tracer-313-15-k-fit-range",
        }
    ]
    if full:
        _append_source(
            destination / "bundle.toml",
            "neutral-dipole-diagnostic",
            "Clough et al. (1973), J. Chem. Phys. 59, 2254; Tripathi (2016), doi:10.5821/dissertation-2117-106297, chapter 7",
            "Fixed physical gas-phase moments; no moment or H2O/MEA PC-SAFT parameter was fitted here.",
        )
        polar_rows.extend(
            {
                "record_id": f"{component_id}-dipole-moment",
                "component_id": component_id,
                "family": "dipole_moment",
                "value": value,
                "unit": "debye",
                "source_id": "neutral-dipole-diagnostic",
                "locator": (
                    "Clough et al. (1973): rotationless H2O moment"
                    if component_id == "water"
                    else "stable gas-phase MEA conformer; single-conformer diagnostic"
                ),
                "domain_id": "mea-tracer-313-15-k-fit-range",
            }
            for component_id, value in NEUTRAL_DIPOLES.items()
        )
    comparison._append_csv(destination / "single.csv", tuple(polar_rows))
    comparison._append_csv(
        destination / "model.csv",
        (
            {
                "record_id": "gross-vrabec-point-multipole",
                "family": "polar_formulation",
                "value": "gross-vrabec-point-multipole",
                "unit": "",
                "source_id": "gross-vrabec-polar-equations",
                "locator": "Gross (2005); Gross and Vrabec (2006); Vrabec and Gross (2009)",
                "domain_id": "mea-tracer-313-15-k-fit-range",
            },
        ),
    )


def _m4_values(model_id: str) -> dict[tuple[str, str], float]:
    rows = {row["parameter_identity"]: row for row in _rows(RESULTS / "pressure_fit_parameters.csv")}
    prefix = model_id.lower()
    return {
        ("protonated-monoethanolamine", "segment_diameter"): float(
            rows["MEAH+::sigma"][f"{prefix}_fitted_value"]
        ),
        ("protonated-monoethanolamine", "dispersion_energy_over_k"): float(
            rows["MEAH+::epsilon_over_k"][f"{prefix}_fitted_value"]
        ),
    }


def _prepare_bundle(destination: Path, model_id: str) -> None:
    if model_id in {"M5Q", "M5"}:
        _prepare_polar_bundle(destination, full=model_id == "M5")
        return
    if model_id == "M3":
        receipt = json.loads((RESULTS / "comparison_receipt.json").read_text(encoding="utf-8"))
        values = tuple(float(value) for value in receipt["m3_fit"]["primary"]["values"])
        comparison._prepare_bundle(destination, model_id, sigma_values=values)
        return
    comparison._prepare_bundle(destination, model_id if model_id in {"M0", "M1", "M2"} else "M2")
    if model_id in {"M4A", "M4B"}:
        comparison._rewrite_single(destination / "single.csv", _m4_values(model_id))


class PolarActivityEvaluator(comparison.ProviderActivityEvaluator):
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
        reference_x[3:] = REFERENCE_FLOOR * REFERENCE_ION_MULTIPLIERS
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
                "mapping": "ln(a_i)=ln(x_i)+ln(phi_i)-ln(phi_i,neutral-pool limit)",
                "neutral_pool_ionic_mole_fraction_floor": REFERENCE_FLOOR,
                "neutral_pool_ion_multipliers": REFERENCE_ION_MULTIPLIERS.tolist(),
                "density_mol_m3": float(finite.molar_density.magnitude),
            },
        )


def _evaluate(
    model: Any,
    unit_registry: Any,
    model_id: str,
    row: dict[str, object],
    bundle: Path,
) -> dict[str, object]:
    evaluator_type = PolarActivityEvaluator if model_id in {"M5Q", "M5"} else comparison.ProviderActivityEvaluator
    evaluator = evaluator_type(
        model,
        unit_registry,
        temperature_k=float(row["temperature_k"]),
        pressure_pa=float(row["pressure_pa"]),
    )
    initial = solve_ideal_speciation(
        float(row["loading"]),
        float(row["mea_mass_fraction"]),
        float(row["temperature_k"]),
    ).mole_fractions
    solved: ActivitySpeciationResult = solve_activity_speciation(
        loading=float(row["loading"]),
        mea_weight_fraction=float(row["mea_mass_fraction"]),
        temperature_K=float(row["temperature_k"]),
        pressure_Pa=float(row["pressure_pa"]),
        evaluator=evaluator,
        initial_mole_fractions=initial,
        max_nfev=500,
    )
    if not solved.success:
        raise RuntimeError(f"reactive solve failed for {model_id}:{row['observation_id']}")
    state = evaluator._liquid_state(solved.mole_fractions)
    if state.fugacity is None:
        raise ValueError("Provider did not return fugacity")
    predicted = float(state.fugacity.value[0].to("pascal").magnitude)
    molar_masses = comparison._molar_masses(bundle)
    density = float(state.molar_density.magnitude) * float(solved.mole_fractions @ molar_masses)
    return {
        "model_id": model_id,
        "observation_id": row["observation_id"],
        "source_key": row["source_key"],
        "mea_mass_fraction": row["mea_mass_fraction"],
        "loading_mol_co2_per_mol_mea": row["loading"],
        "temperature_k": row["temperature_k"],
        "state_pressure_pa": row["pressure_pa"],
        "observed_pco2_pa": row["observed_pco2_pa"],
        "predicted_pco2_pa": predicted,
        "log10_residual": math.log10(predicted / float(row["observed_pco2_pa"])),
        "density_kg_m3": density,
        "ares_polar": state.ares_polar,
        "split": row["split"],
        "role": row["role"],
        "group_id": row["group_id"],
        "measurement_origin": row["measurement_origin"],
        "success": str(solved.success).lower(),
        "max_abs_reaction_balance_residual": solved.max_abs_residual,
        "parameter_fingerprint": model.parameter_fingerprint,
    }


def _retained_hilliard_rows(
    executable_by_id: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for retained in _rows(RESULTS / "pco2_loading_comparison.csv"):
        if retained["model_id"] not in {"M0", "M1", "M3"}:
            continue
        source = executable_by_id[retained["observation_id"]]
        result.append(
            {
                "model_id": retained["model_id"],
                "observation_id": retained["observation_id"],
                "source_key": source["source_key"],
                "mea_mass_fraction": source["mea_mass_fraction"],
                "loading_mol_co2_per_mol_mea": retained["loading"],
                "temperature_k": retained["temperature_k"],
                "state_pressure_pa": retained["pressure_pa"],
                "observed_pco2_pa": retained["observed_pco2_pa"],
                "predicted_pco2_pa": retained["predicted_pco2_pa"],
                "log10_residual": math.log10(
                    float(retained["predicted_pco2_pa"])
                    / float(retained["observed_pco2_pa"])
                ),
                "density_kg_m3": "",
                "ares_polar": 0.0,
                "split": source["split"],
                "role": source["role"],
                "group_id": source["group_id"],
                "measurement_origin": retained["measurement_origin"],
                "success": retained["success"].lower(),
                "max_abs_reaction_balance_residual": retained[
                    "max_abs_reaction_balance_residual"
                ],
                "parameter_fingerprint": retained["parameter_fingerprint"],
            }
        )
    return result


def _retained_pressure_fit_rows() -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for retained in _rows(RESULTS / "pressure_fit_predictions.csv"):
        result.append(
            {
                **retained,
                "density_kg_m3": "",
                "ares_polar": 0.0,
                "success": retained["success"].lower(),
            }
        )
    return result


def main() -> None:
    import epcsaft

    _, executable = pressure_fit._temperature_inventory()
    executable_by_id = {str(row["observation_id"]): row for row in executable}
    predictions = [
        *_retained_hilliard_rows(executable_by_id),
        *_retained_pressure_fit_rows(),
    ]
    fingerprints: dict[str, str] = {}
    replay_checks: dict[str, dict[str, float | str]] = {}
    new_state_evaluations = 0
    jou = [row for row in executable if row["source_key"] == "Jou1995"]
    first_hilliard = next(
        row
        for row in executable
        if row["source_key"] == "Hilliard2008"
        and math.isclose(float(row["mea_mass_fraction"]), 0.30)
    )
    with tempfile.TemporaryDirectory(prefix="mea-polar-comparison-") as temporary:
        root = Path(temporary)
        for model_id in MODEL_ORDER:
            bundle = root / model_id.lower()
            _prepare_bundle(bundle, model_id)
            model = epcsaft.Mixture(
                epcsaft.Parameters.from_bundle(bundle, components=comparison.COMPONENT_IDS)
            )
            fingerprints[model_id] = model.parameter_fingerprint
            retained_fingerprints = {
                str(row["parameter_fingerprint"])
                for row in predictions
                if row["model_id"] == model_id
            }
            if retained_fingerprints and retained_fingerprints != {model.parameter_fingerprint}:
                raise ValueError(f"retained parameter fingerprint drift for {model_id}")
            if model_id in {"M0", "M2"}:
                replay = _evaluate(
                    model,
                    epcsaft.unit_registry,
                    model_id,
                    first_hilliard,
                    bundle,
                )
                retained = next(
                    row
                    for row in predictions
                    if row["model_id"] == model_id
                    and row["observation_id"] == first_hilliard["observation_id"]
                )
                replay_checks[model_id] = {
                    "observation_id": str(first_hilliard["observation_id"]),
                    "retained_predicted_pco2_pa": float(retained["predicted_pco2_pa"]),
                    "merged_eos_predicted_pco2_pa": float(replay["predicted_pco2_pa"]),
                    "absolute_difference_pa": abs(
                        float(retained["predicted_pco2_pa"])
                        - float(replay["predicted_pco2_pa"])
                    ),
                }
                if not math.isclose(
                    float(retained["predicted_pco2_pa"]),
                    float(replay["predicted_pco2_pa"]),
                    rel_tol=2.0e-12,
                    abs_tol=1.0e-10,
                ):
                    raise ValueError(f"merged EOS changed the retained nonpolar result for {model_id}")
                new_state_evaluations += 1
            rows_to_evaluate = (
                jou
                if model_id in {"M0", "M1", "M3"}
                else executable
                if model_id in {"M5Q", "M5"}
                else ()
            )
            for index, row in enumerate(rows_to_evaluate, start=1):
                predictions.append(_evaluate(model, epcsaft.unit_registry, model_id, row, bundle))
                new_state_evaluations += 1
                if index % 8 == 0 or index == len(rows_to_evaluate):
                    print(f"{model_id}: evaluated {index}/{len(rows_to_evaluate)} states", flush=True)

    selections = (
        ("all-30wt", lambda row: math.isclose(float(row["mea_mass_fraction"]), 0.30)),
        ("Hilliard-30wt", lambda row: row["source_key"] == "Hilliard2008" and math.isclose(float(row["mea_mass_fraction"]), 0.30)),
        ("Jou-30wt", lambda row: row["source_key"] == "Jou1995"),
    )
    metrics = [
        _metric(predictions, model_id, selection, predicate)
        for model_id in MODEL_ORDER
        for selection, predicate in selections
    ]
    transfer_models = ("M2", "M4A", "M4B", "M5Q", "M5")
    metrics.extend(
        _metric(predictions, model_id, selection, predicate)
        for model_id in transfer_models
        for selection, predicate in (
            (
                "untouched-17-and-40wt",
                lambda row: not math.isclose(float(row["mea_mass_fraction"]), 0.30),
            ),
            ("all-available", lambda row: True),
        )
    )
    model_rows = [
        {
            "model_id": model_id,
            "label": MODELS[model_id][0],
            "parameter_status": MODELS[model_id][1],
            "parameter_fingerprint": fingerprints[model_id],
            "polar_active": str(model_id in {"M5Q", "M5"}).lower(),
            "co2_quadrupole_debye_angstrom": GROSS_CO2["quadrupole_moment"] if model_id in {"M5Q", "M5"} else "",
            "water_dipole_debye": NEUTRAL_DIPOLES["water"] if model_id == "M5" else "",
            "mea_dipole_debye": NEUTRAL_DIPOLES["monoethanolamine"] if model_id == "M5" else "",
        }
        for model_id in MODEL_ORDER
    ]
    outputs = {
        "m0_m5_pressure_predictions.csv": predictions,
        "m0_m5_pressure_metrics.csv": metrics,
        "m0_m5_model_definitions.csv": model_rows,
    }
    for name, rows in outputs.items():
        _write_csv(RESULTS / name, rows)

    audit = {
        "analysis_role": "executed nested polar-model comparison",
        "status": "M5_EXECUTED_DIAGNOSTIC_ONLY",
        "eos": {
            "version": metadata.version("epcsaft"),
            "build_source": comparison._git_identity(PROVIDER_SOURCE),
            **comparison._installed_provider_identity(),
            "polar_formulation": "gross-vrabec-point-multipole",
            "public_output": "State.ares_polar = a_DD + a_QQ + a_DQ35",
        },
        "models": {
            "M5Q": {
                "definition": "M2 induced association plus the Gross (2005) QQ-consistent CO2 pure set",
                "parameters": GROSS_CO2,
                "scientific_status": "component-consistent CO2 quadrupolar diagnostic",
            },
            "M5": {
                "definition": "M5Q plus physical neutral H2O and MEA dipoles, activating DD and DQ",
                "dipole_moments_debye": NEUTRAL_DIPOLES,
                "scientific_status": "full-equation diagnostic; H2O and MEA PC-SAFT/association parameters were not refitted with DD active",
            },
        },
        "reference_boundary_limit": {
            "reason": "the Provider rejects exact composition-boundary polar states",
            "ionic_mole_fraction_floor": REFERENCE_FLOOR,
            "electroneutral_ion_multipliers": REFERENCE_ION_MULTIPLIERS.tolist(),
            "role": "numerical approximation to the existing neutral-pool reference limit; not a fitted parameter",
        },
        "literature_evidence": {
            "co2_qq_source": "Gross (2005), doi:10.1002/aic.10502, Table A1 first CO2 row",
            "co2_qq_source_sha256": "cfd73fa388bec051b91e4f51e1c2e16ddf2f2140b6c037996ab33290887c5d50",
            "dd_equation_source": "Gross and Vrabec (2006), doi:10.1002/aic.10683",
            "dd_equation_source_sha256": "0f342632efca9425db2a2cc26bac652605ae8c787900f3864d060c59d8e2a2e3",
            "dq_equation_source": "Vrabec and Gross (2009), doi:10.1021/jp072619u",
            "dq_equation_source_sha256": "c58ae9e065580bf9cd80a26f3e8b8eb875b95ef11cb5004e60c42586b507d9a7",
            "water_dipole": "1.8546 D rotationless gas-phase value from Clough et al. (1973)",
            "mea_dipole": "2.27 D stable gas-phase conformer from Tripathi (2016), chapter 7; liquid conformer population remains unresolved",
        },
        "conditions": {
            "temperature_k": 313.15,
            "complete_30wt_rows_per_model": 32,
            "additional_17_and_40wt_rows_for_transfer_models": 12,
            "model_count": len(MODEL_ORDER),
            "new_state_evaluation_count": new_state_evaluations,
            "unified_prediction_row_count": len(predictions),
        },
        "metric_highlights": {
            model_id: {
                selection: next(
                    float(row["log10_rmse"])
                    for row in metrics
                    if row["model_id"] == model_id and row["selection"] == selection
                )
                for selection in (
                    ("all-30wt", "untouched-17-and-40wt")
                    if model_id in {"M2", "M4A", "M4B", "M5Q", "M5"}
                    else ("all-30wt",)
                )
            }
            for model_id in MODEL_ORDER
        },
        "retained_nonpolar_replay_checks": replay_checks,
        "parameter_fingerprints": fingerprints,
        "source_hashes": {
            str(path.relative_to(comparison.ROOT)): comparison._sha256(path)
            for path in (
                comparison.VLE_OBSERVATIONS,
                comparison.PCO2_METROLOGY,
                pressure_fit.SPLIT_MANIFEST,
                RESULTS / "comparison_receipt.json",
                RESULTS / "pressure_fit_parameters.csv",
            )
        },
        "outputs": {name: comparison._sha256(RESULTS / name) for name in outputs},
        "claim_boundary": "The polar equations are executable and compared against admitted pressure data. M5Q is a nested CO2 diagnostic. M5 is not an internally refitted polar parameterization and cannot support parameter promotion or a predictive manuscript claim.",
    }
    (RESULTS / "m5_polar_capability_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
