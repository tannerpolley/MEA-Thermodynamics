from __future__ import annotations

import csv
import json
import math
import shutil
import tempfile
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import qr
from scipy.optimize import lsq_linear

import run_comparison as comparison
from MEA.epcsaft_ionic.speciation_feasibility import (
    ActivitySpeciationResult,
    solve_activity_speciation,
)
from MEA.smith_missen.ideal_speciation import solve_ideal_speciation

ROOT = comparison.ROOT
ANALYSIS = comparison.ANALYSIS
RESULTS = comparison.RESULTS
SPLIT_MANIFEST = ROOT / "data/reference/MEA/manifests/grouped_split_manifest.csv"
PROVIDER_SOURCE = ROOT.parent / "ePC-SAFT-project/ePC-SAFT-eos"
TEMPERATURE_C = 40.0
TEMPERATURE_K = 313.15
SCREEN_ROW_COUNT = 6
IDENTIFIABLE_SINGULAR_RATIO = 0.05
FIT_DIFFERENCE_AFFINE = 0.005
FIT_TRUST_RADIUS_AFFINE = 0.5
FIT_MAX_ITERATIONS = 2
FIT_CONDITION_LIMIT = 100.0


@dataclass(frozen=True)
class Candidate:
    identity: str
    label: str
    kind: str
    record: tuple[str, str] | str
    start: float
    screen_step: float
    screen_scale: float
    unit: str
    fit_eligible: bool
    lower: float | None = None
    upper: float | None = None
    affine_scale: float | None = None


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty result table: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _single_values() -> dict[tuple[str, str], float]:
    return {
        (row["component_id"], row["family"]): float(row["value"])
        for row in _rows(comparison.BASE_BUNDLE / "single.csv")
    }


def _pair_values() -> dict[str, float]:
    return {
        row["record_id"]: float(row["value"])
        for row in _rows(comparison.BASE_BUNDLE / "pair.csv")
    }


def _candidates() -> tuple[Candidate, ...]:
    single = _single_values()
    pair = _pair_values()
    preregistration = json.loads(comparison.PREREGISTRATION.read_text(encoding="utf-8"))
    active = {row["identity"]: row for row in preregistration["active_coordinates"]}
    active_records = {
        "MEAH+::sigma": ("protonated-monoethanolamine", "segment_diameter"),
        "MEAH+::epsilon_over_k": (
            "protonated-monoethanolamine",
            "dispersion_energy_over_k",
        ),
        "MEACOO-::sigma": ("carbamate-anion", "segment_diameter"),
    }
    labels = {
        "MEAH+::sigma": r"MEAH+ sigma",
        "MEAH+::epsilon_over_k": r"MEAH+ epsilon/k",
        "MEACOO-::sigma": r"MEACOO- sigma",
    }
    result: list[Candidate] = []
    for identity, record in active_records.items():
        contract = active[identity]
        start = single[record]
        if not math.isclose(start, float(contract["start"]), rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"preregistered start drift for {identity}")
        scale = float(contract["affine_scale"])
        result.append(
            Candidate(
                identity=identity,
                label=labels[identity],
                kind="single",
                record=record,
                start=start,
                screen_step=FIT_DIFFERENCE_AFFINE * scale,
                screen_scale=scale,
                unit=str(contract["unit"]),
                fit_eligible=True,
                lower=float(contract["bounds"][0]),
                upper=float(contract["bounds"][1]),
                affine_scale=scale,
            )
        )

    diagnostic_single = (
        (
            "MEAH+::dBorn",
            "MEAH+ dBorn",
            ("protonated-monoethanolamine", "born_diameter"),
            "angstrom",
        ),
        (
            "MEACOO-::epsilon_over_k",
            "MEACOO- epsilon/k",
            ("carbamate-anion", "dispersion_energy_over_k"),
            "K",
        ),
        (
            "MEACOO-::dBorn",
            "MEACOO- dBorn",
            ("carbamate-anion", "born_diameter"),
            "angstrom",
        ),
    )
    for identity, label, record, unit in diagnostic_single:
        start = single[record]
        result.append(
            Candidate(
                identity=identity,
                label=label,
                kind="single",
                record=record,
                start=start,
                screen_step=0.005 * start,
                screen_scale=start,
                unit=unit,
                fit_eligible=False,
            )
        )

    diagnostic_pairs = (
        ("kij::CO2-H2O", "kij CO2-H2O", "kij-carbon-dioxide-water"),
        ("kij::CO2-MEA", "kij CO2-MEA", "kij-carbon-dioxide-monoethanolamine"),
        (
            "kij::CO2-MEAH+",
            "kij CO2-MEAH+",
            "kij-carbon-dioxide-protonated-monoethanolamine",
        ),
        (
            "kij::CO2-MEACOO-",
            "kij CO2-MEACOO-",
            "kij-carbon-dioxide-carbamate-anion",
        ),
        (
            "kij::MEAH+-MEACOO-",
            "kij MEAH+-MEACOO-",
            "kij-protonated-monoethanolamine-carbamate-anion",
        ),
    )
    for identity, label, record in diagnostic_pairs:
        result.append(
            Candidate(
                identity=identity,
                label=label,
                kind="pair",
                record=record,
                start=pair[record],
                screen_step=0.002,
                screen_scale=1.0,
                unit="dimensionless",
                fit_eligible=False,
            )
        )
    return tuple(result)


def _rewrite_pair(path: Path, overrides: dict[str, float]) -> None:
    rows = _rows(path)
    fields = tuple(rows[0])
    seen: set[str] = set()
    for row in rows:
        record_id = row["record_id"]
        if record_id in overrides:
            row["value"] = format(overrides[record_id], ".15g")
            row["source_id"] = "mea-model-comparison"
            row["locator"] = f"313.15 K pressure sensitivity: {record_id}"
            seen.add(record_id)
    if seen != set(overrides):
        raise ValueError(f"pair overrides not found: {set(overrides) - seen}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _prepare_bundle(
    destination: Path,
    candidates: tuple[Candidate, ...],
    values: dict[str, float],
) -> None:
    comparison._prepare_bundle(destination, "M2")
    single_overrides = {
        candidate.record: values[candidate.identity]
        for candidate in candidates
        if candidate.kind == "single"
        and not math.isclose(values[candidate.identity], candidate.start)
    }
    pair_overrides = {
        str(candidate.record): values[candidate.identity]
        for candidate in candidates
        if candidate.kind == "pair"
        and not math.isclose(values[candidate.identity], candidate.start)
    }
    if single_overrides:
        comparison._rewrite_single(destination / "single.csv", single_overrides)
    if pair_overrides:
        _rewrite_pair(destination / "pair.csv", pair_overrides)


def _temperature_inventory() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    canonical = _rows(comparison.VLE_OBSERVATIONS)
    metrology = {row["observation_id"]: row for row in _rows(comparison.PCO2_METROLOGY)}
    splits = {
        row["record_id"]: row
        for row in _rows(SPLIT_MANIFEST)
        if row["target_family"] == "vle_pressure"
    }
    inventory: list[dict[str, object]] = []
    executable: list[dict[str, object]] = []
    for row in canonical:
        reported_temperature = (
            row["temperature_canonical_C"] or row["temperature_reported_C"]
        )
        if not reported_temperature:
            continue
        if not math.isclose(float(reported_temperature), TEMPERATURE_C):
            continue
        if not row["CO2_loading"] or not row["CO2_pressure"]:
            continue
        role = metrology[row["observation_id"]]
        split = splits.get(row["active_row_id"]) or splits.get(row["observation_id"])
        eligible = bool(role["target_eligible"] == "yes" and role["state_pressure_pa"])
        inventory.append(
            {
                "observation_id": row["observation_id"],
                "source_key": row["source_key"],
                "mea_mass_fraction": row["MEA_weight_fraction"],
                "loading_mol_co2_per_mol_mea": row["CO2_loading"],
                "observed_pco2_kpa": row["CO2_pressure"],
                "state_pressure_pa": role["state_pressure_pa"],
                "measurement_origin": role["measurement_origin"],
                "target_eligible": role["target_eligible"],
                "eligibility_reason": role["eligibility_reason"],
                "split": split["split"] if split else "unassigned",
                "role": split["role"] if split else "unassigned",
                "group_id": split["group_id"] if split else "unassigned",
                "source_locator": role["source_locator"],
            }
        )
        if not eligible:
            continue
        if split is None:
            raise ValueError(f"eligible row has no frozen split: {row['observation_id']}")
        if role["measurement_origin"] != "calibration_derived_partial_pressure":
            raise ValueError(f"unexpected pCO2 origin: {row['observation_id']}")
        if role["pressure_specification"] != "row_reported_total_pressure":
            raise ValueError(f"unexpected state-pressure role: {row['observation_id']}")
        if not row["total_pressure"] or not math.isclose(
            float(role["state_pressure_pa"]),
            1000.0 * float(row["total_pressure"]),
        ):
            raise ValueError(f"state-pressure unit drift: {row['observation_id']}")
        if not math.isclose(
            float(role["observed_pco2_kpa"]), float(row["CO2_pressure"])
        ):
            raise ValueError(f"observed-pCO2 unit drift: {row['observation_id']}")
        executable.append(
            {
                "observation_id": row["observation_id"],
                "source_key": row["source_key"],
                "mea_mass_fraction": float(row["MEA_weight_fraction"]),
                "loading": float(row["CO2_loading"]),
                "temperature_k": TEMPERATURE_K,
                "pressure_pa": float(role["state_pressure_pa"]),
                "observed_pco2_pa": 1000.0 * float(row["CO2_pressure"]),
                "measurement_origin": role["measurement_origin"],
                "split": split["split"],
                "role": split["role"],
                "group_id": split["group_id"],
                "source_locator": role["source_locator"],
            }
        )
    group_counts: dict[tuple[str, float, str], int] = {}
    for row in executable:
        key = (str(row["source_key"]), float(row["mea_mass_fraction"]), str(row["split"]))
        group_counts[key] = group_counts.get(key, 0) + 1
    expected = {
        ("Hilliard2008", 0.17, "validation"): 6,
        ("Hilliard2008", 0.30, "training"): 24,
        ("Hilliard2008", 0.40, "validation"): 6,
        ("Jou1995", 0.30, "validation"): 8,
    }
    if group_counts != expected:
        raise ValueError(f"313.15 K executable pressure partition drift: {group_counts}")
    return inventory, executable


class Experiment:
    def __init__(self, root: Path, candidates: tuple[Candidate, ...]) -> None:
        self.root = root
        self.candidates = candidates
        self.cache: dict[tuple[tuple[float, ...], str], dict[str, object]] = {}
        self.fingerprints: dict[tuple[float, ...], str] = {}
        self.bundle_counter = 0
        self.state_evaluations = 0

    def _key(self, values: dict[str, float]) -> tuple[float, ...]:
        return tuple(float(values[candidate.identity]) for candidate in self.candidates)

    def evaluate(
        self,
        values: dict[str, float],
        rows: list[dict[str, object]],
    ) -> list[dict[str, object]]:
        key = self._key(values)
        missing = [row for row in rows if (key, str(row["observation_id"])) not in self.cache]
        if missing:
            import epcsaft

            bundle = self.root / f"evaluation-{self.bundle_counter:03d}"
            self.bundle_counter += 1
            _prepare_bundle(bundle, self.candidates, values)
            model = epcsaft.Mixture(
                epcsaft.Parameters.from_bundle(bundle, components=comparison.COMPONENT_IDS)
            )
            self.fingerprints[key] = model.parameter_fingerprint
            for row in missing:
                evaluator = comparison.ProviderActivityEvaluator(
                    model,
                    epcsaft.unit_registry,
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
                state = evaluator._liquid_state(solved.mole_fractions)
                if state.fugacity is None:
                    raise ValueError("Provider did not return component fugacities")
                prediction = float(state.fugacity.value[0].to("pascal").magnitude)
                result = {
                    **row,
                    "predicted_pco2_pa": prediction,
                    "log10_residual": math.log10(prediction / float(row["observed_pco2_pa"])),
                    "success": solved.success,
                    "max_abs_reaction_balance_residual": solved.max_abs_residual,
                    "parameter_fingerprint": model.parameter_fingerprint,
                }
                if not solved.success:
                    raise RuntimeError(f"reactive state failed: {row['observation_id']}")
                self.cache[(key, str(row["observation_id"]))] = result
                self.state_evaluations += 1
            shutil.rmtree(bundle)
        return [self.cache[(key, str(row["observation_id"]))] for row in rows]


def _residuals(results: list[dict[str, object]]) -> np.ndarray:
    return np.asarray([float(row["log10_residual"]) for row in results])


def _screen(
    experiment: Experiment,
    candidates: tuple[Candidate, ...],
    starts: dict[str, float],
    training: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str], dict[str, Any]]:
    ordered = sorted(training, key=lambda row: float(row["loading"]))
    indices = np.rint(np.linspace(0, len(ordered) - 1, SCREEN_ROW_COUNT)).astype(int)
    screen_rows = [ordered[int(index)] for index in indices]
    base = _residuals(experiment.evaluate(starts, screen_rows))
    by_row: list[dict[str, object]] = []
    columns: dict[str, np.ndarray] = {}
    summary: list[dict[str, object]] = []
    for candidate in candidates:
        perturbed = dict(starts)
        perturbed[candidate.identity] += candidate.screen_step
        changed = perturbed[candidate.identity] - starts[candidate.identity]
        derivative = (
            _residuals(experiment.evaluate(perturbed, screen_rows)) - base
        ) / (changed / candidate.screen_scale)
        response = derivative * (changed / candidate.screen_scale)
        columns[candidate.identity] = derivative
        for row, value in zip(screen_rows, derivative, strict=True):
            by_row.append(
                {
                    "parameter_identity": candidate.identity,
                    "observation_id": row["observation_id"],
                    "mea_mass_fraction": row["mea_mass_fraction"],
                    "loading": row["loading"],
                    "standardized_log10_sensitivity": float(value),
                }
            )
        summary.append(
            {
                "parameter_identity": candidate.identity,
                "label": candidate.label,
                "fit_eligible": str(candidate.fit_eligible).lower(),
                "start": candidate.start,
                "unit": candidate.unit,
                "screen_step": candidate.screen_step,
                "screen_scale": candidate.screen_scale,
                "rms_standardized_log10_sensitivity": float(np.sqrt(np.mean(derivative**2))),
                "max_abs_standardized_log10_sensitivity": float(np.max(np.abs(derivative))),
                "mean_standardized_log10_sensitivity": float(np.mean(derivative)),
                "rms_log10_response_for_screen_step": float(
                    np.sqrt(np.mean(response**2))
                ),
            }
        )

    active = [candidate for candidate in candidates if candidate.fit_eligible]
    matrix = np.column_stack([columns[candidate.identity] for candidate in active])
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = max(1, int(np.sum(singular >= IDENTIFIABLE_SINGULAR_RATIO * singular[0])))
    _, _, pivots = qr(matrix, mode="economic", pivoting=True)
    selected = [active[int(index)].identity for index in pivots[:rank]]
    correlations: list[dict[str, object]] = []
    correlation = np.corrcoef(matrix, rowvar=False)
    for row_index, left in enumerate(active):
        for column_index, right in enumerate(active):
            correlations.append(
                {
                    "parameter_a": left.identity,
                    "parameter_b": right.identity,
                    "correlation": float(correlation[row_index, column_index]),
                }
            )
    diagnostics = {
        "screen_observation_ids": [str(row["observation_id"]) for row in screen_rows],
        "active_singular_values": singular.tolist(),
        "rank_threshold_ratio": IDENTIFIABLE_SINGULAR_RATIO,
        "active_screen_rank": rank,
        "selected_fit_parameters": selected,
    }
    return summary, by_row, selected, {**diagnostics, "correlations": correlations}


def _jacobian(
    experiment: Experiment,
    candidates: dict[str, Candidate],
    selected: list[str],
    values: dict[str, float],
    rows: list[dict[str, object]],
) -> tuple[np.ndarray, np.ndarray]:
    base = _residuals(experiment.evaluate(values, rows))
    columns = []
    for identity in selected:
        candidate = candidates[identity]
        assert candidate.affine_scale is not None
        assert candidate.lower is not None and candidate.upper is not None
        step = FIT_DIFFERENCE_AFFINE * candidate.affine_scale
        if values[identity] + step > candidate.upper:
            step = -step
        perturbed = dict(values)
        perturbed[identity] += step
        changed_affine = step / candidate.affine_scale
        columns.append(
            (_residuals(experiment.evaluate(perturbed, rows)) - base) / changed_affine
        )
    return base, np.column_stack(columns)


def _fit(
    experiment: Experiment,
    candidates: tuple[Candidate, ...],
    selected: list[str],
    starts: dict[str, float],
    training: list[dict[str, object]],
) -> tuple[dict[str, float], dict[str, Any]]:
    lookup = {candidate.identity: candidate for candidate in candidates}
    values = dict(starts)
    iterations: list[dict[str, object]] = []
    for iteration in range(FIT_MAX_ITERATIONS):
        residual, jacobian = _jacobian(experiment, lookup, selected, values, training)
        lower = []
        upper = []
        for identity in selected:
            candidate = lookup[identity]
            assert candidate.affine_scale is not None
            assert candidate.lower is not None and candidate.upper is not None
            lower.append(
                max(
                    (candidate.lower - values[identity]) / candidate.affine_scale,
                    -FIT_TRUST_RADIUS_AFFINE,
                )
            )
            upper.append(
                min(
                    (candidate.upper - values[identity]) / candidate.affine_scale,
                    FIT_TRUST_RADIUS_AFFINE,
                )
            )
        step = lsq_linear(
            jacobian,
            -residual,
            bounds=(np.asarray(lower), np.asarray(upper)),
            lsmr_tol="auto",
        ).x
        base_rmse = float(np.sqrt(np.mean(residual**2)))
        accepted = False
        accepted_factor = 0.0
        candidate_rmse = base_rmse
        candidate_values = dict(values)
        for factor in (1.0, 0.5, 0.25, 0.125):
            trial = dict(values)
            for identity, delta in zip(selected, step, strict=True):
                scale = float(lookup[identity].affine_scale)
                trial[identity] += factor * float(delta) * scale
            trial_residual = _residuals(experiment.evaluate(trial, training))
            trial_rmse = float(np.sqrt(np.mean(trial_residual**2)))
            if trial_rmse < candidate_rmse:
                accepted = True
                accepted_factor = factor
                candidate_rmse = trial_rmse
                candidate_values = trial
                break
        iterations.append(
            {
                "iteration": iteration + 1,
                "base_log10_rmse": base_rmse,
                "candidate_log10_rmse": candidate_rmse,
                "accepted": accepted,
                "line_search_factor": accepted_factor,
                "affine_step": {identity: float(delta) for identity, delta in zip(selected, step, strict=True)},
            }
        )
        if not accepted:
            break
        values = candidate_values
        if base_rmse - candidate_rmse < 1.0e-4:
            break

    final_residual, final_jacobian = _jacobian(
        experiment, lookup, selected, values, training
    )
    singular = np.linalg.svd(final_jacobian, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    scaled_bound_distances = {}
    for identity in selected:
        candidate = lookup[identity]
        assert candidate.affine_scale is not None
        assert candidate.lower is not None and candidate.upper is not None
        scaled_bound_distances[identity] = min(
            values[identity] - candidate.lower,
            candidate.upper - values[identity],
        ) / candidate.affine_scale
    diagnostics = {
        "method": "two-step bounded Gauss-Newton with forward finite differences and backtracking",
        "residual": "log10(predicted_pCO2_pa / observed_pCO2_pa)",
        "weights": "equal row weights; source-reported numeric covariance is unavailable",
        "regularization": "none",
        "selected_parameters": selected,
        "iterations": iterations,
        "final_values": {identity: values[identity] for identity in selected},
        "final_log10_rmse": float(np.sqrt(np.mean(final_residual**2))),
        "final_jacobian_singular_values": singular.tolist(),
        "final_jacobian_rank": int(np.linalg.matrix_rank(final_jacobian)),
        "final_jacobian_condition_number": condition,
        "scaled_distance_to_nearest_bound": scaled_bound_distances,
        "numerically_supported": bool(
            np.linalg.matrix_rank(final_jacobian) == len(selected)
            and condition <= FIT_CONDITION_LIMIT
            and min(scaled_bound_distances.values()) > 1.0e-7
        ),
        "confirmation_multistart": "not_run",
        "promotion_allowed": False,
    }
    return values, diagnostics


def _metrics(
    predictions: list[dict[str, object]],
    model_id: str,
    selection: str,
) -> dict[str, object]:
    chosen = [
        row
        for row in predictions
        if row["model_id"] == model_id
        and (selection == "all" or row["split"] == selection)
    ]
    residual = np.asarray([float(row["log10_residual"]) for row in chosen])
    ratio = np.asarray(
        [float(row["predicted_pco2_pa"]) / float(row["observed_pco2_pa"]) for row in chosen]
    )
    return {
        "model_id": model_id,
        "selection": selection,
        "row_count": len(chosen),
        "log10_rmse": float(np.sqrt(np.mean(residual**2))),
        "log10_bias": float(np.mean(residual)),
        "aard_percent": float(100.0 * np.mean(np.abs(ratio - 1.0))),
        "max_abs_log10_residual": float(np.max(np.abs(residual))),
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    inventory, executable = _temperature_inventory()
    training = [row for row in executable if row["split"] == "training"]
    candidates = _candidates()
    starts = {candidate.identity: candidate.start for candidate in candidates}
    with tempfile.TemporaryDirectory(prefix="mea-pressure-fit-") as temporary:
        experiment = Experiment(Path(temporary), candidates)
        sensitivity, sensitivity_by_row, selected, screen = _screen(
            experiment, candidates, starts, training
        )
        fitted, fit = _fit(experiment, candidates, selected, starts, training)
        predictions: list[dict[str, object]] = []
        for model_id, values in (("M2", starts), ("M4", fitted)):
            for row in experiment.evaluate(values, executable):
                predictions.append({"model_id": model_id, **row})
        state_evaluations = experiment.state_evaluations
        fingerprints = {
            model_id: experiment.fingerprints[experiment._key(values)]
            for model_id, values in (("M2", starts), ("M4", fitted))
        }

    fit_parameters = []
    for candidate in candidates:
        fit_parameters.append(
            {
                "parameter_identity": candidate.identity,
                "label": candidate.label,
                "fit_eligible": str(candidate.fit_eligible).lower(),
                "selected_for_fit": str(candidate.identity in selected).lower(),
                "unit": candidate.unit,
                "start": candidate.start,
                "fitted_value": fitted[candidate.identity],
                "change": fitted[candidate.identity] - candidate.start,
                "lower_bound": "" if candidate.lower is None else candidate.lower,
                "upper_bound": "" if candidate.upper is None else candidate.upper,
                "affine_scale": "" if candidate.affine_scale is None else candidate.affine_scale,
            }
        )
    correlation_rows = screen.pop("correlations")
    prediction_rows = []
    for row in predictions:
        prediction_rows.append(
            {
                "model_id": row["model_id"],
                "observation_id": row["observation_id"],
                "source_key": row["source_key"],
                "mea_mass_fraction": row["mea_mass_fraction"],
                "loading_mol_co2_per_mol_mea": row["loading"],
                "temperature_k": row["temperature_k"],
                "state_pressure_pa": row["pressure_pa"],
                "observed_pco2_pa": row["observed_pco2_pa"],
                "predicted_pco2_pa": row["predicted_pco2_pa"],
                "log10_residual": row["log10_residual"],
                "split": row["split"],
                "role": row["role"],
                "group_id": row["group_id"],
                "measurement_origin": row["measurement_origin"],
                "success": row["success"],
                "max_abs_reaction_balance_residual": row[
                    "max_abs_reaction_balance_residual"
                ],
                "parameter_fingerprint": row["parameter_fingerprint"],
            }
        )
    metrics = [
        _metrics(predictions, model_id, split)
        for model_id in ("M2", "M4")
        for split in ("training", "validation", "all")
    ]
    for model_id in ("M2", "M4"):
        for source, fraction in (("Hilliard2008", 0.17), ("Hilliard2008", 0.40), ("Jou1995", 0.30)):
            group = [
                row
                for row in predictions
                if row["model_id"] == model_id
                and row["source_key"] == source
                and math.isclose(float(row["mea_mass_fraction"]), fraction)
            ]
            residual = np.asarray([float(row["log10_residual"]) for row in group])
            ratio = np.asarray(
                [float(row["predicted_pco2_pa"]) / float(row["observed_pco2_pa"]) for row in group]
            )
            metrics.append(
                {
                    "model_id": model_id,
                    "selection": f"validation:{source}:w={fraction:g}",
                    "row_count": len(group),
                    "log10_rmse": float(np.sqrt(np.mean(residual**2))),
                    "log10_bias": float(np.mean(residual)),
                    "aard_percent": float(100.0 * np.mean(np.abs(ratio - 1.0))),
                    "max_abs_log10_residual": float(np.max(np.abs(residual))),
                }
            )

    output_tables = {
        "pressure_row_inventory.csv": inventory,
        "pressure_parameter_sensitivity.csv": sensitivity,
        "pressure_sensitivity_by_row.csv": sensitivity_by_row,
        "pressure_active_correlation.csv": correlation_rows,
        "pressure_fit_parameters.csv": fit_parameters,
        "pressure_fit_predictions.csv": prediction_rows,
        "pressure_fit_metrics.csv": metrics,
    }
    for name, rows in output_tables.items():
        _write_csv(RESULTS / name, rows)
    receipt = {
        "analysis": "313.15 K pressure sensitivity and M4 exploratory fit",
        "status": "nonpromoting_experiment",
        "temperature_selection": {
            "temperature_k": TEMPERATURE_K,
            "selection_basis": "largest executable pCO2 pool in the canonical repository data",
            "executable_rows": 44,
            "training_rows": 24,
            "reserved_validation_rows": 20,
        },
        "model_definition": {
            "M2": "fixed Pabsch CO2-water induced association with preregistered starts",
            "M4": "M2 plus the sensitivity-selected subset of the three preregistered ionic coordinates",
        },
        "screen": screen,
        "fit": fit,
        "provider": {
            "version": metadata.version("epcsaft"),
            **comparison._installed_provider_identity(),
            "build_source": comparison._git_identity(PROVIDER_SOURCE),
        },
        "parameter_fingerprints": fingerprints,
        "state_evaluation_count": state_evaluations,
        "source_hashes": {
            str(path.relative_to(ROOT)): comparison._sha256(path)
            for path in (
                comparison.VLE_OBSERVATIONS,
                comparison.PCO2_METROLOGY,
                SPLIT_MANIFEST,
                comparison.PREREGISTRATION,
            )
        },
        "outputs": {
            name: comparison._sha256(RESULTS / name) for name in output_tables
        },
        "claim_boundary": (
            "Pressure-only, single-temperature sensitivity experiment. Validation rows were not "
            "used for selection or fitting; no parameter promotion or manuscript claim is allowed."
        ),
    }
    (RESULTS / "pressure_fit_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
