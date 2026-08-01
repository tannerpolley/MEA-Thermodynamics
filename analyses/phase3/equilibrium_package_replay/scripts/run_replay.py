from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
COMPARISON_SCRIPTS = ROOT / "analyses/phase3/m0_m3_model_comparison/scripts"
if str(COMPARISON_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(COMPARISON_SCRIPTS))

import run_polar_comparison as variants  # noqa: E402


ANALYSIS = ROOT / "analyses/phase3/equilibrium_package_replay"
RESULTS = ANALYSIS / "results"
REACTIONS = ROOT / "data/reference/MEA/manifests/chemical_reaction_source_contract.json"
SENTINEL = (
    ROOT / "data/reference/MEA/manifests/homogeneous_speciation_sentinel_contract.json"
)
LEGACY = (
    ROOT
    / "analyses/phase3/m0_m3_model_comparison/results/m0_m5_pressure_predictions.csv"
)
TARGET_OBSERVATION = "vle_obs_0137"
MODEL_ORDER = variants.MODEL_ORDER


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_identity(path: Path, ref: str) -> dict[str, str]:
    return {
        "commit": subprocess.check_output(
            ["git", "rev-parse", ref], cwd=path, text=True
        ).strip(),
        "tree": subprocess.check_output(
            ["git", "rev-parse", f"{ref}^{{tree}}"], cwd=path, text=True
        ).strip(),
    }


def _installed_identity(distribution_name: str, header: str) -> dict[str, str]:
    distribution = metadata.distribution(distribution_name)
    record = distribution.read_text("RECORD")
    header_path = distribution.locate_file(header)
    if record is None or not header_path.is_file():
        raise ValueError(f"installed {distribution_name} identity is incomplete")
    return {
        "distribution": distribution_name,
        "version": distribution.version,
        "record_sha256": hashlib.sha256(record.encode()).hexdigest(),
        "public_header_sha256": _sha256(header_path),
    }


def _reaction_ln_k(reaction: dict[str, Any], temperature_k: float) -> float:
    correlation = reaction["correlation"]
    kind = correlation["kind"]
    if kind == "ln_a_plus_b_over_t_plus_c_ln_t_plus_d_t":
        return (
            correlation["a"]
            + correlation["b_k"] / temperature_k
            + correlation["c"] * math.log(temperature_k)
            + correlation["d_per_k"] * temperature_k
        )
    if kind == "ln_a_plus_b_over_t":
        return correlation["a"] + correlation["b_k"] / temperature_k
    if kind == "ln_from_negative_log10_a_over_t_plus_b_plus_c_t":
        return -math.log(10.0) * (
            correlation["a_k"] / temperature_k
            + correlation["b"]
            + correlation["c_per_k"] * temperature_k
        )
    raise ValueError(f"unsupported reaction correlation: {kind}")


def _legacy_rows() -> dict[str, dict[str, str]]:
    with LEGACY.open(newline="", encoding="utf-8") as handle:
        rows = {
            row["model_id"]: row
            for row in csv.DictReader(handle)
            if row["observation_id"] == TARGET_OBSERVATION
        }
    if tuple(model for model in MODEL_ORDER if model in rows) != MODEL_ORDER:
        raise ValueError("the frozen tracer does not contain all eight model variants")
    return rows


def _problem(
    equilibrium: Any,
    bundle: Path,
    contract: dict[str, Any],
    row: dict[str, str],
) -> Any:
    component_ids = tuple(contract["provider_species_order"])
    species = contract["species"]
    element_order = tuple(contract["balance_row_order"])
    elemental = tuple(
        tuple(float(item["formula"][element]) for item in species)
        for element in element_order
    )
    balances = (
        elemental[element_order.index("C")],
        elemental[element_order.index("N")],
    )
    mass_fraction = float(row["mea_mass_fraction"])
    loading = float(row["loading_mol_co2_per_mol_mea"])
    sentinel = json.loads(SENTINEL.read_text(encoding="utf-8"))
    masses = sentinel["molar_mass_basis"]["values"]
    water = (
        (1.0 - mass_fraction)
        / mass_fraction
        * float(masses["MEA"])
        / float(masses["H2O"])
    )
    feed = (loading, 1.0, water, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    totals = tuple(
        math.fsum(balance[index] * feed[index] for index in range(len(feed)))
        for balance in balances
    )

    with (bundle / "single.csv").open(newline="", encoding="utf-8") as handle:
        reported_by_component = {
            item["component_id"]: float(item["value"])
            for item in csv.DictReader(handle)
            if item["family"] == "molar_mass"
        }
    reported = np.asarray(
        [reported_by_component[component_id] for component_id in component_ids]
    )
    elemental_array = np.asarray(elemental)
    elemental_masses, *_ = np.linalg.lstsq(elemental_array.T, reported, rcond=None)
    certified_masses = tuple(
        float(value) for value in elemental_array.T @ elemental_masses
    )

    temperature_k = float(row["temperature_k"])
    common = contract["common_source_standard_state"]
    ln_k = tuple(
        _reaction_ln_k(reaction, temperature_k) + float(offset)
        for reaction, offset in zip(
            contract["reactions"],
            common["source_to_common_ln_k_offsets"],
            strict=True,
        )
    )
    reference_pressure = float(
        contract["provider_transform"]["deterministic_payload"][
            "source_standard_reference_pressure_pa"
        ]
    )
    standard_state = equilibrium.ChemicalStandardState(
        id=common["identity"],
        activity_scale_id=common["identity"],
        log_activity_scale_factors=tuple(
            common["log_activity_scale_factors_by_species"]
        ),
        reference_pressure_pa=reference_pressure,
    )
    return equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=tuple(int(item["charge"]) for item in species),
        molar_masses_kg_per_mol=certified_masses,
        balance_matrix=balances,
        conserved_totals=totals,
        reaction_matrix=tuple(
            tuple(float(value) for value in reaction["stoichiometry"])
            for reaction in contract["reactions"]
        ),
        feed_amounts_mol=feed,
        equilibrium_constants=tuple(
            equilibrium.ChemicalEquilibriumConstant(
                ln_value=value,
                source_id="+".join(reaction["source_record_ids"]),
                reference_id=common["identity"],
                reaction_orientation="products_positive",
                conversion_id="source-standard-state-to-provider-neutral-reference",
                dimensionless=True,
            )
            for reaction, value in zip(contract["reactions"], ln_k, strict=True)
        ),
        strict_interior_amount_floor_mol=1.0e-18,
        source_standard_state=standard_state,
    )


def _failure_class(reason: str) -> str:
    if "neutral-reference ABI contract is incomplete" in reason:
        return "BLOCKED_PROVIDER_NEUTRAL_REFERENCE_ABI"
    if "inverse-packing geometry ABI contract is incomplete" in reason:
        return "BLOCKED_PROVIDER_INVERSE_PACKING_ABI"
    if "primal_solution_not_certified" in reason:
        return "NOT_CERTIFIED_LOCAL_MINIMUM"
    return "BLOCKED_UNCLASSIFIED"


def _evaluate_model(
    epcsaft: Any,
    equilibrium: Any,
    model_id: str,
    row: dict[str, str],
    contract: dict[str, Any],
    work: Path,
) -> dict[str, Any]:
    bundle = work / model_id.lower()
    variants._prepare_bundle(bundle, model_id)
    component_ids = tuple(contract["provider_species_order"])
    model = epcsaft.Mixture(
        epcsaft.Parameters.from_bundle(bundle, components=component_ids)
    )
    result: dict[str, Any] = {
        "model_id": model_id,
        "observation_id": TARGET_OBSERVATION,
        "temperature_k": float(row["temperature_k"]),
        "state_pressure_pa": float(row["state_pressure_pa"]),
        "loading_mol_co2_per_mol_mea": float(row["loading_mol_co2_per_mol_mea"]),
        "legacy_predicted_pco2_pa": float(row["predicted_pco2_pa"]),
        "parameter_fingerprint": model.parameter_fingerprint,
        "status": "NOT_RUN",
        "reason": "",
        "equilibrium_predicted_pco2_pa": "",
        "log10_equilibrium_over_legacy": "",
        "solver_status": "",
        "numerical_status": "",
        "physical_status": "",
        "local_minimum_status": "",
        "elapsed_s": 0.0,
    }
    problem = _problem(equilibrium, bundle, contract, row)
    phase = equilibrium.ProviderPhase(
        model=model,
        expected_parameter_fingerprint=model.parameter_fingerprint,
        admissible_packing_fraction_interval=(1.0e-6, 0.74),
    )
    started = time.perf_counter()
    try:
        solved = equilibrium.chemical_equilibrium(
            phase,
            result["temperature_k"] * epcsaft.unit_registry.kelvin,
            result["state_pressure_pa"] * epcsaft.unit_registry.pascal,
            problem,
        )
        molar_density = math.fsum(solved.amounts_mol) / solved.volume_m3
        state = model.state(
            T=result["temperature_k"] * epcsaft.unit_registry.kelvin,
            rho=molar_density
            * epcsaft.unit_registry.mole
            / epcsaft.unit_registry.meter**3,
            x=solved.mole_fractions,
        )
        if state.fugacity is None:
            raise ValueError("Provider returned no fugacity for a certified state")
        pressure = float(state.fugacity.value[0].to("pascal").magnitude)
        result.update(
            status="CERTIFIED_LOCAL_EQUILIBRIUM",
            equilibrium_predicted_pco2_pa=pressure,
            log10_equilibrium_over_legacy=math.log10(
                pressure / result["legacy_predicted_pco2_pa"]
            ),
            solver_status=solved.diagnostics.solver_status,
            numerical_status=solved.diagnostics.numerical_status,
            physical_status=solved.diagnostics.physical_status,
            local_minimum_status=solved.diagnostics.local_minimum_status,
            equilibrium_mole_fractions=list(solved.mole_fractions),
            ln_k_provider_basis=list(solved.ln_k_provider_basis or ()),
            equilibrium_artifact_identity=asdict(solved.artifact_identity),
        )
    except equilibrium.ChemicalEquilibriumError as error:
        diagnostics = error.diagnostics
        result.update(
            status=_failure_class(str(error)),
            reason=str(error),
            solver_status=diagnostics.solver_status,
            numerical_status=diagnostics.numerical_status,
            physical_status=diagnostics.physical_status,
            local_minimum_status=diagnostics.local_minimum_status,
            chemical_certification_level=diagnostics.chemical_certification_level,
            provider_domain_status=diagnostics.provider_domain_status,
        )
    result["elapsed_s"] = time.perf_counter() - started
    return result


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "model_id",
        "observation_id",
        "temperature_k",
        "state_pressure_pa",
        "loading_mol_co2_per_mol_mea",
        "legacy_predicted_pco2_pa",
        "equilibrium_predicted_pco2_pa",
        "log10_equilibrium_over_legacy",
        "status",
        "reason",
        "solver_status",
        "numerical_status",
        "physical_status",
        "local_minimum_status",
        "elapsed_s",
        "parameter_fingerprint",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields} for row in rows
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--models",
        nargs="+",
        choices=MODEL_ORDER,
        default=list(MODEL_ORDER),
    )
    parser.add_argument(
        "--provider-source",
        type=Path,
        default=ROOT.parent / "ePC-SAFT-project/ePC-SAFT-eos",
    )
    parser.add_argument(
        "--equilibrium-source",
        type=Path,
        default=ROOT.parent / "ePC-SAFT-project/ePC-SAFT-equilibrium",
    )
    parser.add_argument("--provider-ref", default="main")
    parser.add_argument("--equilibrium-ref", default="main")
    parser.add_argument("--provider-wheel-sha256", default="not_recorded")
    parser.add_argument("--equilibrium-wheel-sha256", default="not_recorded")
    args = parser.parse_args()

    import epcsaft
    import epcsaft_equilibrium

    contract = json.loads(REACTIONS.read_text(encoding="utf-8"))
    legacy = _legacy_rows()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="mea-equilibrium-replay-") as temporary:
        work = Path(temporary)
        for model_id in args.models:
            row = _evaluate_model(
                epcsaft,
                epcsaft_equilibrium,
                model_id,
                legacy[model_id],
                contract,
                work,
            )
            rows.append(row)
            print(f"{model_id}: {row['status']} ({row['elapsed_s']:.3f} s)", flush=True)

    statuses = {str(row["status"]) for row in rows}
    if "BLOCKED_UNCLASSIFIED" in statuses:
        overall = "BLOCKED_UNCLASSIFIED"
    elif statuses == {"CERTIFIED_LOCAL_EQUILIBRIUM"}:
        overall = "PARITY_EVALUATED"
    elif any(status.startswith("BLOCKED_") for status in statuses):
        overall = "BLOCKED_SHARED_PROVIDER_CAPABILITY"
    else:
        overall = "NOT_CERTIFIED"
    _write_csv(RESULTS / "equilibrium_replay_comparison.csv", rows)
    receipt = {
        "analysis": "public ePC-SAFT Equilibrium replay",
        "overall_status": overall,
        "target_observation_id": TARGET_OBSERVATION,
        "model_order": list(args.models),
        "provider_source": _git_identity(args.provider_source, args.provider_ref),
        "equilibrium_source": _git_identity(
            args.equilibrium_source, args.equilibrium_ref
        ),
        "installed_provider": _installed_identity(
            "epcsaft", "epcsaft/include/epcsaft/native_sdk_v1.h"
        ),
        "installed_equilibrium": _installed_identity(
            "epcsaft-equilibrium",
            "epcsaft_equilibrium/include/epcsaft/regression/evaluator_v1.h",
        ),
        "wheel_sha256": {
            "epcsaft": args.provider_wheel_sha256,
            "epcsaft-equilibrium": args.equilibrium_wheel_sha256,
        },
        "source_hashes": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (REACTIONS, SENTINEL, LEGACY)
        },
        "results": rows,
        "claim_boundary": (
            "A successful row establishes one certified local fixed-T,P homogeneous "
            "state only. It does not establish global equilibrium, parameter "
            "identifiability, predictive validity, or regression readiness."
        ),
    }
    status_path = RESULTS / "equilibrium_replay_status.json"
    status_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt["output_hashes"] = {
        "equilibrium_replay_comparison.csv": _sha256(
            RESULTS / "equilibrium_replay_comparison.csv"
        )
    }
    status_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if overall == "BLOCKED_UNCLASSIFIED":
        raise RuntimeError("Equilibrium replay encountered an unclassified failure")


if __name__ == "__main__":
    main()
