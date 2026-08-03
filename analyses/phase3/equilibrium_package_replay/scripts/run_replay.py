from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
import zipfile
from dataclasses import asdict
from dataclasses import fields as dataclass_fields
from importlib import metadata
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[4]
COMPARISON_SCRIPTS = ROOT / "analyses/phase3/m0_m3_model_comparison/scripts"
ANALYSIS = ROOT / "analyses/phase3/equilibrium_package_replay"
RESULTS = ANALYSIS / "results"
REACTIONS = ROOT / "data/reference/MEA/manifests/chemical_reaction_source_contract.json"
SOURCE_REFERENCE = ANALYSIS / "source_reference_transfer_contract.json"
SENTINEL = (
    ROOT / "data/reference/MEA/manifests/homogeneous_speciation_sentinel_contract.json"
)
LEGACY = (
    ROOT
    / "analyses/phase3/m0_m3_model_comparison/results/m0_m5_pressure_predictions.csv"
)
TARGET_OBSERVATION = "vle_obs_0137"
MODEL_ORDER = ("M0", "M1", "M2", "M3", "M4A", "M4B", "M5Q", "M5")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_identity(path: Path, ref: str) -> dict[str, str]:
    commit = subprocess.check_output(
        ["git", "rev-parse", ref], cwd=path, text=True
    ).strip()
    remote_branches = subprocess.check_output(
        ["git", "branch", "-r", "--contains", commit], cwd=path, text=True
    ).splitlines()
    return {
        "commit": commit,
        "tree": subprocess.check_output(
            ["git", "rev-parse", f"{ref}^{{tree}}"], cwd=path, text=True
        ).strip(),
        "origin_contains_commit": any(
            branch.strip().startswith("origin/") for branch in remote_branches
        ),
    }


def _mapping_sha256(values: dict[str, str]) -> str:
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def _installed_identity(
    distribution_name: str, package: str, header: str
) -> dict[str, str]:
    distribution = metadata.distribution(distribution_name)
    record = distribution.read_text("RECORD")
    header_path = distribution.locate_file(header)
    package_hashes = {
        str(item): _sha256(distribution.locate_file(item))
        for item in distribution.files or ()
        if str(item).startswith(f"{package}/")
        and "__pycache__" not in str(item)
        and distribution.locate_file(item).is_file()
    }
    if record is None or not header_path.is_file() or not package_hashes:
        raise ValueError(f"installed {distribution_name} identity is incomplete")
    return {
        "distribution": distribution_name,
        "version": distribution.version,
        "record_sha256": hashlib.sha256(record.encode()).hexdigest(),
        "public_header_sha256": _sha256(header_path),
        "package_payload_sha256": _mapping_sha256(package_hashes),
        "installation_root": str(distribution.locate_file("")),
    }


def _wheel_identity(path: Path, package: str, header: str) -> dict[str, str]:
    if not path.is_file():
        raise ValueError(f"wheel does not exist: {path}")
    with zipfile.ZipFile(path) as archive:
        record_names = [
            name for name in archive.namelist() if name.endswith(".dist-info/RECORD")
        ]
        header_names = [name for name in archive.namelist() if name.endswith(header)]
        if len(record_names) != 1 or len(header_names) != 1:
            raise ValueError(f"wheel identity is incomplete: {path}")
        record = archive.read(record_names[0])
        public_header = archive.read(header_names[0])
        package_hashes = {
            name: hashlib.sha256(archive.read(name)).hexdigest()
            for name in archive.namelist()
            if name.startswith(f"{package}/") and not name.endswith("/")
        }
    return {
        "filename": path.name,
        "sha256": _sha256(path),
        "record_sha256": hashlib.sha256(record).hexdigest(),
        "public_header_sha256": hashlib.sha256(public_header).hexdigest(),
        "package_payload_sha256": _mapping_sha256(package_hashes),
    }


def _verify_installed_artifact(
    installed: dict[str, str],
    wheel: dict[str, str],
    expected_sha256: str,
) -> None:
    if wheel["sha256"] != expected_sha256:
        raise ValueError("wheel SHA-256 does not match the required artifact")
    for field in ("package_payload_sha256", "public_header_sha256"):
        if installed[field] != wheel[field]:
            raise ValueError(f"installed artifact does not match wheel {field}")


def _source_reference_capability(
    epcsaft: Any,
    equilibrium: Any,
    reaction_contract: dict[str, Any],
    transfer_contract: dict[str, Any],
) -> dict[str, Any]:
    source = transfer_contract["source_reference"]
    required = transfer_contract["required_public_capability"]
    species_order = [item["name"] for item in reaction_contract["species"]]
    if source["species_order"] != species_order:
        raise ValueError(
            "source-reference species order differs from reaction contract"
        )
    if (
        source["provider_component_order"]
        != reaction_contract["provider_species_order"]
    ):
        raise ValueError(
            "source-reference Provider order differs from reaction contract"
        )
    composition = tuple(float(value) for value in source["mole_fractions"])
    water_index = species_order.index("H2O")
    if (
        len(composition) != len(species_order)
        or not math.isclose(math.fsum(composition), 1.0, rel_tol=0.0, abs_tol=0.0)
        or composition[water_index] != 1.0
        or any(
            value != 0.0
            for index, value in enumerate(composition)
            if index != water_index
        )
    ):
        raise ValueError("source-reference composition is not exact pure water")
    standard_state_fields = tuple(
        field.name for field in dataclass_fields(equilibrium.ChemicalStandardState)
    )
    required_field = str(required["equilibrium_standard_state_field"])
    operation_name = str(required["public_operation_name"])
    operations = {
        "epcsaft": callable(getattr(epcsaft, operation_name, None)),
        "epcsaft.Mixture": callable(getattr(epcsaft.Mixture, operation_name, None)),
        "epcsaft_equilibrium": callable(getattr(equilibrium, operation_name, None)),
    }
    supported = required_field in standard_state_fields and any(operations.values())
    return {
        "required_capability_identity": required["identity"],
        "status": "SUPPORTED" if supported else "UNSUPPORTED",
        "source_reference_contract_identity": transfer_contract["identity"],
        "source_reference_kind": source["kind"],
        "source_reference_composition": list(composition),
        "source_reference_species_order": species_order,
        "equilibrium_standard_state_fields": list(standard_state_fields),
        "required_equilibrium_standard_state_field": required_field,
        "public_operation_name": operation_name,
        "public_operation_availability": operations,
        "reason": (
            "installed public EOS/Equilibrium APIs expose no caller-declared "
            "source-solvent reference transfer"
            if not supported
            else "installed public APIs advertise the required transfer"
        ),
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


def _blocked_source_reference_row(
    model_id: str,
    row: dict[str, str],
    capability: dict[str, Any],
    status: str,
) -> dict[str, Any]:
    return {
        "model_id": model_id,
        "observation_id": TARGET_OBSERVATION,
        "temperature_k": float(row["temperature_k"]),
        "state_pressure_pa": float(row["state_pressure_pa"]),
        "loading_mol_co2_per_mol_mea": float(row["loading_mol_co2_per_mol_mea"]),
        "legacy_predicted_pco2_pa": float(row["predicted_pco2_pa"]),
        "parameter_fingerprint": "",
        "status": status,
        "reason": capability["reason"],
        "equilibrium_predicted_pco2_pa": "",
        "log10_equilibrium_over_legacy": "",
        "solver_status": "NOT_RUN",
        "failure_kind": "unsupported_source_reference_transfer",
        "chemical_certification_level": "NOT_EVALUATED",
        "numerical_status": "NOT_EVALUATED",
        "physical_status": "NOT_EVALUATED",
        "local_minimum_status": "NOT_EVALUATED",
        "first_failed_numerical_criterion": "",
        "first_failed_physical_criterion": "",
        "elapsed_s": 0.0,
        "reaction_affinity_inf_norm": None,
        "pressure_relative_residual": None,
        "kkt_stationarity_inf_norm": None,
        "reduced_hessian_inertia": None,
        "minimum_amount_mol": None,
        "packing_fraction": None,
        "start_accounting": {
            "generated": 0,
            "evaluated": 0,
            "dropped": 0,
            "failed": 0,
            "certified": 0,
        },
    }


def _problem(
    equilibrium: Any,
    bundle: Path,
    contract: dict[str, Any],
    row: dict[str, str],
) -> Any:
    import numpy as np

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


def _failure_class(failure_kind: str) -> str:
    return {
        "physical_domain_failure": "NOT_ADMITTED_PHYSICAL_DOMAIN_FAILURE",
        "exhausted_multistart_search": ("NOT_CERTIFIED_EXHAUSTED_MULTISTART_SEARCH"),
    }.get(failure_kind, "BLOCKED_UNCLASSIFIED")


def _diagnostic_receipt(diagnostics: Any) -> dict[str, Any]:
    receipt = asdict(diagnostics)
    for field in (
        "reduced_hessian",
        "reduced_hessian_nullspace_basis",
        "objective_gradient",
        "constraint_values",
        "constraint_jacobian",
        "lagrangian_gradient",
        "equality_multipliers",
        "lagrangian_hessian",
        "covariant_lagrangian_hessian",
        "kkt_root_jacobian",
    ):
        receipt.pop(field)
    receipt["search"]["basin_count"] = len(diagnostics.search.basins)
    for field in ("attempts", "basins", "budget_prefixes"):
        receipt["search"].pop(field)
    return receipt


def _evaluate_model(
    epcsaft: Any,
    equilibrium: Any,
    model_id: str,
    row: dict[str, str],
    contract: dict[str, Any],
    work: Path,
) -> dict[str, Any]:
    if str(COMPARISON_SCRIPTS) not in sys.path:
        sys.path.insert(0, str(COMPARISON_SCRIPTS))
    import run_polar_comparison as variants

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
        "failure_kind": "",
        "chemical_certification_level": "",
        "numerical_status": "",
        "physical_status": "",
        "local_minimum_status": "",
        "first_failed_numerical_criterion": "",
        "first_failed_physical_criterion": "",
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
        diagnostics = solved.diagnostics
        result.update(
            status="CERTIFIED_LOCAL_EQUILIBRIUM",
            equilibrium_predicted_pco2_pa=pressure,
            log10_equilibrium_over_legacy=math.log10(
                pressure / result["legacy_predicted_pco2_pa"]
            ),
            solver_status=diagnostics.solver_status,
            failure_kind=diagnostics.failure_kind,
            chemical_certification_level=diagnostics.chemical_certification_level,
            numerical_status=diagnostics.numerical_status,
            physical_status=diagnostics.physical_status,
            local_minimum_status=diagnostics.local_minimum_status,
            first_failed_numerical_criterion=(
                diagnostics.first_failed_numerical_criterion or ""
            ),
            first_failed_physical_criterion=(
                diagnostics.first_failed_physical_criterion or ""
            ),
            chemical_diagnostics=_diagnostic_receipt(diagnostics),
            equilibrium_mole_fractions=list(solved.mole_fractions),
            ln_k_provider_basis=list(solved.ln_k_provider_basis or ()),
            equilibrium_artifact_identity=asdict(solved.artifact_identity),
        )
    except equilibrium.ChemicalEquilibriumError as error:
        diagnostics = error.diagnostics
        result.update(
            status=_failure_class(diagnostics.failure_kind),
            reason=str(error),
            solver_status=diagnostics.solver_status,
            failure_kind=diagnostics.failure_kind,
            chemical_certification_level=diagnostics.chemical_certification_level,
            numerical_status=diagnostics.numerical_status,
            physical_status=diagnostics.physical_status,
            local_minimum_status=diagnostics.local_minimum_status,
            first_failed_numerical_criterion=(
                diagnostics.first_failed_numerical_criterion or ""
            ),
            first_failed_physical_criterion=(
                diagnostics.first_failed_physical_criterion or ""
            ),
            provider_domain_status=diagnostics.provider_domain_status,
            chemical_diagnostics=_diagnostic_receipt(diagnostics),
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
        "failure_kind",
        "chemical_certification_level",
        "numerical_status",
        "physical_status",
        "local_minimum_status",
        "first_failed_numerical_criterion",
        "first_failed_physical_criterion",
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
    parser.add_argument("--provider-wheel", type=Path, required=True)
    parser.add_argument("--equilibrium-wheel", type=Path, required=True)
    parser.add_argument("--provider-wheel-sha256", required=True)
    parser.add_argument("--equilibrium-wheel-sha256", required=True)
    args = parser.parse_args()

    import epcsaft
    import epcsaft_equilibrium

    provider_source = _git_identity(args.provider_source, args.provider_ref)
    equilibrium_source = _git_identity(args.equilibrium_source, args.equilibrium_ref)
    installed_provider = _installed_identity(
        "epcsaft", "epcsaft", "epcsaft/include/epcsaft/native_sdk_v1.h"
    )
    installed_equilibrium = _installed_identity(
        "epcsaft-equilibrium",
        "epcsaft_equilibrium",
        "epcsaft_equilibrium/include/epcsaft/regression/evaluator_v1.h",
    )
    provider_wheel = _wheel_identity(
        args.provider_wheel,
        "epcsaft",
        "epcsaft/include/epcsaft/native_sdk_v1.h",
    )
    equilibrium_wheel = _wheel_identity(
        args.equilibrium_wheel,
        "epcsaft_equilibrium",
        "epcsaft_equilibrium/include/epcsaft/regression/evaluator_v1.h",
    )
    _verify_installed_artifact(
        installed_provider, provider_wheel, args.provider_wheel_sha256
    )
    _verify_installed_artifact(
        installed_equilibrium, equilibrium_wheel, args.equilibrium_wheel_sha256
    )
    for installed, source in (
        (installed_provider, args.provider_source),
        (installed_equilibrium, args.equilibrium_source),
    ):
        try:
            Path(installed["installation_root"]).resolve().relative_to(source.resolve())
        except ValueError:
            continue
        raise ValueError("runtime imported a package from a source checkout")
    for installed in (installed_provider, installed_equilibrium):
        installed.pop("installation_root")
        installed["source_checkout_imported"] = False

    contract = json.loads(REACTIONS.read_text(encoding="utf-8"))
    transfer_contract = json.loads(SOURCE_REFERENCE.read_text(encoding="utf-8"))
    capability = _source_reference_capability(
        epcsaft, epcsaft_equilibrium, contract, transfer_contract
    )
    legacy = _legacy_rows()
    RESULTS.mkdir(parents=True, exist_ok=True)
    if capability["status"] == "SUPPORTED":
        raise RuntimeError(
            "a public source-reference transfer is now advertised; explicitly wire "
            "and review that API before enabling chemistry"
        )
    blocked_status = str(transfer_contract["failure_status"])
    rows = [
        _blocked_source_reference_row(
            model_id,
            legacy[model_id],
            capability,
            blocked_status,
        )
        for model_id in args.models
    ]
    for row in rows:
        print(
            f"{row['model_id']}: {row['status']} (solver not called)",
            flush=True,
        )

    statuses = {str(row["status"]) for row in rows}
    if statuses == {blocked_status}:
        overall = blocked_status
    elif "BLOCKED_UNCLASSIFIED" in statuses:
        overall = "BLOCKED_UNCLASSIFIED"
    elif statuses == {"CERTIFIED_LOCAL_EQUILIBRIUM"}:
        overall = "PARITY_EVALUATED"
    elif any(status.startswith("BLOCKED_") for status in statuses):
        overall = "BLOCKED_SHARED_PROVIDER_CAPABILITY"
    elif "CERTIFIED_LOCAL_EQUILIBRIUM" in statuses:
        overall = "PARTIAL_LOCAL_EQUILIBRIUM_CERTIFICATION"
    else:
        overall = "NOT_CERTIFIED"
    _write_csv(RESULTS / "equilibrium_replay_comparison.csv", rows)
    receipt = {
        "analysis": "public ePC-SAFT Equilibrium replay",
        "overall_status": overall,
        "target_observation_id": TARGET_OBSERVATION,
        "model_order": list(args.models),
        "provider_source": provider_source,
        "equilibrium_source": equilibrium_source,
        "installed_provider": installed_provider,
        "installed_equilibrium": installed_equilibrium,
        "wheel_identity": {
            "epcsaft": provider_wheel,
            "epcsaft-equilibrium": equilibrium_wheel,
        },
        "source_reference_transfer": capability,
        "source_hashes": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (REACTIONS, SOURCE_REFERENCE, SENTINEL, LEGACY)
        },
        "results": rows,
        "claim_boundary": (
            "No corrected M5 state is reported unless an installed public API "
            "performs the exact declared-source-to-Provider reference transfer. "
            "A future successful row would establish one certified local fixed-T,P "
            "homogeneous state only, not global equilibrium, parameter validity, "
            "or regression readiness."
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
