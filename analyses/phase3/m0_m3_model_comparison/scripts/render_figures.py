from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]
ANALYSIS = ROOT / "analyses/phase3/m0_m3_model_comparison"
RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"
AMUNDSEN = (
    ROOT
    / "data/reference/MEA/observations/density_viscosity"
    / "Amundsen_2009_density_viscosity.csv"
)
COLORS = {"M0": "#4C78A8", "M1": "#F58518", "M2": "#54A24B", "M3": "#B279A2"}
DISPLAY = {"M0": "M0", "M1": "M1", "M2": "M2", "M3": "M3 (bound-limited)"}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "png", "pdf"):
        path = FIGURES / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        if suffix == "svg":
            path.write_text(
                "\n".join(
                    line.rstrip()
                    for line in path.read_text(encoding="utf-8").splitlines()
                )
                + "\n",
                encoding="utf-8",
            )
    plt.close(fig)


def _pco2_loading() -> None:
    rows = _rows(RESULTS / "pco2_loading_comparison.csv")
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 6.2), sharex=True, sharey=True)
    for axis, model_id in zip(axes.flat, COLORS, strict=True):
        selected = [row for row in rows if row["model_id"] == model_id]
        loading = [float(row["loading"]) for row in selected]
        observed = [float(row["observed_pco2_pa"]) / 1000.0 for row in selected]
        predicted = [float(row["predicted_pco2_pa"]) / 1000.0 for row in selected]
        axis.scatter(
            loading,
            observed,
            facecolors="none",
            edgecolors="black",
            linewidths=1.0,
            s=30,
            label="Hilliard (2008)",
        )
        axis.scatter(
            loading,
            predicted,
            color=COLORS[model_id],
            marker="x",
            s=30,
            label="model",
        )
        axis.set_title(DISPLAY[model_id])
        axis.set_yscale("log")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.2, which="both")
    for axis in axes[-1, :]:
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$p_{CO_2}$ (kPa)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.95),
        ncol=2,
        frameon=False,
    )
    fig.suptitle(
        "313.15 K, 30 wt% MEA: observed and predicted $p_{CO_2}$",
        y=0.995,
    )
    fig.subplots_adjust(top=0.85, hspace=0.28, wspace=0.16)
    _save(fig, "m0_m3_pco2_loading_comparison")


def _species() -> None:
    rows = _rows(RESULTS / "species_predictions.csv")
    species = ("MEACOO-", "HCO3-", "CO2")
    titles = (r"$MEACOO^-$", r"$HCO_3^-$", r"molecular $CO_2$")
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.3), sharex=True)
    for axis, target, title in zip(axes, species, titles, strict=True):
        for model_id in COLORS:
            selected = [
                row
                for row in rows
                if row["model_id"] == model_id and row["species"] == target
            ]
            axis.plot(
                [float(row["loading_mol_co2_per_mol_mea"]) for row in selected],
                [float(row["mole_fraction"]) for row in selected],
                marker="o",
                color=COLORS[model_id],
                label=DISPLAY[model_id],
            )
        if target == "MEACOO-":
            axis.scatter([0.466], [0.0502], marker="x", s=70, color="black", label="observed")
        axis.set_title(title)
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("liquid mole fraction")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False)
    fig.subplots_adjust(top=0.78, wspace=0.28)
    _save(fig, "m0_m3_species_by_loading")


def _density() -> None:
    model_rows = _rows(RESULTS / "loading_predictions.csv")
    source_rows = [
        row
        for row in _rows(AMUNDSEN)
        if row["property"] == "density"
        and row["temperature_C"] == "40"
        and row["mea_mass_fraction"] == "0.30"
        and row["co2_loading_mol_per_mol_mea"]
    ]
    fig, axis = plt.subplots(figsize=(6.2, 3.8))
    for model_id in COLORS:
        selected = [row for row in model_rows if row["model_id"] == model_id]
        axis.plot(
            [float(row["loading_mol_co2_per_mol_mea"]) for row in selected],
            [float(row["density_kg_m3"]) for row in selected],
            marker="o",
            color=COLORS[model_id],
            label=DISPLAY[model_id],
        )
    axis.errorbar(
        [float(row["co2_loading_mol_per_mol_mea"]) for row in source_rows],
        [1000.0 * float(row["value"]) for row in source_rows],
        yerr=[1000.0 * float(row["uncertainty_value"]) for row in source_rows],
        fmt="kx",
        capsize=3,
        label="Amundsen (pressure not reported)",
    )
    axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
    axis.set_ylabel(r"liquid density (kg m$^{-3}$)")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, ncol=2)
    axis.set_title("Density context at 313.15 K and 30 wt% MEA; not fitted")
    _save(fig, "m0_m3_density_context")


def main() -> None:
    _pco2_loading()
    _species()
    _density()
    receipt_path = RESULTS / "comparison_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["outputs"] = {
        name: _sha256(RESULTS / name)
        for name in (
            "model_summary.csv",
            "species_predictions.csv",
            "loading_predictions.csv",
            "pco2_loading_comparison.csv",
        )
    }
    receipt["figure_outputs"] = {
        path.name: _sha256(path)
        for path in sorted(FIGURES.iterdir())
        if path.suffix in {".svg", ".png", ".pdf"}
    }
    receipt["analysis_file_hashes"] = {
        str(path.relative_to(ROOT)): _sha256(path)
        for path in (
            ANALYSIS / "README.md",
            ANALYSIS / "analysis.yaml",
            Path(__file__).with_name("run_comparison.py"),
            Path(__file__),
        )
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
