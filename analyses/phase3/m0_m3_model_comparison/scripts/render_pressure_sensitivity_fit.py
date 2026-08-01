from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[4]
ANALYSIS = ROOT / "analyses/phase3/m0_m3_model_comparison"
RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"
MODELS = {"M2": ("#777777", "M2"), "M4": ("#2A6FBB", "M4 pressure fit")}
FRACTIONS = (0.17, 0.30, 0.40)
plt.rcParams["svg.hashsalt"] = "mea-m4-pressure-sensitivity"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _save(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "png", "pdf"):
        metadata = (
            {"Date": None}
            if suffix == "svg"
            else {"CreationDate": None, "ModDate": None}
            if suffix == "pdf"
            else {"Software": "MEA-Thermodynamics"}
        )
        path = FIGURES / f"{stem}.{suffix}"
        fig.savefig(
            path,
            dpi=300,
            bbox_inches="tight",
            metadata=metadata,
        )
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


def _fit_by_composition() -> None:
    rows = _rows(RESULTS / "pressure_fit_predictions.csv")
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.1), sharey=True)
    for axis, fraction in zip(axes, FRACTIONS, strict=True):
        observed = [
            row
            for row in rows
            if row["model_id"] == "M2"
            and np.isclose(float(row["mea_mass_fraction"]), fraction)
        ]
        for source in sorted({row["source_key"] for row in observed}):
            source_rows = [row for row in observed if row["source_key"] == source]
            validation = source_rows[0]["split"] == "validation"
            axis.scatter(
                [float(row["loading_mol_co2_per_mol_mea"]) for row in source_rows],
                [float(row["observed_pco2_pa"]) / 1000.0 for row in source_rows],
                marker="x" if source == "Jou1995" else "o",
                facecolors="none" if source != "Jou1995" else None,
                edgecolors="#222222" if source != "Jou1995" else None,
                color="#222222" if source == "Jou1995" else None,
                linewidths=1.0,
                s=28,
                label=f"{source.replace('2008', '').replace('1995', '')} ({'validation' if validation else 'training'})",
                zorder=3,
            )
        for model_id, (color, label) in MODELS.items():
            model_rows = [
                row
                for row in rows
                if row["model_id"] == model_id
                and np.isclose(float(row["mea_mass_fraction"]), fraction)
            ]
            for source in sorted({row["source_key"] for row in model_rows}):
                selected = sorted(
                    (row for row in model_rows if row["source_key"] == source),
                    key=lambda row: float(row["loading_mol_co2_per_mol_mea"]),
                )
                axis.plot(
                    [float(row["loading_mol_co2_per_mol_mea"]) for row in selected],
                    [float(row["predicted_pco2_pa"]) / 1000.0 for row in selected],
                    color=color,
                    linestyle="--",
                    linewidth=1.6,
                    alpha=0.7 if source == "Jou1995" else 1.0,
                    label=label if source == sorted({row["source_key"] for row in model_rows})[0] else None,
                )
        role = (
            "training + Jou validation"
            if np.isclose(fraction, 0.30)
            else "validation"
        )
        axis.set_title(f"{100 * fraction:.0f} wt% MEA\n{role}")
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
        axis.set_yscale("log")
        axis.grid(alpha=0.2, which="both")
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$p_{CO_2}$ (kPa)")
    handles = []
    labels = []
    for axis in axes:
        axis_handles, axis_labels = axis.get_legend_handles_labels()
        handles.extend(axis_handles)
        labels.extend(axis_labels)
    unique = dict(zip(labels, handles, strict=True))
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.89),
        ncol=5,
        frameon=False,
        fontsize=8.5,
    )
    fig.suptitle("313.15 K pressure fit and untouched concentration/source checks", y=0.98)
    fig.subplots_adjust(top=0.68, wspace=0.12)
    _save(fig, "m4_pressure_fit_by_composition")


def _sensitivity() -> None:
    summary = _rows(RESULTS / "pressure_parameter_sensitivity.csv")
    fit = {row["parameter_identity"]: row for row in _rows(RESULTS / "pressure_fit_parameters.csv")}
    ordered = sorted(
        summary,
        key=lambda row: float(row["rms_log10_response_for_screen_step"]),
    )
    active_ids = [
        row["parameter_identity"]
        for row in summary
        if row["fit_eligible"] == "true"
    ]
    correlations = _rows(RESULTS / "pressure_active_correlation.csv")
    matrix = np.asarray(
        [
            [
                float(
                    next(
                        row["correlation"]
                        for row in correlations
                        if row["parameter_a"] == left and row["parameter_b"] == right
                    )
                )
                for right in active_ids
            ]
            for left in active_ids
        ]
    )
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.0, 4.5), gridspec_kw={"width_ratios": [1.55, 1]})
    y = np.arange(len(ordered))
    colors = [
        "#2A6FBB"
        if fit[row["parameter_identity"]]["selected_for_fit"] == "true"
        else "#E07A2D"
        if row["fit_eligible"] == "true"
        else "#888888"
        for row in ordered
    ]
    left.hlines(
        y,
        0.0,
        [float(row["rms_log10_response_for_screen_step"]) for row in ordered],
        color=colors,
        linewidth=1.2,
    )
    left.scatter(
        [float(row["rms_log10_response_for_screen_step"]) for row in ordered],
        y,
        color=colors,
        s=35,
        zorder=3,
    )
    left.set_yticks(y, [row["label"] for row in ordered])
    left.set_xlabel("RMS log-pressure response to declared perturbation")
    left.set_title("Six-state screening")
    left.grid(axis="x", alpha=0.2)
    left.spines[["top", "right"]].set_visible(False)
    left.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="", color="#2A6FBB", label="selected"),
            Line2D([], [], marker="o", linestyle="", color="#E07A2D", label="active, not selected"),
            Line2D([], [], marker="o", linestyle="", color="#888888", label="diagnostic only"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=8,
    )

    image = right.imshow(matrix, vmin=-1.0, vmax=1.0, cmap="coolwarm")
    short = ["MEAH+ sigma", "MEAH+ epsilon/k", "MEACOO- sigma"]
    right.set_xticks(range(len(short)), short, rotation=35, ha="right")
    right.set_yticks(range(len(short)), short)
    for row_index in range(len(short)):
        for column_index in range(len(short)):
            right.text(
                column_index,
                row_index,
                f"{matrix[row_index, column_index]:.2f}",
                ha="center",
                va="center",
                color="white" if abs(matrix[row_index, column_index]) > 0.55 else "black",
                fontsize=8,
            )
    right.set_title("Active-coordinate correlation")
    fig.colorbar(image, ax=right, fraction=0.046, pad=0.04)
    fig.suptitle("Which parameters can 313.15 K pressure data distinguish?", y=1.02)
    fig.subplots_adjust(wspace=0.55)
    _save(fig, "m4_pressure_parameter_sensitivity")


def _residuals() -> None:
    rows = _rows(RESULTS / "pressure_fit_predictions.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.55), sharex=True, sharey=True)
    markers = {0.17: "s", 0.30: "o", 0.40: "^"}
    for axis, model_id in zip(axes, MODELS, strict=True):
        selected = [row for row in rows if row["model_id"] == model_id]
        for fraction in FRACTIONS:
            group = [
                row
                for row in selected
                if np.isclose(float(row["mea_mass_fraction"]), fraction)
            ]
            for split, color in (("training", "#2A6FBB"), ("validation", "#D1495B")):
                subset = [row for row in group if row["split"] == split]
                if not subset:
                    continue
                axis.scatter(
                    [float(row["loading_mol_co2_per_mol_mea"]) for row in subset],
                    [float(row["log10_residual"]) for row in subset],
                    marker=markers[fraction],
                    facecolors="none",
                    edgecolors=color,
                    s=30,
                    label=f"{100 * fraction:.0f} wt%, {split}",
                )
        axis.axhline(0.0, color="#333333", linewidth=0.8)
        axis.set_title(MODELS[model_id][1])
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
        axis.grid(alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$\log_{10}(p_{model}/p_{obs})$")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=3, frameon=False)
    fig.suptitle("Residual shape before and after the pressure-only fit", y=1.18)
    fig.subplots_adjust(wspace=0.12)
    _save(fig, "m4_pressure_fit_residuals")


def main() -> None:
    _fit_by_composition()
    _sensitivity()
    _residuals()
    receipt_path = RESULTS / "pressure_fit_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    receipt["figure_outputs"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(FIGURES.glob("m4_pressure_*"))
        if path.suffix in {".svg", ".png", ".pdf"}
    }
    receipt["outputs"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(RESULTS.glob("pressure_*.csv"))
    }
    receipt["analysis_file_hashes"] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (
            ANALYSIS / "README.md",
            ANALYSIS / "analysis.yaml",
            Path(__file__).with_name("run_comparison.py"),
            Path(__file__).with_name("run_pressure_sensitivity_fit.py"),
            Path(__file__),
        )
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
