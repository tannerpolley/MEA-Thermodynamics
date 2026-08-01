from __future__ import annotations

import hashlib
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import render_pressure_sensitivity_fit as rendering

RESULTS = rendering.RESULTS
FIGURES = rendering.FIGURES
MODEL_ORDER = ("M0", "M1", "M2", "M3", "M4A", "M4B", "M5Q", "M5")
MODEL_LABELS = {
    "M0": "base",
    "M1": "neutral-parent transfer",
    "M2": "induced association",
    "M3": "two-row ionic fit",
    "M4A": "Hilliard pressure fit",
    "M4B": "Hilliard + Jou fit",
    "M5Q": r"M2 + CO$_2$ QQ",
    "M5": "full DD + QQ + DQ",
}
COLORS = {
    "M0": "#6F6F6F",
    "M1": "#8C6BB1",
    "M2": "#2A6FBB",
    "M3": "#7A9E2A",
    "M4A": "#E07A2D",
    "M4B": "#D1495B",
    "M5Q": "#00A6A6",
    "M5": "#8B1E3F",
}
plt.rcParams["svg.hashsalt"] = "mea-m0-m5-polar-comparison"

_rows = rendering._rows
_save = rendering._save


def _metric_lookup() -> dict[tuple[str, str], float]:
    return {
        (row["model_id"], row["selection"]): float(row["log10_rmse"])
        for row in _rows(RESULTS / "m0_m5_pressure_metrics.csv")
    }


def _observations(axis: plt.Axes, rows: list[dict[str, str]]) -> None:
    hilliard = [row for row in rows if row["source_key"] == "Hilliard2008"]
    jou = [row for row in rows if row["source_key"] == "Jou1995"]
    axis.scatter(
        [float(row["loading_mol_co2_per_mol_mea"]) for row in hilliard],
        [float(row["observed_pco2_pa"]) / 1000.0 for row in hilliard],
        marker="o",
        facecolors="none",
        edgecolors="#222222",
        linewidths=0.9,
        s=24,
        zorder=3,
    )
    axis.scatter(
        [float(row["loading_mol_co2_per_mol_mea"]) for row in jou],
        [float(row["observed_pco2_pa"]) / 1000.0 for row in jou],
        marker="x",
        color="#222222",
        linewidths=0.9,
        s=24,
        zorder=3,
    )


def _pressure_panels() -> None:
    rows = _rows(RESULTS / "m0_m5_pressure_predictions.csv")
    rows = [row for row in rows if np.isclose(float(row["mea_mass_fraction"]), 0.30)]
    metrics = _metric_lookup()
    fig, axes = plt.subplots(2, 4, figsize=(12.2, 6.6), sharex=True, sharey=True)
    for axis, model_id in zip(axes.flat, MODEL_ORDER, strict=True):
        selected = [row for row in rows if row["model_id"] == model_id]
        _observations(axis, selected)
        for source in ("Hilliard2008", "Jou1995"):
            source_rows = sorted(
                (row for row in selected if row["source_key"] == source),
                key=lambda row: float(row["loading_mol_co2_per_mol_mea"]),
            )
            axis.plot(
                [float(row["loading_mol_co2_per_mol_mea"]) for row in source_rows],
                [float(row["predicted_pco2_pa"]) / 1000.0 for row in source_rows],
                color=COLORS[model_id],
                linestyle="--",
                linewidth=1.45,
                alpha=1.0 if source == "Hilliard2008" else 0.65,
            )
        axis.set_title(
            f"{model_id}: {MODEL_LABELS[model_id]}\nRMSE = {metrics[(model_id, 'all-30wt')]:.3f}",
            fontsize=9.5,
        )
        axis.set_yscale("log")
        axis.grid(alpha=0.18, which="both")
        axis.spines[["top", "right"]].set_visible(False)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$p_{CO_2}$ (kPa)")
    fig.legend(
        handles=(
            Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor="#222222", linestyle="", label="Hilliard observed"),
            Line2D([], [], marker="x", color="#222222", linestyle="", label="Jou observed"),
            Line2D([], [], color="#555555", linestyle="--", label="model"),
        ),
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=3,
        frameon=False,
    )
    fig.suptitle("313.15 K, 30 wt% MEA: all reactive-model variants", y=1.08)
    fig.subplots_adjust(hspace=0.35, wspace=0.14)
    _save(fig, "m0_m5_pco2_30wt")


def _residual_panels() -> None:
    rows = _rows(RESULTS / "m0_m5_pressure_predictions.csv")
    rows = [row for row in rows if np.isclose(float(row["mea_mass_fraction"]), 0.30)]
    fig, axes = plt.subplots(2, 4, figsize=(12.2, 6.2), sharex=True, sharey=True)
    for axis, model_id in zip(axes.flat, MODEL_ORDER, strict=True):
        selected = [row for row in rows if row["model_id"] == model_id]
        for source, marker in (("Hilliard2008", "o"), ("Jou1995", "x")):
            source_rows = [row for row in selected if row["source_key"] == source]
            axis.scatter(
                [float(row["loading_mol_co2_per_mol_mea"]) for row in source_rows],
                [float(row["log10_residual"]) for row in source_rows],
                marker=marker,
                facecolors="none" if marker == "o" else COLORS[model_id],
                edgecolors=COLORS[model_id] if marker == "o" else None,
                color=COLORS[model_id] if marker == "x" else None,
                s=26,
                linewidths=0.9,
            )
        ordered = sorted(selected, key=lambda row: float(row["loading_mol_co2_per_mol_mea"]))
        axis.plot(
            [float(row["loading_mol_co2_per_mol_mea"]) for row in ordered],
            [float(row["log10_residual"]) for row in ordered],
            color=COLORS[model_id],
            linestyle="--",
            linewidth=0.9,
            alpha=0.45,
        )
        axis.axhline(0.0, color="#333333", linewidth=0.75)
        axis.set_title(f"{model_id}: {MODEL_LABELS[model_id]}", fontsize=9.5)
        axis.grid(alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"loading (mol $CO_2$/mol MEA)")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\log_{10}(p_{model}/p_{obs})$")
    fig.suptitle("Residual curvature across all 30 wt% model variants", y=1.02)
    fig.subplots_adjust(hspace=0.30, wspace=0.14)
    _save(fig, "m0_m5_residuals_30wt")


def _metric_scatter() -> None:
    metrics = _metric_lookup()
    x = np.arange(len(MODEL_ORDER))
    fig, axis = plt.subplots(figsize=(9.4, 4.4))
    for selection, marker, label, color in (
        ("all-30wt", "o", "all 30 wt%", "#222222"),
        ("Hilliard-30wt", "s", "Hilliard", "#2A6FBB"),
        ("Jou-30wt", "^", "Jou", "#D1495B"),
    ):
        axis.scatter(
            x,
            [metrics[(model_id, selection)] for model_id in MODEL_ORDER],
            marker=marker,
            color=color,
            s=45,
            label=label,
            zorder=3,
        )
    axis.set_yscale("log")
    axis.set_xticks(x, MODEL_ORDER)
    axis.set_ylabel(r"RMSE of $\log_{10}(p_{model}/p_{obs})$")
    axis.set_xlabel("model variant")
    axis.set_title("Pressure error by source at 313.15 K and 30 wt% MEA")
    axis.grid(axis="y", alpha=0.2, which="both")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, ncol=3)
    _save(fig, "m0_m5_metric_scatter")


def main() -> None:
    _pressure_panels()
    _residual_panels()
    _metric_scatter()
    audit_path = RESULTS / "m5_polar_capability_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["figure_outputs"] = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(FIGURES.glob("m0_m5_*.*"))
    }
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
