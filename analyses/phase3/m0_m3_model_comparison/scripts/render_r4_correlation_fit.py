from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from MEA.common.config import REPO_ROOT


ANALYSIS = REPO_ROOT / "analyses/phase3/m0_m3_model_comparison"
RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"
MODEL_STYLES = {
    "M5_literature_R4": {"label": "Literature R4", "color": "#9C3D10", "marker": "s"},
    "M5_fitted_R4": {"label": "Fitted R4", "color": "#006D8F", "marker": "D"},
}
SOURCE_STYLES = {
    "Hilliard2008": {"label": "Hilliard", "color": "#1B6CA8", "marker": "o"},
    "Jou1995": {"label": "Jou", "color": "#D97706", "marker": "s"},
    "Xu2011": {"label": "Xu", "color": "#4C956C", "marker": "^"},
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _finish(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "png", "pdf"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _style(axis: plt.Axes) -> None:
    axis.grid(color="#D8D8D8", linewidth=0.6, alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _parity(rows: list[dict[str, str]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), sharex=True, sharey=True)
    all_values = []
    for axis, (model_id, model_style) in zip(axes, MODEL_STYLES.items(), strict=True):
        selected = [row for row in rows if row["model_id"] == model_id]
        for source, source_style in SOURCE_STYLES.items():
            source_rows = [row for row in selected if row["source_key"] == source]
            if not source_rows:
                continue
            observed = np.asarray(
                [_float(row, "observed_pco2_pa") for row in source_rows]
            )
            predicted = np.asarray(
                [_float(row, "predicted_pco2_pa") for row in source_rows]
            )
            all_values.extend(observed)
            all_values.extend(predicted)
            axis.scatter(
                observed,
                predicted,
                s=28,
                marker=source_style["marker"],
                facecolors="none",
                edgecolors=source_style["color"],
                linewidths=1.0,
                label=source_style["label"],
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(model_style["label"])
        axis.set_xlabel(r"Measured $p_{CO_2}$ (Pa)")
        _style(axis)
    lower = 10.0 ** np.floor(np.log10(min(all_values)))
    upper = 10.0 ** np.ceil(np.log10(max(all_values)))
    for axis in axes:
        axis.plot(
            [lower, upper], [lower, upper], color="black", linewidth=0.9, linestyle="--"
        )
        axis.set_xlim(lower, upper)
        axis.set_ylim(lower, upper)
    axes[0].set_ylabel(r"Modeled liquid $CO_2$ fugacity (Pa)")
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle(
        "Reactive M5 pressure comparison: 116 evaluated of 121 admitted states",
        y=1.01,
    )
    fig.tight_layout()
    _finish(fig, "r4_correlation_fit_pco2")


def _residual_structure(rows: list[dict[str, str]]) -> None:
    rows = [row for row in rows if row["model_id"] == "M5_fitted_R4"]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.0))
    x_fields = (
        ("temperature_c", r"Temperature ($^{\circ}$C)"),
        ("loading_mol_co2_per_mol_mea", r"$CO_2$ loading (mol/mol MEA)"),
        ("observed_pco2_pa", r"Measured $p_{CO_2}$ (Pa)"),
        ("mea_mass_fraction", "Unloaded MEA mass fraction"),
    )
    for axis, (field, label) in zip(axes.flat, x_fields, strict=True):
        for source, source_style in SOURCE_STYLES.items():
            selected = [row for row in rows if row["source_key"] == source]
            if not selected:
                continue
            axis.scatter(
                [_float(row, field) for row in selected],
                [_float(row, "log10_pressure_residual") for row in selected],
                s=28,
                marker=source_style["marker"],
                facecolors="none",
                edgecolors=source_style["color"],
                linewidths=1.0,
                label=source_style["label"],
            )
        axis.axhline(0.0, color="black", linewidth=0.9)
        axis.set_xlabel(label)
        axis.set_ylabel(r"$\log_{10}(p_{model}/p_{data})$")
        if field == "observed_pco2_pa":
            axis.set_xscale("log")
        _style(axis)
    axes[0, 1].legend(frameon=False, fontsize=9)
    fig.suptitle("Residual structure after the diagnostic R4 fit")
    fig.tight_layout()
    _finish(fig, "r4_correlation_fit_residual_structure")


def _carbamate(rows: list[dict[str, str]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.1), sharex=True, sharey=True)
    maximum = 0.0
    for axis, (model_id, model_style) in zip(axes, MODEL_STYLES.items(), strict=True):
        selected = [
            row
            for row in rows
            if row["model_id"] == model_id and row["carbamate_target_mole_fraction"]
        ]
        for source, source_style in SOURCE_STYLES.items():
            source_rows = [row for row in selected if row["source_key"] == source]
            if not source_rows:
                continue
            observed = [
                _float(row, "carbamate_target_mole_fraction") for row in source_rows
            ]
            predicted = [
                _float(row, "predicted_carbamate_mole_fraction") for row in source_rows
            ]
            maximum = max(maximum, *observed, *predicted)
            axis.scatter(
                observed,
                predicted,
                s=30,
                marker=source_style["marker"],
                facecolors="none",
                edgecolors=source_style["color"],
                linewidths=1.0,
                label=source_style["label"],
            )
        axis.set_title(model_style["label"])
        axis.set_xlabel(r"Interpolated Böttinger $x_{MEACOO^-}$")
        _style(axis)
    for axis in axes:
        axis.plot(
            [0.0, maximum], [0.0, maximum], color="black", linewidth=0.9, linestyle="--"
        )
        axis.set_xlim(0.0, maximum * 1.04)
        axis.set_ylim(0.0, maximum * 1.04)
    axes[0].set_ylabel(r"Modeled $x_{MEACOO^-}$")
    axes[1].legend(frameon=False, fontsize=9)
    fig.suptitle("Carbamate comparison excluded from the pressure objective", y=1.01)
    fig.tight_layout()
    _finish(fig, "r4_correlation_fit_carbamate")


def _sensitivity() -> None:
    reaction = _read(RESULTS / "r4_reaction_sensitivity_rows.csv")
    eos = _read(RESULTS / "r4_eos_sensitivity_rows.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
    reactions = sorted({row["reaction_id"] for row in reaction})
    reaction_p = [
        np.median(
            [
                abs(_float(row, "dlog10_pco2_dlnk"))
                for row in reaction
                if row["reaction_id"] == key
            ]
        )
        for key in reactions
    ]
    reaction_x = [
        np.median(
            [
                abs(_float(row, "dlog10_meacoo_dlnk"))
                for row in reaction
                if row["reaction_id"] == key
            ]
        )
        for key in reactions
    ]
    locations = np.arange(len(reactions))
    width = 0.36
    axes[0].bar(
        locations - width / 2, reaction_p, width, color="#006D8F", label=r"$p_{CO_2}$"
    )
    axes[0].bar(
        locations + width / 2,
        reaction_x,
        width,
        color="#D97706",
        label=r"$x_{MEACOO^-}$",
    )
    axes[0].set_xticks(locations, reactions)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Median absolute log-response")
    axes[0].set_title("Reaction-constant sensitivities")
    axes[0].legend(frameon=False, fontsize=9)
    _style(axes[0])

    parameters = sorted({row["parameter"] for row in eos})
    pressure = [
        np.median(
            [
                abs(_float(row, "dlog10_observable_daffine_coordinate"))
                for row in eos
                if row["parameter"] == key and row["observable"] == "pco2"
            ]
        )
        for key in parameters
    ]
    carbamate = [
        np.median(
            [
                abs(_float(row, "dlog10_observable_daffine_coordinate"))
                for row in eos
                if row["parameter"] == key
                and row["observable"] == "carbamate_mole_fraction"
            ]
        )
        for key in parameters
    ]
    locations = np.arange(len(parameters))
    axes[1].bar(
        locations - width / 2, pressure, width, color="#006D8F", label=r"$p_{CO_2}$"
    )
    axes[1].bar(
        locations + width / 2,
        carbamate,
        width,
        color="#D97706",
        label=r"$x_{MEACOO^-}$",
    )
    axes[1].set_xticks(locations, [value.replace("::", "\n") for value in parameters])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Median absolute scaled response")
    axes[1].set_title("Selected EOS-parameter sensitivities")
    axes[1].legend(frameon=False, fontsize=9)
    _style(axes[1])
    fig.tight_layout()
    _finish(fig, "r4_correlation_fit_sensitivity")


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
        }
    )
    rows = _read(RESULTS / "r4_correlation_fit_rows.csv")
    _parity(rows)
    _residual_structure(rows)
    _carbamate(rows)
    _sensitivity()


if __name__ == "__main__":
    main()
