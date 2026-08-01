# M0--M4 reactive-model comparison

This standalone analysis compares five deliberately limited MEA--H2O--CO2
model variants at 313.15 K. It does not alter the manuscript or promote a
parameter set.

| Model | Definition |
|---|---|
| M0 | Current nine-species regression-input bundle. |
| M1 | M0 with MEAH+ segment diameter and dispersion energy copied from neutral MEA. |
| M2 | M0 with fixed CO2--water induced association following Pabsch et al. (2020). |
| M3 | M2 with MEAH+ and MEACOO- segment diameters fitted to the frozen two-row pressure/speciation tracer. |
| M4 | M2 with a sensitivity-selected subset of the preregistered ionic coordinates fitted to the 313.15 K Hilliard pressure campaign. |

M1 is a controlled neutral-parent transfer experiment, not a literature
reproduction. M2 retains the current CO2 and water parameters and applies the
Pabsch induced-association topology and combining rules without fitting an
association parameter. M3 is a two-observation diagnostic fit. It is neither
independent validation nor evidence that the ion parameters are identifiable
over the full experimental domain.

The executed M3 candidate is retained even when it fails a preregistered
acceptance check. In particular, a parameter at either bound makes M3 a
bound-limited diagnostic rather than an accepted fit.

## M4 pressure sensitivity experiment

The repository contains 44 pCO2 observations at 313.15 K with both eligible
metrology and a complete model-state pressure, more than at any other
temperature. The frozen grouped split assigns 24 Hilliard observations at
30 wt% MEA to training. It reserves 20 observations for checks that do not
enter parameter selection or fitting: six Hilliard rows at 17 wt% MEA, six
Hilliard rows at 40 wt% MEA, and eight Jou rows at 30 wt% MEA.

The local screen includes the three preregistered coordinates (MEAH+ segment
diameter, MEAH+ dispersion energy, and MEACOO- segment diameter). It also
reports, but does not fit, both Born diameters, the fixed MEACOO- dispersion
energy, and five selected binary interaction parameters. Active-coordinate
rank is determined from six loading-spanning training states. The retained
coordinates are fitted to all 24 training pressures with two bounded
Gauss--Newton steps, equal row weights, no regularization, and no use of the
reserved rows.

The screen retains MEAH+ segment diameter and MEAH+ dispersion energy. The
MEACOO- diameter direction is strongly confounded with MEAH+ dispersion
energy (correlation -0.90). The fitted values are 3.44617 angstrom and
144.590 K, respectively. Training log10-RMSE decreases from 0.2957 to 0.2771,
but reserved log10-RMSE increases from 0.4094 to 0.4528. The fit mainly lowers
the pressure curve and does not remove its loading-dependent residual shape.
These values therefore remain diagnostic and must not replace the regression
input bundle.

The model comparison uses the public Provider fixed-density EOS and an
MEA-owned reactive-speciation solver. Loading, reaction constants, balances,
activity convention, observation identities, and residual definitions remain
application-owned. The two M3 residuals are

\[
r_1=\log_{10}(p_{CO_2}^{model}/574\,\mathrm{Pa}),\qquad
r_2=\log_{10}(x_{MEACOO^-}^{model}/0.0502).
\]

Run with an environment containing the clean `epcsaft` 0.2 provider:

```bash
PYTHONPATH=src python analyses/phase3/m0_m3_model_comparison/scripts/run_comparison.py
PYTHONPATH=src python analyses/phase3/m0_m3_model_comparison/scripts/render_figures.py
PYTHONPATH=src:analyses/phase3/m0_m3_model_comparison/scripts \
  python analyses/phase3/m0_m3_model_comparison/scripts/run_pressure_sensitivity_fit.py
PYTHONPATH=src python \
  analyses/phase3/m0_m3_model_comparison/scripts/render_pressure_sensitivity_fit.py
```

The result CSV files are the exact plotted data. The Amundsen density series is
context only: its source table does not report pressure, so it is not used in
the M3 objective. The four-panel pressure comparison uses the 24 Hilliard rows
at 313.15 K and 30 wt% MEA for which both calibrated CO2 partial pressure and
row-reported total pressure are available. Each model is evaluated at each
row's reported total pressure; sources without a complete pressure state are
not mixed into this figure.

The M4 CSV tables retain the complete 313.15 K row inventory, the exact
screening perturbations and sensitivities, the fitted values, every plotted
prediction, and training/reserved metrics. The sensitivity figure plots the
finite log-pressure response to each declared perturbation; its horizontal
scale is therefore comparable across parameter families even though their
derivative coordinate scales differ. It is a response to those explicit step
choices, not a scale-free ranking of global parameter importance.
