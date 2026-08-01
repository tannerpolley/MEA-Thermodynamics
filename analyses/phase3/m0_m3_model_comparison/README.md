# M0--M3 reactive-model comparison

This standalone analysis compares four deliberately limited MEA--H2O--CO2
model variants at 313.15 K. It does not alter the manuscript or promote a
parameter set.

| Model | Definition |
|---|---|
| M0 | Current nine-species regression-input bundle. |
| M1 | M0 with MEAH+ segment diameter and dispersion energy copied from neutral MEA. |
| M2 | M0 with fixed CO2--water induced association following Pabsch et al. (2020). |
| M3 | M2 with MEAH+ and MEACOO- segment diameters fitted to the frozen two-row pressure/speciation tracer. |

M1 is a controlled neutral-parent transfer experiment, not a literature
reproduction. M2 retains the current CO2 and water parameters and applies the
Pabsch induced-association topology and combining rules without fitting an
association parameter. M3 is a two-observation diagnostic fit. It is neither
independent validation nor evidence that the ion parameters are identifiable
over the full experimental domain.

The executed M3 candidate is retained even when it fails a preregistered
acceptance check. In particular, a parameter at either bound makes M3 a
bound-limited diagnostic rather than an accepted fit.

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
```

The result CSV files are the exact plotted data. The Amundsen density series is
context only: its source table does not report pressure, so it is not used in
the M3 objective. The four-panel pressure comparison uses the 24 Hilliard rows
at 313.15 K and 30 wt% MEA for which both calibrated CO2 partial pressure and
row-reported total pressure are available. Each model is evaluated at each
row's reported total pressure; sources without a complete pressure state are
not mixed into this figure.
