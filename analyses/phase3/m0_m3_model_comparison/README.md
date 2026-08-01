# M0--M4B reactive-model comparison

This standalone analysis compares six deliberately limited MEA--H2O--CO2
model variants at 313.15 K. It does not alter the manuscript or promote a
parameter set.

| Model | Definition |
|---|---|
| M0 | Current nine-species regression-input bundle. |
| M1 | M0 with MEAH+ segment diameter and dispersion energy copied from neutral MEA. |
| M2 | M0 with fixed CO2--water induced association following Pabsch et al. (2020). |
| M3 | M2 with MEAH+ and MEACOO- segment diameters fitted to the frozen two-row pressure/speciation tracer. |
| M4A | M2 with a sensitivity-selected subset of the preregistered ionic coordinates fitted to the 313.15 K Hilliard pressure campaign. |
| M4B | The M4A coordinates fitted jointly to the 313.15 K Hilliard and Jou 30 wt% pressure campaigns with equal total weight per source. |

M1 is a controlled neutral-parent transfer experiment, not a literature
reproduction. M2 retains the current CO2 and water parameters and applies the
Pabsch induced-association topology and combining rules without fitting an
association parameter. M3 is a two-observation diagnostic fit. It is neither
independent validation nor evidence that the ion parameters are identifiable
over the full experimental domain.

The executed M3 candidate is retained even when it fails a preregistered
acceptance check. In particular, a parameter at either bound makes M3 a
bound-limited diagnostic rather than an accepted fit.

## M4 pressure sensitivity experiments

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
coordinates are first fitted to all 24 Hilliard training pressures as M4A,
using two bounded Gauss--Newton steps, equal row weights, and no
regularization. M4B then fits the same coordinates to the 24 Hilliard and eight
Jou 30 wt% rows. Its per-row residual multipliers are
\(\sqrt{2/3}\) for Hilliard and \(\sqrt{2}\) for Jou, giving the two sources
equal total squared weight despite their different row counts.

The screen retains MEAH+ segment diameter and MEAH+ dispersion energy. The
MEACOO- diameter direction is strongly confounded with MEAH+ dispersion
energy (correlation -0.90). The fitted values are 3.44617 angstrom and
144.590 K, respectively, for M4A. Training log10-RMSE decreases from 0.2957 to 0.2771,
but reserved log10-RMSE increases from 0.4094 to 0.4528. The fit mainly lowers
the pressure curve and does not remove its loading-dependent residual shape.
M4B intentionally consumes the Jou source holdout and therefore cannot support
the source-transfer check used for M4A. The 17 and 40 wt% Hilliard campaigns
remain untouched concentration validation for both fits. All fitted values
remain diagnostic and must not replace the regression input bundle.

M4B moves the MEAH+ diameter to 3.14061 angstrom and its dispersion energy to
336.295 K. The equal-source 30 wt% objective log10-RMSE is 0.3183, compared
with 0.3342 for M2 and 0.3715 for M4A. The gain comes mainly from recovering
the Jou campaign: its unweighted log10-RMSE decreases from 0.3687 for M2 and
0.4464 for M4A to 0.3408. The Hilliard 30 wt% error returns from M4A's 0.2771
to 0.2941, nearly the M2 value of 0.2957. On the untouched 17 and 40 wt%
campaigns, M4B has log10-RMSE 0.4463 versus 0.4344 for M2. Joint-source fitting
therefore chooses a different ionic-parameter compromise but still does not
remove the loading-dependent curvature or improve concentration transfer.

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

The M4A/M4B CSV tables retain the complete 313.15 K row inventory, the exact
screening perturbations and sensitivities, the fitted values, every plotted
prediction, and group-specific metrics. The sensitivity figure plots the
finite log-pressure response to each declared perturbation; its horizontal
scale is therefore comparable across parameter families even though their
derivative coordinate scales differ. It is a response to those explicit step
choices, not a scale-free ranking of global parameter importance.

## M5 polar-capability audit

M5 would add the Gross--Vrabec dipole--dipole, quadrupole--quadrupole, and
dipole--quadrupole residual Helmholtz contributions to the neutral H2O, CO2,
and MEA model. It is not executable in the installed EOS used here. The current
dipole and polarizability inputs feed the mixture-permittivity and Born path;
they do not activate a residual polar attraction. EOS issues 51--54 define the
pending aggregate `ares_polar` implementation.

Gross (2005) supplies a directly usable PCP-SAFT CO2 candidate fitted with the
quadrupole contribution active: \(m=1.5131\), \(\sigma=3.1869\) angstrom,
\(\epsilon/k=163.33\) K, and \(|Q|=4.4\) D angstrom. The audited sources do not
supply an equivalently coherent dipolar-associating refit for both water and
MEA. Schick et al. (2023), for example, retains ordinary associating water,
uses induced association for CO2, and activates a dipolar term for NMP rather
than MEA. A future M5 comparison must therefore refit or independently qualify
the water and MEA neutral parameters with the polar equations active; combining
the CO2 set with the current nonpolar water and MEA sets would not be a
controlled model comparison.

The resulting staged route is M5a (CO2 QQ diagnostic after EOS issue 52), M5b
(neutral water and MEA refit with DD after issue 53), and M5c (DQ cross terms
after issue 54). The first stage can isolate whether explicit CO2 quadrupolar
physics changes the pressure-curve shape. It cannot, by itself, qualify an
all-polar MEA--H2O--CO2 parameterization.
