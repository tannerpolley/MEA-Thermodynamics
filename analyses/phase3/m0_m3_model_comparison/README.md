# M0--M5 reactive-model comparison

This standalone analysis compares eight deliberately limited MEA--H2O--CO2
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
| M5Q | M2 with the Gross (2005) QQ-consistent CO2 pure-component set and the measured 4.4 D angstrom quadrupole. |
| M5 | M5Q with fixed physical H2O and MEA dipoles, activating the complete neutral DD, QQ, and DQ contribution. |

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
PYTHONPATH=src:analyses/phase3/m0_m3_model_comparison/scripts \
  python analyses/phase3/m0_m3_model_comparison/scripts/run_polar_comparison.py
PYTHONPATH=src python \
  analyses/phase3/m0_m3_model_comparison/scripts/render_polar_comparison.py
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

## M5 polar experiments

The clean EOS at commit `b2638deb64772f2353f0382d0dc5a3210889a827`
implements the Gross--Vrabec DD, QQ, and DQ35 contributions and exposes their
sum as `State.ares_polar`. M5Q uses the Gross (2005) CO2 set fitted with QQ
active: \(m=1.5131\), \(\sigma=3.1869\) angstrom,
\(\epsilon/k=163.33\) K, and \(|Q|=4.4\) D angstrom. M5 adds fixed dipoles of
1.8546 D for water and 2.27 D for the stable gas-phase MEA conformer reported
by Tripathi (2016, doi:10.5821/dissertation-2117-106297, chapter 7). The MEA
value is a conformer-specific diagnostic, not a fitted liquid-phase effective
moment.

At 30 wt% MEA, M5Q has a log10-pressure RMSE of 0.7029, compared with 0.3156
for M2. It is close at the lowest Hilliard state (7.38 Pa predicted versus
5.7 Pa observed) but increasingly overpredicts pressure with loading, reaching
74.28 kPa versus 28.3 kPa at loading 0.591. The 17 and 40 wt% transfer RMSE is
0.5862, also worse than M2's 0.4344. The QQ-consistent CO2 set is therefore not
transferable into the existing reactive mixture model without revisiting its
unlike interactions and induced-association treatment.

Full M5 is intentionally more stringent: it activates DD and DQ while retaining
the existing water and MEA PC-SAFT and association parameters. It produces a
30 wt% RMSE of 2.7456 and overpredicts pressure by roughly two to three orders
of magnitude across the loading span. Its `ares_polar` ranges from -2.55 to
-2.02 over the executed states, whereas M5Q remains between approximately
-2.6e-7 and zero. The result is direct evidence that the old effective
dispersion and association parameters cannot simply be combined with physical
dipoles. Water and MEA must be refitted with DD active, and the CO2--water and
CO2--MEA mixture treatment must then be reassessed before a polar model can be
judged on predictive data.

The polar EOS rejects exact composition-boundary states. The existing
neutral-pool activity reference is therefore evaluated as an electroneutral
infinite-dilution limit with a declared base ionic mole-fraction floor of
\(10^{-12}\). A sentinel check from \(10^{-8}\) to \(10^{-12}\) changed the
full-M5 pressure by less than 1% at the lowest Hilliard loading; this numerical
effect is negligible relative to the model discrepancy. No reacting state uses
the floor.

M5Q and M5 are diagnostic comparisons only. M5Q changes both the CO2 pure set
and the QQ equation as required for component consistency. M5 does not supply
the missing dipolar-associating refits. Neither result permits parameter
promotion or a manuscript claim.

## R4 temperature-correlation screen

The application-owned R4 screen fits
`ln K4 = A + B/T + C ln(T) + D T` to the admitted Hilliard 40 C and Jou
40--120 C pressure campaigns. Jou is no longer reserved solely for validation:
the grouped split now contains 72 pressure-training rows and 49 reserved rows.
The objective gives each source/temperature group equal total weight and equal
weight to rows within a group. Loading remains a state input, not a target.

Every modeled state uses the installed public pure-water source-reference
transfer and the fixed-T,P homogeneous chemical-equilibrium solver. The
Gauss--Newton step uses exact implicit state sensitivities to the Provider-basis
R4 constant. Four training rows that did not certify in the final-wheel
preflight remain explicit failed predictions and do not enter the fit. They are
listed in `results/r4_state_failures.csv`; no replacement value or relaxed
criterion is used.

The four-coefficient form is evaluated because the five Jou temperatures add
temperature leverage unavailable to the earlier three-row screen. Its scaled
Jacobian rank and condition number are retained with the fitted coefficients.
A numerically improved pressure curve does not make A, B, C, and D separately
identifiable when that matrix remains ill-conditioned.

The single exact-sensitivity update reduces the group-normalized training
log10-pressure RMSE from 1.6584 to 1.3953 and the unweighted reserved RMSE from
2.1777 to 1.5336. The resulting coefficients are
\(A=32032.62\), \(B=-987556.14\,\mathrm{K}\), \(C=-5435.37\), and
\(D=7.49855\,\mathrm{K}^{-1}\). They are not a defensible replacement
correlation: the scaled Jacobian condition number is \(3.43\times10^5\)
(\(1.26\times10^9\) in the raw basis), and the four large coefficients cancel
to produce the fitted values over only five temperatures. No coefficient is
promoted.

The residual evidence also rejects an R4-only explanation. After the update,
the pressure residual has a Pearson correlation of 0.67 with loading and a
slope of 3.77 log10 units per mol CO2/mol MEA. Median absolute pressure
sensitivities to R2, R4, and R5 are 0.382, 0.434, and 0.434 log10 units per
unit change in ln K, respectively; R4 and R5 are therefore effectively
indistinguishable to this pressure observable. On the declared affine scales,
the median pressure responses to MEAH+ dispersion energy, MEAH+ diameter, and
MEACOO- diameter are 5.97, 1.27, and 0.58 log10 units. The source-specific
fitted RMSE values remain 1.84 for Hilliard, 1.25 for Jou, and 0.72 for Xu.
Pressure data alone cannot separate the R4 temperature correlation from R2/R5
chemistry and the selected ionic EOS parameters.

Böttinger carbamate values are interpolated only where temperature,
composition, and loading support that comparison. They never enter the
pressure objective. The retained figures show pressure parity, residuals
against temperature, loading, pressure, and composition, carbamate behavior,
and reaction-versus-EOS sensitivity. All results remain diagnostic because the
R4 and R5 source correlations are extrapolated above 323.15 K, the M5 neutral
parameters were not refitted with the polar terms active, and several states
remain outside a certified Provider/Equilibrium path.

Run in an environment containing the exact candidate EOS and Equilibrium
wheels:

```bash
PYTHONPATH=src:analyses/phase3/m0_m3_model_comparison/scripts \
  python analyses/phase3/m0_m3_model_comparison/scripts/run_r4_correlation_fit.py
PYTHONPATH=src \
  python analyses/phase3/m0_m3_model_comparison/scripts/render_r4_correlation_fit.py
```
