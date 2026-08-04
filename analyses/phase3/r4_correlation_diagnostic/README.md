# Multisource R4 correlation diagnostic

This analysis tests whether refitting the four coefficients in

\[
\ln K_4(T)=A+B/T+C\ln(T)+DT
\]

can explain the admitted MEA--H2O--CO2 pressure data when the full-polar M5
EOS model and the public homogeneous chemical-equilibrium solver are held
fixed. It is excluded from the manuscript and does not promote parameters.

The frozen split contains 72 pressure-training rows and 49 reserved rows.
Training uses the Hilliard 313.15 K campaign and the Jou 313.15--393.15 K
campaigns with equal total weight per source/temperature group. Loading is a
state input, not a regression target. Four training states and one reserved
state fail the installed Provider/Equilibrium acceptance path and remain in
`results/r4_state_failures.csv` rather than receiving replacement values.

One exact-sensitivity, column-scaled Gauss--Newton update reduces the
group-normalized training log10-pressure RMSE from 1.6584 to 1.3953 and the
reserved RMSE from 2.1777 to 1.5336. The fitted coefficients are numerically
unstable: the scaled Jacobian condition number is 3.43e5 and the raw condition
number is 1.26e9. Residuals remain strongly loading-dependent, while R4 and R5
have nearly indistinguishable pressure sensitivities. Selected ionic EOS
coordinates also move pressure by comparable or larger amounts. These results
reject an R4-only explanation; they do not establish a new reaction
correlation.

The retained evidence is deliberately small:

- `results/r4_correlation_fit_rows.csv`: measured and modeled pressure rows;
- `results/r4_correlation_fit_parameters.csv`: literature and fitted values;
- `results/r4_correlation_fit_metrics.csv`: grouped residual metrics;
- `results/r4_reaction_sensitivity_rows.csv` and
  `results/r4_eos_sensitivity_rows.csv`: confounding evidence;
- `results/r4_state_failures.csv`: rejected states;
- `results/r4_correlation_fit_receipt.json`: compact run and identity summary;
- three PNG figures for pressure parity, residual structure, and sensitivity.

Run in an isolated environment containing the exact candidate EOS and
Equilibrium wheels:

```bash
PYTHONPATH=src python analyses/phase3/r4_correlation_diagnostic/scripts/run_analysis.py
PYTHONPATH=src python analyses/phase3/r4_correlation_diagnostic/scripts/render_figures.py
```

The full numerical campaign is expensive and is not part of the default test
suite. Repository tests protect source partitions and contracts; this command
owns the integrated scientific reproduction.
