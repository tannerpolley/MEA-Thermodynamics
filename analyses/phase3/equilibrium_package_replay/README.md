# Public Equilibrium-package replay

This standalone analysis asks whether the public `epcsaft_equilibrium`
fixed-temperature, fixed-pressure homogeneous reacting-phase solver can replace
the MEA-owned SciPy speciation solver used in the M0--M5 comparison. It uses
the same eight parameter variants and the same frozen tracer state:

- 313.15 K;
- 7326.7 Pa total liquid-state pressure;
- 30 wt% unloaded MEA;
- 0.466 mol CO2 per mol MEA;
- nine species and the five source-bound reactions in the canonical reaction
  contract.

The reference value for each model is the retained prediction for
`vle_obs_0137` in the existing comparison. The literature constants use pure
water as the solvent reference, with every nonwater solute at infinite
dilution in water on the common aqueous-molality scale. The Provider's
salt-free equimolar CO2/MEA/water neutral reference is a different
computational gauge. An exact, charge-neutral transfer between those
references is therefore required before the reacting-phase solver may run.

## Run

Build or install exact `epcsaft` and `epcsaft-equilibrium` wheels, then run:

```bash
PYTHONPATH=src:analyses/phase3/m0_m3_model_comparison/scripts \
python analyses/phase3/equilibrium_package_replay/scripts/run_replay.py \
  --models M0 M5 \
  --provider-wheel /path/to/epcsaft.whl \
  --equilibrium-wheel /path/to/epcsaft_equilibrium.whl \
  --provider-wheel-sha256 <sha256> \
  --equilibrium-wheel-sha256 <sha256>
```

Before evaluating chemistry, the script matches the installed distributions to
the exact wheel RECORD and public-header bytes, confirms that neither package
was imported from a source checkout, and inspects the installed public API for
the required transfer. If that capability is absent, it writes
`BLOCKED_UNSUPPORTED_SOURCE_REFERENCE_TRANSFER` and does not construct a model
or call the solver.

The default run tests all eight models. Use `--models M0` for a short
capability probe before committing to the expensive multi-start campaign.

## Current result

The retained audit uses EOS commit
`02104702e822e1f062bf829f0fe2280e801bbbc4` and Equilibrium candidate commit
`6604555a4b0c4efb733281bfa00c7f5efdefd772`. The Equilibrium commit is local and
is not contained in an `origin/*` branch. Its continuation and certification
interfaces do not add a caller-declared source-solvent reference. The installed
EOS callback also exposes its fixed computational reference without a general
public source-reference transfer operation.

Both M0 and M5 therefore report
`BLOCKED_UNSUPPORTED_SOURCE_REFERENCE_TRANSFER`. The solver is not called, so
there are no reaction affinities, pressure or KKT residuals, reduced-Hessian
inertia, minimum amount, packing fraction, or start accounting to report. The
previous M0 certificate and M5 diagnostic composition remain diagnostic
evidence only; neither is a corrected source-bound M5 result.
