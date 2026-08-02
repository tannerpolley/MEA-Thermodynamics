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
`vle_obs_0137` in the existing comparison. The Equilibrium lane deliberately
uses the newer source-adjudicated reaction constants and exact Provider
neutral-reference transformation. Exact agreement with the old lane is
therefore not presumed: a successful run would separate solver replacement
effects from the known reaction/reference-state correction.

## Run

Build or install exact `epcsaft` and `epcsaft-equilibrium` wheels, then run:

```bash
PYTHONPATH=src:analyses/phase3/m0_m3_model_comparison/scripts \
python analyses/phase3/equilibrium_package_replay/scripts/run_replay.py
```

The script records installed package identities, source commits, model
fingerprints, elapsed time, typed failure diagnostics, and—when a state is
certified—the Equilibrium composition and liquid CO2 fugacity. It fails closed
if an upstream callback or certification gate rejects the state.

The default run tests all eight models. Use `--models M0` for a short
capability probe before committing to the expensive multi-start campaign.

## Current result

The retained campaign uses EOS
`02104702e822e1f062bf829f0fe2280e801bbbc4` and Equilibrium
`7f60cbc3619b2036ea4fbe0f5a3109d63703410c`. The latter is an exact local
feature commit, not a merged or released dependency. Six of the eight models
return certified homogeneous local equilibria:

| Model | Result | Liquid CO2 fugacity (Pa) | log10(new/retained) |
| --- | --- | ---: | ---: |
| M0 | certified local equilibrium | 2.28693e-4 | -7.055 |
| M1 | certified local equilibrium | 3.14608e-4 | -6.753 |
| M2 | physical domain not admitted | -- | -- |
| M3 | certified local equilibrium | 2.42687e-5 | -7.374 |
| M4A | certified local equilibrium | 2.80879e-5 | -7.575 |
| M4B | certified local equilibrium | 5.46497e-4 | -6.318 |
| M5Q | certified local equilibrium | 4.15670e-4 | -7.159 |
| M5 | multistart search exhausted without certification | -- | -- |

M2 is rejected during Provider start-pressure bisection before a candidate
state is available. M5 evaluates 21 starts under a declared budget of 25, but
no candidate passes all first- and second-order gates. Its retained terminal
state satisfies material balance and Provider-domain checks but fails physical
KKT stationarity, reaction affinity, and pressure-residual gates.

The certified values differ from the older MEA-owned solver by six to eight
orders of magnitude because this lane also applies the source-adjudicated
reaction constants and exact neutral-reference transformation. The comparison
therefore does not establish numerical parity or predictive accuracy. Each
successful row establishes one local homogeneous fixed-temperature,
fixed-pressure equilibrium only; global equilibrium, parameter validity, and
regression readiness remain outside this analysis.
