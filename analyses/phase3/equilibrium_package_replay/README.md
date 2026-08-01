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

With EOS `b2638deb64772f2353f0382d0dc5a3210889a827` and Equilibrium
`9b7876d486f10d52b567a3e640ba98beafcb8fbe`, all eight model bundles load, but
the first chemical solve fails before optimization because the installed
Provider does not advertise a neutral-reference callback for this
three-neutral, six-ion system. The retained results therefore classify the
current public route as `BLOCKED_SHARED_PROVIDER_CAPABILITY`; no numerical
parity claim is made.

An additional diagnostic using Provider
`b59c1c7` passed the neutral-reference and inverse-packing interfaces and
reached a numerically and physically valid solution, but Equilibrium did not
certify a local minimum (`primal_solution_not_certified`). This confirms that
the current blocker is not a reason to copy chemistry back downstream. The
next useful work belongs in the generic Provider/Equilibrium contracts. Only
after the value path returns certified states is it worth exercising the
generic evaluator and Jacobian path in `ePC-SAFT-regression`.
