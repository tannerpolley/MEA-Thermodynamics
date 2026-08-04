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

The retained installed-only audit uses EOS commit
`0fbe76038b5cab22b116f074359447a66ce9da9d` and Equilibrium commit
`e6fd5df6f3d3d156ad09488f3eecc6529cff83a3`; both commits are present on their
respective origins. The exact wheels have SHA-256 values
`1622162b929cb8cd1a10d7c582a6b913babb8580a9f2188ee4fdf324d92f2772` and
`397f0745fc692d33ea3a2d855a33346c516038bfd221798bc05f8ca02fde9b77`.
The installed public EOS performs the declared pure-water source-reference
transfer at the trial temperature and pressure, and the retained result carries
its immutable transfer receipt.

M0 and M5 both report `LOCAL_EQUILIBRIUM`; numerical, physical, first-order,
and second-order checks pass. For M5, the reaction-affinity norm is
\(2.15\times10^{-14}\), the relative pressure residual is
\(5.34\times10^{-11}\), the KKT residual is \(2.84\times10^{-13}\), and the
certified reduced-Hessian inertia is (6, 0, 0). Its packing fraction is 0.4771
and its minimum species amount is \(2.73\times10^{-10}\) mol. The canonical
M5 mole fractions are

```text
[2.276679e-5, 1.943569e-2, 8.755756e-1, 5.269521e-2,
 4.008662e-2, 1.175949e-2, 4.244862e-4, 3.066710e-11,
 1.275447e-7]
```

The search generated 23 starts for each model and evaluated 18 after five
Provider-domain rejections. It found 16 certified M0 starts and 12 certified
M5 starts; two M0 and six M5 evaluated starts did not certify. These results
establish homogeneous fixed-T,P local equilibria only. They do not establish
global equilibrium, vapor-liquid equilibrium, parameter validity, regression
readiness, or predictive accuracy.
