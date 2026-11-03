# Radiation method comparison, 2026-11-03

Controlled comparison of existing explicit, energy-iterated, and longitudinally
paired program settings. No numerical kernel or production default is changed.

The completed campaign, derivations, frozen inputs, figures, and limitations are
in `/public/home/wfollowm/APRL_Workspace/results/radiation_method_validation_20261103/`.
Read `RESULT_REVIEW.md` there first. The confirmed engineering target is about
5 percent; neither the absolute radiation calibration nor the cheaper methods
have received a general physical acceptance.

The prepare/Slurm launch files preserve campaign inputs and command provenance.
Do not rerun them merely to regenerate a report. The original FEL batch exited
nonzero during its old full-record uniform-time assertion, after all four PIC
processes had exited zero. Its frozen evidence is preserved. The repaired
analysis is in each case's `analysis-regularized-r1/`, using identical retained
native time intervals without missing-sample interpolation.

Closed-data report regeneration uses the physics Python environment:

```bash
python build_report.py /public/home/wfollowm/APRL_Workspace/results/radiation_method_validation_20261103/fel-r4
python validate_campaign.py
```

`check_retarded.py` and `check_lw_normalization.py` independently check the
theoretical reference, not the PIC solution. `model_scope.py` estimates the
missing transverse-space-charge envelope effect without fitting the PIC.
`spatial_audit.py` and the fixed-reference-mode analyses diagnose spatial
contamination; they are not general near-field subtraction tools.

All three candidate methods use the same existing binary and physical input.
The no-feedback case is a diagnostic control and is not a candidate coupled
method. Reported costs include initialization, strict ledgers, and output on
the same one-node allocation. They do not estimate diagnostic-off production
cost or certify SASE shot noise.
