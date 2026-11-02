# Standard Yee lab energy/momentum diagnostic

This verification case reads the immutable production field library. It changes
no Maxwell, particle, current-deposition, reference-supply or CPML kernel.
The target is the -1.053453% residual of the completed oblique vacuum pulse in
`results/lab_energy_balance_20261102/vacuum_propagation_r4/yee-oblique`.

Prepare a new output directory:

```
python3 prepare.py --workspace /public/home/wfollowm/APRL_Workspace --output /public/home/wfollowm/APRL_Workspace/results/lab_energy_balance_20261102/yee_sbp_r1
```

The generated job reserves one exclusive 128-CPU node, uses 16 MPI x 8 threads,
checks physical core binding, and runs only the previous Yee oblique pulse.
All inputs/programs are hashed before submission. Results and frozen sources
stay separate from these versioned templates. Existing results are never reused
as writable work directories. The original 600/60-second MPI-free teardown
watchdog and allocation-local self,vader,tcp transport are retained.

Let E_n, B_minus, B_plus be actual stored fields, B_bar=(B_minus+B_plus)/2.
With native quadrature weights, define

```
U_sync = eps/2 <E_n,E_n> + 1/(2mu) <B_bar,B_bar>
H_Yee  = eps/2 <E_n,E_n> + 1/(2mu) <B_minus,B_plus>
V_B    = 1/(8mu) <B_plus-B_minus,B_plus-B_minus>
U_sync = H_Yee + V_B
```

H is the time-paired discrete energy associated with these Yee equations, not
a newly supplied physical source. The lab time-translation charge is formed
as Gamma(H + beta cP_z), with the same native paired momentum geometry used by
the independent three-dimensional summation-by-parts proof.

For the transition n -> n+1, the independently evaluated surface terms are:

```
Q_H = dt * native surface quadrature of (E_n+E_new)/2 x B_plus / mu
Q_P = [Q_E(E_n)+Q_E(E_new)]/2 + Q_B(B_plus)
```

The electric stress must average the quadratic expressions in time. In general
Q_E((E_n+E_new)/2) differs by -Q_E(E_new-E_n)/4. Native energy/stress ports also
multiply the matching staggered components before transverse quadrature; they
are not all products of six fields averaged onto a single face centre.

Gauss E and B defects are measured independently. The observed old residual is
compared with the sum of magnetic stock time variance, energy-flux spatial
pairing, stress spatial pairing, stress time pairing, and the matched residual.
No term is inferred from the old residual. Omitting a measured term must fail
the reconstruction. Canonical closure does not establish continuum accuracy,
validate Cowan, separate radiation, or explain the full 19174 particle gain.

`check_geometry.cpp` checks the surface formulas and finite-step energy and
momentum identities on 64 arbitrary three-dimensional fields before the wave
run. Those are program/math checks, kept distinct from actual wave evidence.
