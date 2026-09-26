# Cowan-z Maxwell kernel and CPML boundary

## Controlled-dispersion update

The primary solver keeps the ordinary staggered Yee E/B storage and changes
only Faraday's update:

```text
Dt B = -curl'(E)
Dt E =  c^2 curl(B) - J/epsilon0
```

For a derivative in direction `i`, `curl'` uses `D_i S_i`. `S_i` is a
nine-point smoothing stencil in the two directions transverse to `i`. The
smoothed values are evaluated in place and are never retained as three
full-domain field copies.

Let `lambda` be no greater than any grid spacing and define

```text
r_i = lambda^2 / delta_i^2
D   = r_x r_y + r_y r_z + r_z r_x .
```

The implemented Cowan coefficients are

```text
beta_i     = r_i/8 * (1 - r_x r_y r_z / D)
delta_ij   = r_i r_j * (1/16 - r_i r_j/(8 D))
alpha_i    = 1 - 2 beta_j - 2 beta_k - 4 delta_jk .
```

The last relation makes every smoothing stencil sum to one, so a constant
field is neither damped nor amplified. With these coefficients the vacuum
dispersion relation factors into a bounded product. Analytic stability
requires `c*dt <= lambda` and `lambda <= delta_i` for every axis.

For the z-priority solver, `lambda=dz`. Therefore `dx>=dz` and `dy>=dz` are
necessary together with the Cowan coefficient formula and `c*dt<=dz`; the two
mesh inequalities alone are not a general stability proof. The code chooses
`c*dt=dz`, which additionally makes resolved vacuum propagation exactly
dispersion-free along z. This is the analytic stability boundary, so
long-time, particle-coupled and boundary-coupled regression remains mandatory.

For a cubic grid all `r_i=1`. Substitution gives

```text
beta  = (1/8)(1-1/3) = 1/12
delta = 1/16 - 1/24  = 1/48
alpha = 1 - 4 beta - 4 delta = 7/12 .
```

These three values are only the cubic special case. An anisotropic production
mesh obtains different coefficients, printed in the startup log.

Reference: B. M. Cowan et al., [Generalized algorithm for control of numerical
dispersion in explicit time-domain electromagnetic
simulations](https://doi.org/10.1103/PhysRevSTAB.16.041303), 2013.

## CFS-CPML coupling

The boundary is an unsplit complex-frequency-shifted convolutional PML. For a
derivative in direction `i`, the stretched derivative is evaluated as

```text
D_i -> D_i/kappa_i + psi_i
psi_i(n) = b_i psi_i(n-1) + c_i D_i(n) .
```

For Cowan Faraday updates, `D_i(n)` means the complete modified difference
`D_i S_i(E)`. Ampere updates use the ordinary staggered `D_i(B)`. Thus the PML
does not silently change the boundary cells back to a Yee magnetic stencil.

The conductivity, kappa and optional frequency-shift profiles are polynomially graded.
The requested `target_reflection` determines the continuous-profile maximum
conductivity; it is a design parameter rather than a guarantee of the measured
discrete reflection for every wavelength, angle and layer count.

The supplied FEL examples use zero frequency shift (`alpha_fraction=0`). In
the current near-axis forward-pulse regression this reflected substantially
less energy than the initially tested nonzero shift. Nonzero alpha remains
available for targeted low-frequency or evanescent-field studies, but is not
assumed to be a universal improvement.

There are two convolution histories for each field component, twelve in
total, but each history is allocated only where its derivative coordinate lies
inside a PML slab. No disabled axis and no non-PML interior point receives
history storage. The outermost surface is terminated by PEC.

Version 1 limitations:

- TF/SF incident seed waves are rejected while CPML is active. Supporting this
  combination requires a closed injection surface and complete correction
  terms inside the CPML interior; these are not implemented in version 1.
- Each z CPML slab must fit on its endpoint MPI rank. The program checks this
  before particle input and reports how to change ranks, z cells or thickness.
- CPML absorbs fields. Particle escape uses the CPML-aware carrier treatment
  specified in [PARTICLE_OPEN_BOUNDARY.md](PARTICLE_OPEN_BOUNDARY.md); it does
  not introduce a physical absorbing material model.

Reference: J. A. Roden and S. D. Gedney, [Convolution PML: An efficient FDTD
implementation of the CFS-PML for arbitrary
media](https://doi.org/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A),
2000.
