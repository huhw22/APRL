# Gaussian synchrotron benchmark for four longitudinal shapes

This campaign compares `linear`, `quintic-compensated`, `linear-notch`
(7-point), and `linear-binomial3` (3-point) with the same
`conservative-source` particle coupling.

The beam is a product Gauss-Hermite representation of a cold Gaussian bunch:
20 um rms in x/y, 20 nm rms longitudinally, gamma 1174.171072, and 10,000
represented electrons. It traverses a 2.4 mm hard-edge uniform dipole of 3 m
radius. At the selected emission point (the dipole entrance) its tangent is
exactly +z. The screen aperture is centred on the automatically propagated
intersection of that tangent, x=y=0 at z=0.3524 m.

The current conservative-source implementation requires both a moving window
and frozen lab-plane injection. Therefore each shape is simulated as a
strictly matched pair:

1. a zero-field straight control;
2. the finite bend, with every other input and numerical setting identical.

The analysis subtracts the two **complex Fourier E and B fields** before
forming the Poynting spectrum. The independent Lienard-Wiechert reference uses
the same finite record and subtracts the identical incoming straight
worldline. This removes the common frozen-release transient without a fitted
constant, background curve, amplitude normalization, or scalar power
subtraction.

The registered comparison is the signed 8--12 eV energy in the central
20 x 20 um aperture. The theory includes the exact finite hard-edge arc,
finite observation distance and aperture, incoming/outgoing velocity fields,
and Gaussian transverse and longitudinal form factors. The deterministic
weighted PIC source represents the coherent continuum term. The physical
incoherent shot term is reported as a correction estimate but is not inserted
into the simulation or theory curve.

This first campaign is an absolute single-grid comparison. It does not claim
physical accuracy acceptance without a later spatial/time resolution study.
