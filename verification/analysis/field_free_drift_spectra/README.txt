Existing full-run field spectra: undulator at 4.5 m vs free drift at 0.5 m

Scope
Read the existing lab E/B output from production jobs 18932 and 18938. No PIC,
near-field subtraction, mean removal, notch filtering or time resampling.
The original input files are read-only. Full scans run on Slurm compute nodes;
plots and numerical review read only the small reduced NPZ/CSV/JSON files.

Sources
APRL_Workspace/results/elegant_field_only_20261027/production (18932, 4.5 m)
APRL_Workspace/results/elegant_free_drift_20261027/production (18938, 0.5 m)
Both use build/detector-sampling-20261027/aprl, with the corrected MPI detector
interpolation. These are not the historical 18893/20261021 spectrum sources.

Method
Use the accepted ../retirement_field/analyze_spectra.py kernel, frozen with
this wrapper in the result directory. For every physical detector pixel,
compute F_k = dt * sum_n w_n * field_n * exp(-2 pi i k n/N).
Then sum the squared amplitudes or E/B cross products over pixels. Never FFT
a spatially averaged field to estimate full-plane energy.

T=N*dt; epsilon=h*f in eV, h=4.135667696e-15 eV s; d_epsilon=h/T.
s_k=2 on interior positive-frequency bins and 1 on DC/even-N Nyquist.
Field PSD D_Ej = s_k * sum_pixels |F_Ej|^2 / (N_pixels*T*h).
It integrates over epsilon to the spatial/time mean square of that component.
Magnetic fields are multiplied by c before this calculation for comparable
electric-field units. PSD is not a radiation-energy density.

Native whole-plane energy spectrum:
dW/d_epsilon = s_k*dA/(mu0*h) * sum_pixels Re(F_Ex*conj(F_By)-F_Ey*conj(F_Bx)).
It integrates to dt*sum_samples,pixels[(Ex By-Ey Bx)*dA/mu0], retaining sign.
Use the same spatial pixel for each cross product, without cross-pixel terms.
Native fields inherit different Yee time layers and their Lorentz mixing;
this is the spectrum of the stored native product, not certified simultaneous
Poynting flux or a separated pure-radiation spectrum. A lab-B-only phase shift
is not applied. Neither E^2 alone nor a two-run difference establishes radiation.

The rectangular window is primary. The auxiliary Hann window is divided by
sqrt(mean(w^2)) on all six components. This normalizes mean-square window gain;
for a nonstationary pulse it need not reproduce the rectangular total energy.
Both separately pass Parseval. No curve is scaled to match another run.
Both runs retain their full record and native frequency bins. Their T differs;
mean-square PSD divides by the corresponding T, while dW/d_epsilon does not.
Band powers, if shown, divide band energy by that run's T and are record-average
powers, not instantaneous or repetition-rate averaged source power.

Physical aperture: remove the configured 10 CPML cells on each transverse edge,
leaving 140x140 pixels at 5 um spacing, area 0.49 mm^2. Check complete=1,
committed samples == field/time extents, lab/raw/corrected-MPI metadata, every
field scalar finite, cadence uniformity, Parseval, and unchanged file stat.
18932 batch exit124 is retained explicitly even when its closed field passes;
field validity does not certify a successful MPI teardown of the whole job.

Interpretation limits
The two runs share input/held/boost/resolution/core. Turning off the undulator
also changes the automatic longitudinal extent and natural Boris substeps.
Surviving detector cohorts differ because of transverse CPML exits. Free drift
still contains the original input modulation, startup and self-consistent
motion; it is a background comparison, not a certified zero-radiation solution.
No quiet post-pulse tail or continuum/grid convergence is established here.

Execution
Prepare a fresh result directory with frozen scripts/analyze_spectra.py,
scripts/reduce.py, scripts/run.sbatch and scripts/SHA256SUMS. Submit run.sbatch
with that directory and workspace as positional arguments. It reserves one
exclusive 128-CPU node and uses 8 bound CPUs, 128 GiB memory and 35 minutes.
The binding/self-test precedes both reductions; each reduction has 900 s cap.
Nodes are not pinned. The plot.py script consumes the reductions afterwards.
