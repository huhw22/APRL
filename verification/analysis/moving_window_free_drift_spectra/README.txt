Completed 0.5 m free-drift spectrum: moving window vs earlier fixed window

Inputs
New: results/moving_window_comparison_20261029/free_drift_0p5m (18993).
Old: results/field_free_drift_spectra_20261028/free_drift (18938 saved spectra).
The earlier raw field remains in results/elegant_free_drift_20261027/production.
Read neither the running 4.5 m output nor any deleted absorption-run data.

Preparation
Use a fresh result directory containing scripts/ and old_free_drift/.
Copy these reduce.py/run.sbatch/plot.py/review.py into scripts/.
Freeze ../retirement_field/analyze_spectra.py as scripts/analyze_spectra.py;
freeze ../field_free_drift_spectra/plot.py as scripts/shared_plot.py and
../field_free_drift_spectra/independent_review.py under its original basename.
Copy the old saved summary.json, spectra.npz and spectra.csv into old_free_drift/.
Before reuse, check the old raw stat, summary hash and original source audit.
Check the new complete/committed transaction and successful physical stop.
Record source hashes and the comparison of lab time samples and aperture.
Generate scripts/SHA256SUMS with paths relative to the new result root.

Execution
Submit run.sbatch with the prepared result directory and workspace arguments.
It allocates one exclusive 128-CPU node, actually binds 8 CPUs for FFT workers,
uses 128 GiB memory, and has a 20-minute hard limit with a 1000-second reduction
limit. There is no PIC simulation. The large new field is read only on Slurm.
An odd/even FFT test and binding check precede the read.

After reduction, plot.py --input-dir RESULT and review.py --root RESULT read
only small NPZ/JSON products; they do not reopen the large fields. Inspect PNG
and PDF exports and check plotted band metrics against INDEPENDENT_REVIEW.json.

Analysis definition
Unchanged accepted per-pixel time FFT, then full physical-plane area sum.
140x140 physical pixels, 5 um pitch, 0.49 mm^2; transverse CPML edges excluded.
Keep DC in every integral. Rectangular window is primary; Hann/sqrt(mean(w^2))
is an alternate record weighting, not an uncertainty interval or tail recovery.
The native signed E/B cross-spectrum is already area-integrated in J/eV;
integrate once over energy bins, then divide once by N*dt for record-average W.
Do not multiply by area twice, use E^2 alone as radiation, or FFT a plane average.
The E/cB PSD is a pixel/time mean-square density and has a different meaning.
All half-open band sums use lower <= bin centre < upper. Native bins remain.

The fixed-slab comb-width parameter is not meaningful for the moving window.
reduce.py uses a legacy placeholder only within the unchanged FFT producer,
then stores slab_cells=null and sampling_comb_spacing_eV=null. This argument
has no effect on data, normalization, Fourier transforms, or bin integration.
No peak, DC, trend, background or time-layer correction is applied.

For these two complete runs all 30019 lab timestamps are bitwise equal, so
the full record is also the common record. There is no fitted alignment.
Changing initial reference and moving boundaries remains a physical method
comparison, not bitwise replay. Stored lab E/B mix native Yee time layers;
the resulting total-field flux is not separated pure radiation. Nyquist is
only a sampling bound, not a physical-resolution certification.

Versioned scripts affect postprocessing only. They do not alter simulation
core, cards, binaries, defaults, or submit any additional production runs.
