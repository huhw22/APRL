# Retirement detector field plots

Read an existing laboratory field-plane file after particle retirement. This
workflow does not rerun APRL, subtract a ballistic reference, or assume that the
remaining Maxwell field has been proved to be pure radiation.

`../elegant_4m/analyze_fields.py --field-only --field-file FILE --run-dir RUN
--output-dir OUTPUT --no-plots` streams only the committed prefix, requires a
clean HDF5 completion marker, checks every field value for finiteness, and saves
small NPZ/CSV/JSON reductions. `RUN/run.yaml` must belong to `FILE`: it supplies
the transverse CPML thickness used for the physical aperture. Verify their
identity before running. The source HDF5 and input card are read only.

Use `run.sbatch` for a large scan. Copy the reducer and this batch script into
a fresh result directory's `scripts/`, and create `scripts/SHA256SUMS` relative
to that result directory. The batch arguments are the production run directory,
fresh result directory, and detector HDF5 file, all absolute paths. It reserves
one exclusive 128-CPU node under the workspace policy; the sequential scan uses
one bound CPU and one library thread, with a 20-minute wall limit. Binding and
the reducer's analytic plane-wave check run before the scan. The saved
allocation and binding log distinguish reservation from actual CPU use.

`plot_fields.py` renders the small derived files as PNG and PDF. Read its command
line help for output arguments. Rendering can be repeated without rescanning
the detector file. Use `--sampling-note TEXT` to carry a known diagnostic issue
into every figure caption and the plot audit, without altering any data.

Definitions:

- Coordinates and fields are laboratory-frame SI values; displayed distances
  are micrometres and electric fields are MV/m.
- A peak map is `max_t sqrt(Ex² + Ey² + Ez²)` at each pixel. Its pixels need
  not peak at the same time. A time-RMS map uses trapezoidal integration over
  the whole recorded time interval, including the quiet portions.
- A spatial-RMS trace is `sqrt(mean_xy(Ex² + Ey² + Ez²))` within the physical
  aperture at each time. Component traces use the same spatial averaging.
- The central waveform is the vector average of the four nearest cell centres
  for the even grid, not a direct sample at the geometric origin.
- Time zero is the first recorded detector sample. The absolute value is
  retained in the summary and timeseries files.
- The stored peak frame is selected by the **full recorded aperture**. Check
  its location before describing it as the physical aperture's peak frame.
- Native signed flux is `(Ex By - Ey Bx)/mu0`. The saved B has Yee half-step
  staggering and the laboratory transformation mixes staggered components;
  the derived flux is not an exact simultaneous-field radiation measurement.
- Finite aperture and a complete file do not prove complete signal capture.
  Neither a transverse field pattern nor `E ≈ cB` alone proves near-field
  removal.

The 18893 analysis is in the workspace's `results/field_18893_20261027/`.
Its separate `DATA_STATUS_REVIEW.json` explains the post-run count-check failure
and distinguishes that failure from detector-file completeness.

## 18893 interpretation review

The complete 34,215-frame scan passed finite-value and reduction checks. It
also exposed a detector interpolation defect: all 160 prominent narrow peaks
occur at sample index `14 mod 182`, with 182 cells per MPI slab. Geometry and
source review confirm that these are the first samples after plane ownership
changes. Local half-z-staggered components clamp at the local last cell rather
than including the needed neighbouring rank's value. This path exists in both
production core `ac58e67` and the checkpoint branch's `deb212e`; this analysis
does not modify it. See the result directory's `SAMPLING_REVIEW.txt` and
`SPIKE_PERIOD_REVIEW.json` for exact geometry, source locations and evidence.

The raw peak is 2.164706 GV/m and native peak power is 21.163494 MW, including
these suspect samples. They are not accepted physical radiation peaks. The
154.07 fs recorded-window native integral is 664.47 nJ; its bias has not been
quantified. Signal remains strong at the end (8.826 MW native power), and that
last sample does not trigger the identified boundary clamp. Neither the slow
offset nor the terminal signal has been attributed solely to the narrow-peak
defect. Do not infer pure radiation or complete pulse coverage from these plots.

## Raw time spectra

`analyze_spectra.py` reads a committed, complete raw laboratory field file once,
then transforms every physical-aperture pixel in time. It averages squared
Fourier amplitudes across the aperture, never the electric vectors before the
transform. Electric and `cB` component spectra have units `(V/m)^2/eV` and
integrate to the discrete space/time mean square. They are not pure-radiation
energy spectra. The separately saved signed native flux spectrum has units
`J/eV` and inherits the detector's known spatial and temporal sampling issues.

The principal rectangular-window result retains DC, the native signal, all
spikes and all frequencies; no mean removal, detrending, background subtraction
or angular selection is applied. An additional Hann window is divided by
`sqrt(mean(window**2))` on all six components. Its different temporal weighting
tests endpoint leakage; it cannot recover missing pulse history. The FFT period
is `N*dt`, with correct odd/even one-sided factors, and both windows must pass
component and signed-flux Parseval checks. Native timestamps are accepted only
when a uniform grid fits to within 1e-4 sample; they are not interpolated.

`spectra.sbatch` scans the preserved original ELEGANT field and the 18893 field
sequentially. It reserves one exclusive 128-CPU node and 128 GiB for up to
45 minutes, using eight bound CPU cores for SciPy FFT and one BLAS thread. The
physical field cube is about 30 GiB per run; only one run is resident at a time.
Place frozen `analyze_spectra.py` and `spectra.sbatch` in a new result directory's
`scripts/`, generate `scripts/SPECTRA_SHA256SUMS` relative to that directory,
then submit with the absolute result-directory and workspace arguments. A
binding check and odd/even analytic-spectrum self-test precede the scans.

`plot_spectra.py` renders the small results. The two configurations differ in
detector position, startup, initial-field treatment and longitudinal grid;
spectral differences are not an isolated absorption comparison. The theoretical
sampling ceiling is not a physical convergence guarantee. The known MPI
interface artifact produces different comb spacings for 182 versus 178 cells
per rank; these are marked rather than removed.

The two-run analysis is in `results/field_spectra_20261027/`, including a
metadata-only source audit and a separately labelled central-waveform preview.
