# SPARC completed-output bandwidth scan

This analysis scans nominal resonance frequency half-widths of 1%, 2%, 5%, 10%,
20%, and 50% on all eleven completed detector planes. It reads each physical
28-by-28 aperture once and reuses each Fourier transform across widths.

The primary time window uses the later half of the 1694 source slices. Within
that half, the known parabolic profile selects I/I_peak >= 0.8; 10% of this
selected interval is removed at each edge. Thresholds 0.6, 0.7, 0.8 and 0.9,
and the prior full-half edge guards, are retained in a sensitivity CSV.
The window is independent of field amplitude and fixed across frequency widths.

Outputs distinguish full-record Hann energy, unwindowed filtered power envelopes,
and a separate Hann spectrum of the central time interval. Spectral energies
use the established RMS-normalized Hann window. Broader integrated positive
spectral energy is monotonic; power within a finite time interval can change
nonmonotonically because of interference between admitted frequencies.

Run the low-cost Fourier normalization check with:

```bash
/public/home/wfollowm/micromamba/envs/physics/bin/python3 analyze_bandwidth_scan.py --self-test
```

The compute analysis must run on Slurm. Pass a closed campaign, a new output
directory and this template directory to analyze.sbatch. Existing output
directories are refused. A successful analysis writes COMPLETE.json, verifies
Parseval and the existing eleven-plane 5% results, and saves reduced spectra
and envelope NPZ files to allow later time-window statistics without H5 rereads.

The first live campaign is results/sparc_sase_11segment_2slip_20261111_r1,
job 19735. No running momentum-matched field H5 is accessed. This
post-processing does not alter the PIC configuration, beam, or binary.
