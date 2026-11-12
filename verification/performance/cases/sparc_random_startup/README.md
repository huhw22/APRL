# Independent-six-dimensional SPARC startup control

This restores the accepted independent physical-electron sampler, and uniformly
thins its fixed-N conceptual source. Only sampled records are stored. There are
no beamlet copies, phase corrections, empirical moment fits or added Fawley
noise. Source IDs and all six original phase-space coordinates are retained.
One w512 sample is nested in the w64 sample; a second physical seed supplies
another w512 case. All samples have the same total charge.

The full first77-period undulator, its one-period end ramps,151.9MeV parameters,
56x56 grid at the same40um spacing,CPML6,20 longitudinal cells per boosted wavelength,EC+F3,80 requested
particle samples/period, reference supply and600/60s MPI-free watchdog remain.
The source is the central231 optical periods of the current full parabola;
the fixed central77 generation bins are the primary PIC cohort. This finite
source crop is shared with GENESIS and explicitly retained as a limitation.

GENESIS receives exactly those events with one4one=false,nbins=1,shotnoise=false
through native importbeam. Variable slice populations and weighted currents are
retained. Import/write roundtrips must preserve every active 6D/current array
bitwise before PIC can run. This avoids one4one sorting and phase resetting;
no quiet loader or transverse cloning is involved. Both SC-off and longrange
controls use the previously validated ADI-only corrected executable.

prepare.py NEW_DIRECTORY freezes source code/binaries. prepare.sbatch generates
and validates inputs on one exclusive128-core node; main.sbatch is a dependent
three-task array, each16MPIx8 on one exclusive node. No nodes are pinned, no
automatic retry is used, and incomplete outputs are never analyzed.

The raw BF, common-phase BF projection, fixed-core PIC populations and±5%
Rect/Hann PIC band energy are retained at7/19/37/57/77 periods. GENESIS fields
use the same±560um aperture and±5% Rect band. Eulerian GENESIS local statistics
are labeled separately from fixed-ID PIC cohorts. Further matched-statistic
review is required before interpreting local-growth differences. No power or
cross-weight amplitude matching is applied. Raw PIC E/B include bound fields.

These are weighted independent samples, not unit-electron noise: entrance BF
noise follows sum(w²)/sum(w)². A first-segment comparison of two physical seeds
and two weights cannot establish ensemble SASE growth or eleven-segment
saturation. Result files keep physical_acceptance=false until reviewed.

The initial r1 attempt used40x40 and stored the physical macro-weight as the
input relative weight. Its initialization checks rejected a few Gaussian tails
inside CPML and a noninteger weight-summation error at the higher particle count;
no physical stepping occurred. The replacement uses56x56 to retain every random
event, and unit relative weights plus the unchanged physical total electron count.
GENESIS current uses the actual physical macro-weight. Spectral comparison still
uses±560um; PIC statistical charge uses actual charge/e, not raw relative weights.
No clipping, reweighting of individual particles, moment fit or core change occurs.
