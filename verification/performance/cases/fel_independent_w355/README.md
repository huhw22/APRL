# Independent w355 FEL campaign

User request: stop eleven-segment w64 Puffin job19775; submit six-segment
Puffin/APRL with identical independent particles and magnets, and restart the
same four-slip eleven-segment Puffin case at w355. APRL retains the previously
requested EC+three-point and EC+momentum-matching methods in a serial array.

`prepare.py OUTPUT` freezes static plans and copies the validated Philox fixed-N
generator and existing software. No C++/Fortran changes, builds, default switch,
old input overwrite, phase clones, Fawley loading, moment fitting or extra noise.
Generate each full conceptual mother source once, globally sample IDs uniformly
without replacement, and store only those independent six-dimensional events.
The physical count is round(integrated charge/e); equal weight is N/Nmacro,
approximately355. Stored APRL relative weights are all1.

The six-segment source keeps the original2.7762252654ps full parabola and peak
current53.275A (about98.6031pC). Eleven segments keep the cancelled run's
5.5524505308ps four-slip parabola (about197.2062pC),81-square transverse grid,
40um cells,20 longitudinal intervals/wavelength,80 particle steps/period,
diffraction,0.3-reference-frequency filter,periodic redistribution and outputs.
Only sample count/weight and resource time limit change in eleven-segment Puffin.

Puffin uses the existing qualified continuous`aprl-sparc11`external-field route:
original one-period vector-potential entrance/exit ramp and derivative,off-axis
fields,0.4m gaps,finite quadrupoles and original segment phases. Six segments stop
before the seventh magnet,with no intermediate`correctTrans`. Both six-segment
programs use the same original eleven-element magnetic list. Puffin entrance
events are the exact ballistic1cm transport of the shared APRL z=0 events;
Puffin omits the PIC self-field in that drift. This distinction remains explicit.

Submit`prepare.sbatch OUTPUT`,then each production card with`afterok:PREP_JOB`.
Preparation checks the actual native seven-column Puffin imports bitwise and
their first smooth-entrance steps,zero transverse duplicates and6D/bin moments.
It also runs40-step full-input APRL checks of both methods and bitwise observer
off/on invariance. Those deliberate duration-cutoff checks are marked physically
incomplete; only their diagnostic gates pass. Bad preparation blocks production.

APRL uses the previously validated full-ledger binary,16MPIx8,56-square/40um/
CPML6 to retain Gaussian tails without clipping. Existing boost discrete,lab
momentum/stress,whole-Cowan,worldtube and particle side/end transport remain on.
Its worldtube faces are now+/-840um;common-aperture spectral comparisons may
still use+/-520um. Do not mislabel the new worldtube as the old520um aperture.
Observation padding adds128 optical cycles and the reference mesh is recomputed
from the actual full input. This padding changes no bunch charge or beam length.
APRL detector planes use Puffin's actual integration-step exit coordinates.

All jobs use normal scheduling,no pinned nodes,no requeue or automatic retry,
the existing self/vader/tcp network and600/60s shutdown/watchdog rules. Puffin
six/eleven use one128-rank/two256-rank allocation;APRL serial array members each
use one128-core allocation. Physical hard limits are22h/70h/70h respectively;
the new eleven-segment allocation is72h to accommodate periodic overhead.

After submission and available startup checks,stop proactive monitoring until
the user wakes the task. Never read running H5. Retain all raw fields,particles
and diagnostics; closed-run analysis waits for the user. Neither successful
preparation nor normal runtime completion implies physical acceptance.

`compare_closed_six.py CAMPAIGN OUTPUT [PREVIOUS_REVIEW]` compares the closed
independent six-segment Puffin, EC+F3 and EC+momentum-matching outputs. It reuses
the established central Eulerian
window, RMS-Hann +/-5% physical band, +/-520um aperture and local10-lambda RMS
BF definitions from `puffin_six_same_ledger/growth_compare.py`; its source
qualification and plots explicitly refer to the new independent input.
Archive that helper and `puffin_eleven_four_slip/common.py` beside the analysis
script when running a frozen copy. It never opens eleven-segment scientific
outputs. The optional previous two-method review is accepted only after its
artifact hashes, original output identities, input, window and source parameters
match. Its twelve baseline rows are reused exactly; only the new momentum
particle and field statistics are calculated. The previous review stays intact.

For the completed19796_0/19796_1 cases the physical solvers returned0 at the
configured stop with28199 steps, while their frozen driver expected28198 and asserted.
The comparison requires that exact known discrepancy and driver failure,
successful solver/finalize markers, all detector/ledger completion records,
worldtube and every rank audit at the actual28199 steps, plus original frozen
hashes. It does not alter the original failed batch or frozen completion card,
and is not a general waiver for failed jobs. `summarize_closed_six.py REVIEW`
checks the resulting tables, time denominator, aperture, charge coverage and
output identities and documents this distinction. Run it after visually
reviewing both PNG figures; its report records that review.

`eleven_saturation.py ELEVEN_ROOT OUTPUT` reads only the normally completed
single-node19799 run. It retains the fractional high-current-window rule from
the six-segment comparison, giving a0.993252546ps gate for the four-slip bunch.
The eleven exit spectra use the same physical+/-5% band and+/-520um aperture.
Per-period native power is a separate full-grid/all-retained-frequency quantity;
its charge median is interpolated within one-wavelength slices. Native local
maxima are candidates with explicit10%-prominence and0.5m-separation criteria,
not automatic proof of a universal saturation length.

`present_eleven_saturation.py REVIEW` checks the eleven full-charge exits,
all991 regular period states, CSV, denominator/aperture/BF/Parseval, and original
output identities. It presents the aperture-band exit curve and dense native
curve in separate panels, with their different power definitions visible.
After inspecting the final saturation and bunching PNGs, record that visual
check in the delivery validation; PDF files are same-source exports.
Physical dump steps must match exactly. Iteratively accumulated zbar differs
from analytic exit coordinates by tens of pm late in this run; its position
check allows1nm, far below the350um physical integration displacement. This
does not change the stored coordinates or resample the fields.
