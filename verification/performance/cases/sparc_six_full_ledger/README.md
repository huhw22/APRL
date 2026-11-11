# Six-segment SPARC complete energy-ledger pair

The user requested the current 1,734,656-particle input without regenerating,
reweighting or lengthening the bunch, six undulators, and two methods:
conservative-source + linear-binomial3 and conservative-source +
quintic-compensated. Grid spacing, 40-by-40 transverse cells, CPML6, boost,
subcycling and the first six original magnetic elements are preserved.

The original eleven-segment momentum job 19708 was cancelled by explicit user
request. Its output is retained and must be treated as incomplete.

## Independent physical and discrete records

- Full basic boost ledger, every Maxwell step, plus the strict discrete audit.
- Lab field energy, longitudinal momentum, Maxwell stress, and the complete
  Cowan/CPML/inlet/reference-current/remap audit.
- Equal-time E/B lab energy-current integrals on every longitudinal cell and
  every step; four transverse faces and the independent null-time caps can be
  integrated together with entry and exit ports.
- Existing substep particle-plane diagnostics provide scalar endpoint lab
  energies (rather than recomputing gamma from interpolated three-momentum).
- The additional read-only side-particle diagnostic records outward and inward
  transverse crossings from committed substeps, after the same terminal cutoff
  used by the normal diagnostics. Trial orbits never write records.

The physical observation tube starts at z=0.03 m, downstream of the bounded
frozen-source cloud, and reaches the sixth exit at z=14.946 m. Its transverse
half-width is 520 um, one cell inward from the prescribed surface-current sheets
at 560 um. Source exclusion, particle time-cap occupancy, side transport and all
surface coverage must be checked after completion. No surface term may be
inferred from the energy residual. Magnetostatic lab devices have zero physical
electric work in this tube. Cowan/time-discretization channels remain separate
from the physical lab balance.

## Preparation and launch

Run prepare.py with a new result directory. It verifies the original input SHA,
copies the frozen worldtube observer source, and overlays only diagnostic
initialization, the particle-local audit trace and serial clipped event commits.
The validated original physics archive is reused. The source changes do not
write to the particle, field, CPML, current or solver state.

preflight.sbatch builds a separate executable, runs side-event manufactured
checks, verifies 16-by-8 binding, and compares 40-step same-input diagnostic off/on
outputs bitwise. The original physical stop remains downstream of the complete
six-segment line. Therefore this intentionally short preflight uses the existing
duration-exhaustion mechanism: its physical complete flags remain false, and the
driver requires the closed-output readiness marker and exactly 40 steps. It is a
software/diagnostic check, not a complete six-segment scientific result.

main.sbatch takes the result directory and either energy-binomial3 or
energy-momentum. Each production allocation must depend afterok on the preflight
job, and the runner independently verifies the preflight gate, binary and all
frozen inputs. Complete production output requires the configured physical stop,
26923 steps, all H5 complete markers and every rank's audit completion marker.

The first successful setup is results/sparc_six_full_ledger_20261112_r3.
Earlier r1 stopped before PIC because binding-probe layout arguments were missing;
r2 rejected a short reference stop upstream of the full line. Both are retained
as setup failures. Their dependent production jobs were cancelled without running.
The r3 scripts preserve the full physical stop and use an explicitly incomplete
40-step preflight. No failure is counted as scientific evidence.

On 2026-11-12, preflight 19742 completed 0:0: side-event unit checks, 16-by-8
binding, both observer-enabled methods and the three-point diagnostic off/on
bitwise comparison passed. Jobs 19743 (node6) and 19744 (node2) each passed
binding and advanced at least 100 steps with the complete observers enabled.
STARTUP_CONFIRMED.json in the result directory records this limited startup
evidence. Complete physical energy closure remains unassessed until the data
are closed and all physical surface terms and coverage are checked.
