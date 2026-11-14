# Retired-results review

Read-only inventory and user-confirmable retirement plan for historical APRL
results. Inventory, planning, archiving and validation never delete, rename or
move original files. They do not submit simulations or alter frozen campaigns.

Run `inventory.py WORKSPACE OUTPUT`, `audit_loading.py WORKSPACE OUTPUT`,
`build_plan.py WORKSPACE OUTPUT`, `archive_and_report.py WORKSPACE OUTPUT`, then
`validate.py WORKSPACE OUTPUT`. Inventory reads lstat metadata only; the loading
audit explicitly reads three closed historical input files, never running
scientific output. The evidence archive copies only selected small documents,
plots and logs and verifies their bytes.

Deletion candidates are specific H5/bin/raw files. Original inputs, actual
initial import dumps, code, plots, tables, configuration, logs and all ledger/
lab-audit directories remain. The old complete six-segment particle/energy
planes stay; its worldtube raw integral arrays are a separately selectable B05.
Current campaigns, their frozen dependencies and declared input aliases are
excluded, including protected symlink targets and hardlinked inodes.

Reclaim estimates count allocated512-byte blocks only when every hardlink of
the inode is selected. Unknown external links are conservatively worth zero.
The CSV/JSONL records device,inode,nlink,size,allocation andmtime. Recheck all
identities and current job dependencies after the user selects groups and
before any later deletion. Suggested groups are not approval.

This plan preserves useful diagnostics even where the old physics claim was
withdrawn. A failed Slurm status, large macro weight, cancelled task or newer
replacement alone does not prove that old physical output is wrong.

After the user explicitly approved A01–A08,B01,B03,`execute_approved.py` was
added for that exact dated manifest only. Its separate preflight verifies the
reviewed manifest hashes,current frozen dependencies,archive hashes and every
selected file identity. The execute mode unlinks only the approved per-file
allowlist using directory descriptors and O_NOFOLLOW,with before/after journal
records and a final preserved-file/dependency/archive check. It never removes
whole directories. Do not reuse this fixed authorization for another plan or
additional groups. The original review validation remains a historical record.

On 2026-11-15 the user separately authorized a broader cleanup to bring
APRL_Workspace below 500 decimal GB. `compact_workspace.py` implements only
that dated campaign in `results/workspace_compaction_20261115_r1`. This retires
obsolete valid raw outputs as well as failed tests; retirement does not retract
their derived scientific conclusions. The previous group's authorization is
not reused. The new JSON authorization, metadata inventory, exact file list,
plan hashes and preflight are required before `execute` can run once.

This policy preserves the entire active w355 campaign, immutable dependencies,
declared and recognizable original inputs (including native electron dump 0),
code and programs, plots, CSV/JSON, cards, logs, compact regression data, the
selected canonical physical controls and the old six-segment worldtube/ledger.
For the older independent 60M benchmark it keeps original inputs, the final
substep particle plane, final optical plane and ledger; intermediate planes
are retired. Large old H5/bin/raw/NPY/NPZ outputs are eligible elsewhere.
All aliases of preserved inodes stay; reclaim requires observing and selecting
every hardlink. No directories or symbolic links are removed.

Run `compact_workspace.py plan`, inspect its exact manifest and protected
scope, then `preflight`, then `execute`. Execution rechecks each inode through
the same guarded unlink helper, journals both sides of every removal, verifies
preserved large files and rehashes the active frozen dependencies. Verify the
actual workspace size with `du -x -B1` afterwards: the inventory's projected
regular-file total excludes directory overhead and ongoing output growth.
The campaign is not an automatic retention service; a future cleanup requires
fresh task context and authorization.
