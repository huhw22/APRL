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
