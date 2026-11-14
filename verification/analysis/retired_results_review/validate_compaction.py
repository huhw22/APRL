"""Independently validate the completed 2026-11-15 compaction; never delete."""
import collections
import hashlib
import json
import os
from pathlib import Path
import sqlite3

from execute_approved import current_frozen_checks, match, verify_archive

W = Path('/public/home/wfollowm/APRL_Workspace')
R = W / 'results/workspace_compaction_20261115_r1'
CURRENT = W / 'results/fel_independent_w355_20261114_r1'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plan = read(R / 'PLAN.json')
    complete = read(R / 'EXECUTION_COMPLETE.json')
    assert complete['complete']
    rows = [json.loads(line) for line in (R / 'DELETE_FILES.jsonl').read_text().splitlines()]
    chosen = {row['path']: row for row in rows}
    assert len(rows) == len(chosen) == plan['delete_files'] == complete['deleted_files']
    journal = [json.loads(line) for line in (R / 'DELETE_JOURNAL.jsonl').read_text().splitlines()]
    assert len(journal) == 2 * len(rows)
    reclaimed = 0
    for row, intent, deleted in zip(rows, journal[::2], journal[1::2]):
        assert intent == dict(event='intent', **row)
        assert deleted['event'] == 'deleted'
        assert {key: deleted[key] for key in row} == row
        reclaimed += deleted['reclaimed_bytes']
        assert not os.path.lexists(row['path'])
    by_inode = collections.defaultdict(list)
    for row in rows:
        by_inode[row['dev'], row['inode']].append(row)
    assert all(len(group) == group[0]['nlink'] for group in by_inode.values())
    assert reclaimed == sum(group[0]['allocated_bytes'] for group in by_inode.values())
    assert reclaimed == plan['planned_reclaim_bytes'] == complete['reclaimed_bytes']
    assert not any(Path(path).is_relative_to(CURRENT) for path in chosen)
    assert not any(Path(path).name.endswith('_electrons_0.h5') for path in chosen)

    # This checks existence for every old, unselected regular file. It does not
    # open scientific H5 content or impose static metadata on growing outputs.
    retained_count = 0
    with sqlite3.connect(R / 'inventory.sqlite') as db:
        for (path,) in db.execute("SELECT path FROM files WHERE kind='file'"):
            if path in chosen or Path(path).is_relative_to(CURRENT):
                continue
            assert os.path.lexists(path), ('Unselected historical file missing', path)
            retained_count += 1
    kept = read(R / 'PRESERVED_STATIC_BEFORE.json')
    for row in kept:
        match(row, Path(row['path']).lstat())
    frozen = current_frozen_checks(W, R)
    assert frozen == read(R / 'CURRENT_FROZEN_BEFORE.json')
    assert frozen == read(R / 'CURRENT_FROZEN_AFTER.json')
    archive_count = verify_archive(W / 'results/retired_results_review_20261114_r1')
    package = W / 'results/accepted_physics_transfer_20261112'
    package_count = 0
    for line in (package / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert sha(package / name) == digest
        package_count += 1
    preflight = read(R / 'PREFLIGHT.json')
    assert sha(R / 'PLAN.json') == preflight['plan_sha256']
    assert sha(R / 'DELETE_FILES.jsonl') == preflight['filelist_sha256']
    result = dict(
        passed=True, deleted_files=len(rows), journal_records=len(journal),
        exact_plan_and_journal_match=True, all_selected_paths_absent=True,
        independent_reclaimed_bytes=reclaimed, unselected_historical_files_exist=retained_count,
        preserved_static_large_files_unchanged=len(kept),
        current_frozen_dependencies_unchanged=True,
        current_frozen_unique_inodes_hashed=frozen['unique_inodes_hashed'],
        previous_error_archive_sha256_passed=archive_count,
        accepted_physics_package_sha256_passed=package_count,
        native_initial_electron_dumps_preserved=True, no_running_scientific_output_opened=True,
        actual_du_check_is_separate=True,
    )
    (R / 'VALIDATION_AFTER_DELETE.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
