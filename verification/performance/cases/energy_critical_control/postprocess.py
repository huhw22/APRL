"""Run the paired review only after both bounded PIC cases have closed."""
import hashlib
import json
from pathlib import Path
import sys
import traceback


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main(output):
    output = Path(output).resolve()
    plan = json.loads((output/'plan.json').read_text())
    frozen = json.loads((output/'FROZEN.json').read_text())

    def check_identity():
        for path, expected in frozen.items():
            if digest(path) != expected:
                raise RuntimeError('Postprocessing identity changed: ' + path)

    check_identity()
    root = Path(plan['campaign'])
    missing = []
    for mode in ['off', 'on']:
        complete = root/mode/'PIC_COMPLETE.json'
        result = root/mode/'RESULT.json'
        if not complete.exists() or not result.exists():
            missing.append(mode)
            continue
        completion = json.loads(complete.read_text())
        data = json.loads(result.read_text())
        if not (completion.get('analysis_complete') is True
                and completion.get('stage') == 'full'
                and data.get('program_integrity_pass') is True
                and data.get('job') == str(plan['jobs'][mode])):
            missing.append(mode)
    if missing:
        save(output/'RESULT.json', dict(postprocessing_complete=False,
            status='skipped_missing_or_unaccepted_case', cases=missing,
            scientific_acceptance=False, raw_h5_opened=False))
        print('Paired review skipped: case closure is missing or unaccepted.', flush=True)
        return
    for mode in ['off', 'on']:
        identity = json.loads((root/mode/'FROZEN.json').read_text())
        for path, expected in identity.items():
            if digest(path) != expected:
                raise RuntimeError('Completed case identity changed: ' + path)
    sys.path.insert(0, str(root/'scripts'))
    from analyze import compare as compare_bunching
    sys.path.insert(0, str(output/'scripts'))
    from mode_review import compare as compare_modes
    compare_bunching(root)
    mode_output = compare_modes(root)
    check_identity()
    save(output/'RESULT.json', dict(postprocessing_complete=True,
        scientific_acceptance=False, manual_figure_review_pending=True,
        bunching=str(root/'COMPARISON.json'), modes=str(mode_output/'RESULT.json'),
        figures=[str(root/'bunching_comparison.png'),
                 str(mode_output/'mode_comparison.png')],
        interpretation='Check reproduced fast startup in the off control, measured mode suppression, and substep BF together before interpreting causality.'))
    print('Paired BF and measured-mode review complete.', flush=True)


if __name__ == '__main__':
    try:
        main(sys.argv[1])
    except Exception:
        error = traceback.format_exc()
        Path(sys.argv[1], 'FAILURE.txt').write_text(error)
        raise
