"""Freeze only an explicitly declared duration prefix; never launch a job."""
import argparse
from freeze_campaign import main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root')
    parser.add_argument('--binary', required=True)
    parser.add_argument('--build-record', required=True)
    parser.add_argument('--profile', type=int, default=0)
    args = parser.parse_args()
    main(args.root, args.binary, args.build_record, args.profile, bounded_prefix=True)
