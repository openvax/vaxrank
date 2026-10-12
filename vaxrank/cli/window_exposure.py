"""Assess every SLP window from saved scored antigen evidence."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from mhctools import format_vaccine_window_report, vaccine_window_rows

from ..core_exposure import assess_dataset_windows
from ..epitope_dataset import EpitopeDataset


def main(argv=None):
    parser = argparse.ArgumentParser(prog='vaxrank-window-exposure', description=__doc__)
    parser.add_argument('--input-epitopes', required=True, help='Enriched native --output-epitopes dataset with explicit antigen contexts')
    parser.add_argument('--profile', required=True, help='Scenario JSON: times, sourced rates and optional candidate-specific profiles')
    parser.add_argument('--min-length', type=int, required=True)
    parser.add_argument('--max-length', type=int, required=True)
    parser.add_argument('--n-term', required=True, choices=('free', 'acetylated', 'unknown'))
    parser.add_argument('--c-term', required=True, choices=('free', 'amidated', 'unknown'))
    parser.add_argument('--output-dir', required=True, help='New directory for request, audit JSON, per-core CSV and report')
    args = parser.parse_args(argv)
    try:
        output = Path(args.output_dir)
        if output.exists():
            raise FileExistsError('Output already exists: %s' % output)
        profile_data = Path(args.profile).read_bytes()
        dataset = EpitopeDataset.load(args.input_epitopes)
        request, result = assess_dataset_windows(
            dataset, json.loads(profile_data), min_length=args.min_length, max_length=args.max_length,
            n_term=args.n_term, c_term=args.c_term)
        result['input_files'] = [{'path': path, 'sha256': hashlib.sha256(data).hexdigest()}
                                 for path, data in ((args.profile, profile_data),
                                     (args.input_epitopes, Path(args.input_epitopes).read_bytes()))]
        rows = vaccine_window_rows(result)
        report = format_vaccine_window_report(result)
        report += '\n## Source coverage\n\nAll admitted target-containing windows were evaluated before ranking or count limits. '
        report += 'Targets use supplied core intervals when documented, otherwise the saved scored ligands; no shorter binding core was inferred.\n\n'
        for gap in result['source_coverage_gaps']:
            report += '- %s: %s (%s).\n' % (gap['prediction_id'].replace('\n', ' '), gap['sequence'], gap['reason'])
        if not result['source_coverage_gaps']:
            report += 'No missing source context or saved-score gaps.\n'
        output.mkdir(parents=True, exist_ok=False)
        for filename, value in (('request.json', request), ('exposure.json', result)):
            with (output / filename).open('x') as stream:
                json.dump(value, stream, indent=2, allow_nan=False)
                stream.write('\n')
        with (output / 'cores.csv').open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else ['candidate_id', 'status'])
            writer.writeheader()
            writer.writerows(rows)
        (output / 'report.md').write_text(report)
    except (KeyError, TypeError, ValueError, OSError, RuntimeError) as error:
        parser.error(str(error))
    print(output / 'report.md')


if __name__ == '__main__':
    main()
