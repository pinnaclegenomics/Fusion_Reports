#!/usr/bin/env python3
"""Fail the Nextflow task if the batch generator caught any sample/PDF error."""
import argparse
import csv
import json
from pathlib import Path
from Generate_Fusion_reports_html_PDF import discover_samples


def validate(inputs, reports):
    expected = set(discover_samples(inputs))
    if not expected:
        raise RuntimeError('No recognizable caller results found.')
    with (reports / 'batch_run_summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len(rows) != len(expected) or {r['Sample'] for r in rows} != expected:
        raise RuntimeError('Batch summary does not match the discovered samples.')
    for row in rows:
        if row['Status'] != 'OK' or row['PDF_Status'] != 'OK':
            raise RuntimeError(f"Report failed for {row['Sample']}: {row['Error']}")
        for suffix in ('html', 'pdf'):
            output = reports / f"{row['Sample']}_fusion_report.{suffix}"
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(f'Missing or empty report: {output}')
    if not (reports / 'fusion_report_all_samples.html').is_file():
        raise RuntimeError('Combined HTML report is missing.')
    viewers = sorted(reports.glob('igv/*/*/interactive.html'))
    if len(viewers) != len(expected):
        raise RuntimeError('An interactive IGV viewer is missing.')
    import importlib.metadata
    versions = {name: importlib.metadata.version(name)
                for name in ['pandas', 'matplotlib', 'weasyprint', 'PyYAML', 'pypdf']}
    (reports / 'software_versions.json').write_text(json.dumps(versions, indent=2) + '\n')
    print(f'Validated HTML, PDF, and interactive IGV exports for {len(expected)} samples.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--reports', type=Path, required=True)
    args = parser.parse_args()
    validate(args.inputs, args.reports)
