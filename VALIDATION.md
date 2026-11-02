# Validation — 18 September 2026

Executed locally with Nextflow 26.04.6, Java 24.0.2, the existing Python 3.12
report environment, and the supplied ID1/ID2 caller files.

Passed:

- Full Nextflow reporting task: both sample HTML/PDFs, combined HTML, supporting
  attachments and interactive IGV exports generated and published together.
- Top-10 tables with top-5 IGV events: each sample had 10 rows per caller and
  5 IGV events per caller.
- Top-5 tables with top-10 IGV events: each sample had 5 rows per caller and
  10 IGV events per caller.
- Optional inputs: provided YAML/logo; omitted YAML/logo; default references;
  explicit references/background; omitted BAM directory and empty BAM directory.
- Each full CAB PDF had 10 pages with 19 valid sidebar links per page; the
  interactive HTML link resolved relative to the published report directory.
- Resume reused the successful reporting task.
- A PDF-library failure correctly failed the Nextflow task through the output
  validator, even though the underlying generator continued producing HTML.

Not execution-tested: Conda environment creation, Slurm submission, and actual
IGV BAM loading/screenshot rendering. Sample BAM/index files were unavailable.
Screenshot panels remain explicitly pending in these runs.
