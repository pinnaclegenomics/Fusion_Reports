# CAB fusion reports

Keep `Generate_Fusion_reports_html_PDF.py` and `fusion_igv.py` together, and use the generator in place of the previous script.
From your project folder, run:

```bash
python Generate_Fusion_reports_html_PDF.py \
  --inputs inputs \
  --outdir Fusion_Results_CAB_top10 \
  --top-n 10 \
  --config inputs/Lee.config.yaml \
  --logo inputs/stdRNAseqReportCABlogo.png
```

Use `--top-n 5` for five fusion pairs per caller (the default). Both settings
preserve the original ranking and deduplication. If a caller has fewer
candidates, the report includes those available. Use separate output folders
to retain both versions; HTML/PDF filenames are the same for both modes.

## Appearance and cover information

- CAB orange `#FF7F00`, taken from the supplied R Markdown template.
- Orange headings, table headers, selected sample tabs, and chart bars;
  gray distinguishes a second read-support series.
- Embedded CAB logo in individual and combined HTML reports and individual PDFs.
- A PDF cover page with PI, project, task, project lead, department, SRM,
  reference genome, analysis team, contact, start/completion dates, and
  the actual report generation time.
- Only populated YAML fields are displayed. A `PIPELINE` value is optional;
  no pipeline version is invented. Internal input/output paths are not displayed.
- Both `PROJECT_LEAD_NAME` and the existing `PPOJECT_LEAD_NAME` key are accepted.
  If both are supplied, the correctly spelled key takes precedence.
- `--config` and `--logo` are optional to preserve existing command usage.

Project information comes from the supplied YAML, including its dates and
project lead, rather than from the reference PDF. The three callers remain
Arriba, FusionCatcher, and STAR-Fusion.

## Dependencies

```bash
pip install pandas matplotlib weasyprint PyYAML pypdf
```

WeasyPrint also requires platform-specific native libraries (including Pango).
Use the same working WeasyPrint environment as the original report generator.
If PDF conversion fails, the HTML and TSV reports are still produced, and
`batch_run_summary.tsv` records the error.

## Output files

- `fusion_report_all_samples.html`: combined report with sample tabs.
- `<sample>_fusion_report.html` and `<sample>_fusion_report.pdf`.
- `<sample>_top5_summary.tsv` or `<sample>_top10_summary.tsv`.
- `all_samples_top5_summary.tsv` or `all_samples_top10_summary.tsv`.
- `batch_run_summary.tsv`.

The `CAB_preview_top10` folder contains example reports generated from the
provided ID1 and ID2 inputs using the supplied YAML and logo.

## Supporting PDFs and references

The script automatically finds `<sample>.pdf` (for example, `ID1.pdf` and
`ID2.pdf`) in `--inputs`. Each sample report links to its matching PDF.
If a sample PDF is absent, no broken link is added.

`Fusion_reference.bib` in the input folder is detected automatically. You can
select a different file with `--references inputs/Fusion_reference.bib`.
The supplied mixed resource-prose/BibTeX file is rendered as Additional
Resources, Acknowledgements, and References at the end of each sample report,
and once at the bottom of the combined HTML. Publication DOI links are clickable.
The HTML also includes a download link for the unchanged source reference file.

Supporting PDFs and the source reference file are copied into `attachments/`
inside the output folder. Keep that folder alongside the reports when moving
or sharing them. HTML links are relative; generated PDF links use relative
remote-PDF actions. Opening those actions depends on the PDF viewer's support
for local linked files (use a desktop PDF reader if a browser viewer blocks them).

The displayed resource link beginning `ttps://` is normalized to `https://`.
Citation metadata is otherwise taken from your file without bibliographic correction.
The script supports braced, quoted and numeric BibTeX fields; advanced BibTeX
macros/string concatenation are not supported.

## IGV review

See `README_IGV.md` for `--igv-top-n`, BAM discovery, screenshot generation,
and the interactive HTML viewer. The previews use top-10 report tables and
top-5 IGV selections per caller. Actual BAMs were unavailable, so screenshots
are explicitly pending in these previews.

## Shared three-tool background

The information adapted from `Gene_Fusion_Aug2026.docx` is now a fixed section
before the sample results. It includes workflow stages, split-read and
paired-read evidence, a comparison of STAR-Fusion/Arriba/FusionCatcher, and
individual tool profiles. It is identical across sample reports and appears
once in the combined report.

Keep `fusion_tool_background.json` alongside the two Python files. It is the
single shared text source for future reports; no Word-reading dependency is
needed at report-generation time. To use a revised shared template, edit that
file or pass `--tool-background /path/to/fusion_tool_background.json`.

The source document's broad performance claims are presented as qualitative
background rather than measured results or clinical guarantees. Workflow
stages are distinguished from processing actually performed by this generator.
The FastQC description distinguishes quality assessment from read trimming,
and assembly is described as caller/workflow-dependent. The stray “Export to
Sheets” line in the source document is omitted.

## HTML and PDF Contents navigation

Individual reports and the combined report now use a left-side Contents panel,
matching the floating overview navigation in the supplied CAB RNA-seq HTML
reference. It includes Overview, workflow/evidence, caller comparison/profiles,
sample results, IGV review, and references. The active section is highlighted
in CAB orange. Selecting a sample's sidebar link activates that sample's tab.
Links to sections can be bookmarked using their URL fragment.

On narrow screens, Contents collapses into an expandable menu. The sidebar
uses embedded CSS/JavaScript with no added library dependency.

Each PDF page also has a left-side Contents panel with clickable section links.
PDF pages are 255 × 297 mm to preserve the original body text and table sizes
alongside the panel. Run commands and input requirements are unchanged.
