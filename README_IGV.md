# Fusion report with IGV

Keep `Generate_Fusion_reports_html_PDF.py` and `fusion_igv.py` together.
The existing CAB appearance, YAML cover, sample PDF links, and references are retained.

## Input layout

```text
inputs/
  ID1_arriba.tsv
  ID1_fusioncatcher.txt
  ID1_star-fusion.tsv
  ID1.pdf
  ID2_arriba.tsv
  ID2_fusioncatcher.txt
  ID2_star-fusion.tsv
  ID2.pdf
  Lee.config.yaml
  Fusion_reference.bib
  stdRNAseqReportCABlogo.png
bam/
  ID1.bam
  ID1.bam.bai
  ID2.bam
  ID2.bam.bai
```

Use coordinate-sorted BAMs with their matching BAI or CSI indexes, aligned to
the same reference build as the caller coordinates. The script reads BAM headers
to validate sorting, contig names, and coordinate bounds; it does not verify
that an index was built from that BAM. BAMs are not copied into the report.

Also supported: one BAM under `bam/ID1/`, or names such as `ID1.sorted.bam`
and `ID1_Aligned.sortedByCoord.out.bam`. Ambiguous matches are reported rather
than choosing a BAM arbitrarily. Override discovery with, for example,
`--bam-pattern '{sample}/*.sorted.bam'`.

## Select the number of events

```bash
python Generate_Fusion_reports_html_PDF.py \
  --inputs inputs \
  --outdir Fusion_Results_IGV \
  --top-n 10 \
  --igv-top-n 5 \
  --bam-dir bam \
  --igv-genome hg38 \
  --config inputs/Lee.config.yaml \
  --logo inputs/stdRNAseqReportCABlogo.png
```

- `--top-n 5` or `10`: report table candidates per caller.
- `--igv-top-n 5` or `10`: IGV candidates per caller, independently of tables.
  Defaults to 5, so up to 15 events across three callers per sample.
- `--igv-genome`: IGV reference genome ID. Defaults to `SPECIES` in YAML,
  otherwise `hg38`. Set it explicitly when the data use another assembly.
- `--igv-flank 250`: bases on each side of each breakpoint. Both partner regions
  are shown in a multi-locus view.

The selection uses the same per-caller rankings and fusion-pair deduplication
as the report. A pair shared by callers can appear in each caller's list.
For multiple rows of a pair, coordinates come from its highest-ranked row.
FusionCatcher ranks pairs with its existing aggregated maxima; its representative
coordinate row is selected by spanning pairs, then unique reads. The export
is a gene-pair shortlist, not a complete enumeration of breakpoint isoforms.
Malformed or missing coordinates are reported, not invented or backfilled.

## Generate actual IGV screenshots

IGV Desktop must be installed in the environment used to run the report.
On Linux with a graphical session:

```bash
python Generate_Fusion_reports_html_PDF.py \
  --inputs inputs --outdir Fusion_Results_IGV \
  --top-n 10 --igv-top-n 5 --bam-dir bam --igv-genome hg38 \
  --config inputs/Lee.config.yaml --logo inputs/stdRNAseqReportCABlogo.png \
  --igv-command '/path/to/IGV/igv.sh'
```

On a Linux compute node without a display, use an installed virtual display:

```text
--igv-command 'xvfb-run -a /path/to/IGV/igv.sh'
```

Replace the launcher path with your actual installation. The script appends
`-b <generated_batch_file>` and starts a separate IGV process. It runs each
sample sequentially, waits up to 15 minutes per sample, and records the log in
that sample's IGV export folder. The batch closes its IGV process when finished.
Images are embedded in the HTML and PDF reports after IGV finishes.

Alternatively, with the BAM and index already in place:

1. Generate the report to export the batch files.
2. Run each sample's linked `igv_batch.txt` with IGV Desktop (`igv.sh -b PATH`).
3. Run the report again using the same options and output directory.

The script reuses screenshots from the matching export folder. The folder's
identifier includes the events, genome, flank, and BAM/index path/size/timestamp
so changed inputs cannot silently reuse older screenshots. After moving the
folder to another computer, rerun the generator to refresh absolute batch paths.

Without BAM/index files or rendered screenshots, reports explicitly show
screenshots as pending. No artificial alignment image is substituted.

## Interactive HTML

Each sample section includes an **IGV Visualisation** heading with an
**Open <sample> interactive.html** hyperlink and an expandable inline viewer. Select the sample's local BAM and matching BAI/CSI
index, then click **Load sample alignments**. Choose an event from the list
to view both partner loci. Files are read locally by the browser; they are not
uploaded by this viewer. Internet is required for the pinned IGV JavaScript
library and reference genome resources. Browser file-security settings may
restrict local reports; use a browser that permits local file selection.

PDFs contain static screenshots when available, the selected-event table, and
the same interactive.html hyperlink under **IGV Visualisation**. This link opens
the separate HTML viewer; interaction takes place in the browser. Keep the
`igv/` folder alongside the reports when moving or sharing them. The PDF link
uses a relative local-file action; a desktop PDF reader may be needed if a
browser PDF viewer blocks opening it.

At present, actual BAM loading and IGV Desktop
rendering remain unverified for ID1/ID2 because their BAM/index files were not
available in this session.

## Exported IGV files

Under `outdir/igv/<sample>/<selection>/`:

- `interactive.html`: local-file IGV viewer.
- `events.tsv`: caller, rank, fusion pair, breakpoints, and IGV loci.
- `events.json`: selection and input identity metadata.
- `breakpoints.bed`: two one-base BED records per event, one per partner.
- `igv_batch.txt`: IGV Desktop commands for all selected events.
- `snapshots/`: actual IGV PNG output.
- `igv.log`: launcher output when `--igv-command` is used.

Input breakpoints and displayed IGV loci are 1-based. BED records are converted
to zero-based, half-open intervals `[position - 1, position)`.
Keep the entire report directory together when sharing HTML and its exports.

Implementation follows IGV's official command executor and IGV JavaScript API:
- https://github.com/igvteam/igv/blob/main/src/main/java/org/igv/batch/CommandExecutor.java
- https://github.com/igvteam/igv.js

## What is igv.sh?

`igv.sh` is the Unix launcher distributed with IGV Desktop. It prepares the
Java runtime/classpath and starts IGV, forwarding command-line options such
as `-b` to IGV. Use the launcher from the same installed IGV distribution;
its Java options and directory layout are version-dependent.

The report generator does not create IGV itself. It creates `igv_batch.txt`,
which contains commands such as:

```text
new
genome hg38
load /path/to/bam/ID1.bam
snapshotDirectory /path/to/report/igv/ID1/snapshots
goto 18:63344160-63344660 18:63125403-63125903
snapshot arriba_01.png
exit
```

The generated batch also loads the breakpoint track and configures paired-read
viewing. The generator invokes the installed launcher as:

```bash
/path/to/IGV/igv.sh -b /path/to/generated/igv_batch.txt
```

An optional wrapper, `run_igv.sh`, is included in this package. It is not a
replacement IGV installation. Configure it by setting the installed directory:

```bash
export IGV_HOME=/path/to/IGV
bash run_igv.sh -b /path/to/generated/igv_batch.txt
```

To use that wrapper in the report command:

```text
--igv-command 'bash /path/to/run_igv.sh'
```

Its core action is `exec bash "$IGV_HOME/igv.sh" "$@"`: it passes the arguments
unchanged to the original launcher. On a headless Linux node with Xvfb installed,
use `--igv-command 'xvfb-run -a bash /path/to/run_igv.sh'` instead.
