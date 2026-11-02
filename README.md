<link rel="icon" type="image/png" href="[Nextflow_Fusion_Calling_Workflow_Diagram.png](https://github.com/pinnaclegenomics/Fusion_Reports/blob/8e66e8cb41bb9a7c7da82fce76b92bd0f6dc17c4/Nextflow%20Fusion%20Calling%20Workflow%20Diagram.png)">

# CAB fusion reports — Nextflow

This DSL2 pipeline generates your CAB HTML and PDF fusion reports from existing
Arriba, FusionCatcher and STAR-Fusion results. It retains the YAML information,
logo, fixed three-tool background, references, supporting sample PDF links,
IGV exports, interactive.html hyperlinks, and HTML/PDF Contents panels.
It does not run alignment or the fusion callers themselves.

The existing Python generator is bundled under `scripts/`. Nextflow runs it as
one batch reporting task, preserving the combined multi-sample report. Samples
are processed sequentially inside that task; this version does not parallelize
samples. Nextflow manages execution, staging, caching, resources and publishing.
A validation step fails the task if any sample HTML/PDF or IGV viewer is missing,
or the generator reports a sample/PDF failure.

## 1. Requirements

Install Nextflow (24.10 or newer) with a compatible Java runtime according to
[the official installation guide](https://www.nextflow.io/docs/latest/install.html).
This package was execution-tested with Nextflow 26.04.6 and Java 24.

Choose one report environment:

- **Conda:** install Conda, then use `-profile conda`. Nextflow creates the
  environment from `environment.yml`, including Python and native PDF libraries.
- **Existing environment:** activate the Python environment that already runs
  the supplied reports, and omit `-profile conda`. It must provide pandas,
  matplotlib, weasyprint, PyYAML, pypdf, Pango and fonts.

IGV Desktop, its compatible Java runtime, and (on headless Linux) Xvfb are
separate optional installations for screenshots. The Conda file covers report
rendering only. Interactive HTML uses the existing IGV.js viewer and requires
internet for its library and genome resources.

## 2. Arrange the input files

```text
Fusion_Project_Scripts_Testing/
├── inputs/
│   ├── ID1_arriba.tsv
│   ├── ID1_fusioncatcher.txt
│   ├── ID1_star-fusion.tsv
│   ├── ID1.pdf
│   ├── ID2_arriba.tsv
│   ├── ID2_fusioncatcher.txt
│   ├── ID2_star-fusion.tsv
│   ├── ID2.pdf
│   ├── Lee.config.yaml
│   ├── Fusion_reference.bib
│   └── stdRNAseqReportCABlogo.png
└── bam/                       # optional without screenshots
    ├── ID1.bam
    ├── ID1.bam.bai
    ├── ID2.bam
    └── ID2.bam.bai
```

BAMs must be coordinate-sorted and have matching BAI/CSI indexes. Use the same
genome build as the fusion calls. The generator also supports sample
subdirectories and its existing BAM discovery patterns. Input files must be
available on the execution machine or shared cluster filesystem.

## 3. Run HTML and PDF reporting

Extract this package and enter its `CAB_fusion_nextflow` directory. Edit the
paths in `params.example.yaml`, then run:

```bash
nextflow run main.nf -profile conda -params-file params.example.yaml
```

Or specify your project paths directly (the variable below can be changed):

```bash
FUSION_PROJECT=/Volumes/groups/cab/projects/automapper/common/ykumar/Transcriptomic_Group/Projects/Fusion_Project_Scripts_Testing

nextflow run main.nf -profile conda \
  --inputs "$FUSION_PROJECT/inputs" \
  --report_config "$FUSION_PROJECT/inputs/Lee.config.yaml" \
  --logo "$FUSION_PROJECT/inputs/stdRNAseqReportCABlogo.png" \
  --top_n 10 \
  --igv_top_n 5 \
  --igv_genome hg38 \
  --outdir "$FUSION_PROJECT/Fusion_Results_Nextflow"
```

`--top_n` and `--igv_top_n` each accept **5 or 10 per caller per sample**.
The defaults are 5. Select the correct genome for your data.
Nextflow parameter names use underscores. `--report_config` is the project
YAML; Nextflow's `-c` option is reserved for Nextflow configuration files.
Use `nextflow run main.nf --help` for the parameter overview.

## 4. Include actual IGV screenshots

Add the following parameters to the reporting command:

```bash
  --bam_dir "$FUSION_PROJECT/bam" \
  --igv_command 'xvfb-run -a /path/to/IGV/igv.sh'
```

Replace the IGV path with the installed launcher. With a graphical desktop,
use `/path/to/IGV/igv.sh` without Xvfb. For nested BAM files, optionally add:

```bash
  --bam_pattern '{sample}/*.sorted.bam'
```

The launcher must be accessible on the worker node. Each sample's screenshot
run can take up to 15 minutes; increase Nextflow's task time for large batches.
Without BAM/index files or successful rendering, screenshots are marked pending;
interactive viewers and breakpoint exports are still generated. The pipeline
validates reports, not the presence of actual screenshots: inspect the IGV status
and `igv.log` before interpreting screenshot completeness.

## 5. Find and share the results

Reports are published as a single directory to preserve all relative links:

```text
<outdir>/reports/
├── fusion_report_all_samples.html
├── ID1_fusion_report.html
├── ID1_fusion_report.pdf
├── ID2_fusion_report.html
├── ID2_fusion_report.pdf
├── all_samples_top10_summary.tsv
├── batch_run_summary.tsv
├── software_versions.json
├── attachments/
└── igv/
    └── <sample>/<selection>/
        ├── interactive.html
        ├── events.tsv
        ├── events.json
        ├── breakpoints.bed
        ├── igv_batch.txt
        └── snapshots/          # when screenshots were generated
```

Share the entire `reports/` directory. BAMs are not bundled; viewers select
local BAM/index files in the browser. PDF interactive.html links may require a
desktop PDF reader. Desktop IGV batch files contain absolute execution paths;
after deleting Nextflow work directories or moving to another machine, regenerate
the batch files before running them manually.

## 6. Resume and cluster execution

Repeat the same command with `-resume` to reuse successful cached work:

```bash
nextflow run main.nf -profile conda -params-file params.example.yaml -resume
```

Keep Nextflow's `.nextflow/` and `work/` directories for resume. Changes to the
staged scripts, options or input contents invalidate the batch task. Deep cache
hashing includes BAM contents and can take time with large inputs. Use a new
`--outdir` for each distinct cohort/top-N setting to avoid confusing older files
with a new run. Inputs and BAM directories should not contain the output or
Nextflow work directories.

For Slurm with Conda:

```bash
nextflow run main.nf -profile conda,slurm \
  -params-file params.example.yaml -c conf/slurm.example.config
```

Edit the partition/account in `conf/slurm.example.config` first. All input,
software and work paths must be visible to compute nodes. Defaults are 2 CPUs,
8 GB memory and 8 hours for the entire cohort; adjust them in your site config.
For an existing cluster environment, use `-profile slurm` and load the required
modules through your site config. Other schedulers can be selected with a
Nextflow executor configuration.

## Parameters

| Parameter | Default | Purpose |
|---|---|---|
| `inputs` | required | Folder of existing caller outputs and optional sample PDFs |
| `outdir` | `results_fusion` | Parent directory for published `reports/` |
| `top_n` | 5 | Fusion pairs shown per caller: 5 or 10 |
| `igv_top_n` | 5 | IGV events per caller: 5 or 10 |
| `report_config` | none | Project information YAML |
| `logo` | none | CAB logo file |
| `references` | auto-detect | Reference file; defaults to `inputs/Fusion_reference.bib` |
| `tool_background` | bundled JSON | Override the shared three-tool background |
| `bam_dir` | none | BAM/index directory |
| `bam_pattern` | auto-discovery | Relative BAM glob containing `{sample}` |
| `igv_genome` | YAML SPECIES or hg38 | IGV reference genome |
| `igv_flank` | 250 | Bases flanking each breakpoint |
| `igv_command` | none | Installed IGV launcher command |
| `python` | python3 | Python executable path, not a shell command |

The implementation uses standard Nextflow [process inputs and outputs](https://www.nextflow.io/docs/latest/process.html).
See `README_IGV.md` and `README_Fusion_reports.md` for detailed report behavior.
