#!/usr/bin/env python3
"""
Batch RNA-Seq Fusion HTML Report Generator
==========================================

Input folder example:
    inputs/
      ID1_arriba.tsv
      ID1_fusioncatcher.txt
      ID1_star-fusion.tsv
      ID2_arriba.tsv
      ID2_fusioncatcher.txt
      ID2_star-fusion.tsv

Output:
    outdir/
      fusion_report_all_samples.html      <-- combined HTML with one tab/sheet per sample
      ID1_fusion_report.html              <-- individual sample HTML
      ID1_fusion_report.pdf               <-- individual sample PDF
      ID1_top5_summary.tsv
      ID2_fusion_report.html
      ID2_fusion_report.pdf
      ID2_top5_summary.tsv
      ...
      all_samples_top5_summary.tsv
      batch_run_summary.tsv

Select the maximum number of fusion pairs per tool (default: 5):
    python Generate_Fusion_reports_html_PDF.py -i inputs -o results_top5 --top-n 5
    python Generate_Fusion_reports_html_PDF.py -i inputs -o results_top10 --top-n 10

With --top-n 10, summary filenames use top10 instead of top5.
If fewer unique fusion pairs are available, all available pairs are included.
Ranking and deduplication are unchanged.

CAB branding and project information:
    python Generate_Fusion_reports_html_PDF.py -i inputs -o results --top-n 10 \
        --config inputs/Lee.config.yaml --logo inputs/stdRNAseqReportCABlogo.png

The optional YAML supplies cover-page information; missing fields are omitted.
Both PROJECT_LEAD_NAME and the legacy PPOJECT_LEAD_NAME key are supported.
The optional PNG/JPEG logo is embedded in HTML for portable reports.
Matching <sample>.pdf files and Fusion_reference.bib are detected in inputs.
Use --references PATH to choose another reference file. Supporting files are
copied to outdir/attachments; keep that folder with the reports when sharing.

The shared fusion_tool_background.json template is included before results in
every report. Keep it beside this script together with fusion_igv.py.

Install:
    pip install pandas matplotlib weasyprint PyYAML pypdf

PDF generation uses WeasyPrint. If WeasyPrint is not installed, HTML reports
are still generated and the batch summary records the PDF-generation error.
"""

import argparse
import base64
import html
import io
import json
from html.parser import HTMLParser
import re
import shutil
import unicodedata
from urllib.parse import quote, unquote, urlparse
from datetime import datetime
from pathlib import Path
from fusion_igv import build_igv, report_html as igv_report_html

import pandas as pd
import matplotlib.pyplot as plt

try:
    from weasyprint import HTML as WeasyHTML
    WEASYPRINT_AVAILABLE = True
    WEASYPRINT_IMPORT_ERROR = ""
except Exception as exc:
    WeasyHTML = None
    WEASYPRINT_AVAILABLE = False
    WEASYPRINT_IMPORT_ERROR = str(exc)


# CAB palette from the reference R Markdown template.
CAB_ORANGE = "#FF7F00"
CAB_GRAY = "#555555"
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=[CAB_ORANGE, CAB_GRAY])
plt.rcParams["axes.titlecolor"] = "#B85B00"


def load_tool_background(path=None):
    """Render one shared, data-only background template for every sample report."""
    path = Path(path) if path is not None else Path(__file__).with_name("fusion_tool_background.json")
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)
    esc = lambda text: html.escape(str(text))
    def table(block):
        header = "".join("<th>" + esc(cell) + "</th>" for cell in block["columns"])
        rows = "".join("<tr>" + "".join("<td>" + esc(cell) + "</td>" for cell in row) + "</tr>" for row in block["rows"])
        return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{rows}</tbody></table></div>'
    def title(text):
        return '<h3 class="section-title">' + esc(text) + '</h3>'
    workflow, metrics, comparison, profiles = (data[key] for key in ("workflow", "metrics", "comparison", "profiles"))
    evidence = ''.join('<p><strong>' + esc(item['label']) + ':</strong> ' + esc(item['text']) + '</p>' for item in metrics['items'])
    profile_blocks = []
    for item in profiles['items']:
        profile_blocks.append('<article class="tool-profile"><h4>' + esc(item['name']) + '</h4><p>' + esc(item['description']) + '</p>' +
            ''.join('<p><strong>' + label + ':</strong> ' + esc(item[key]) + '</p>' for key, label in
                    [('strength', 'Strength'), ('limitation', 'Review consideration'), ('report_ranking', 'Ranking in this report')]) + '</article>')
    return f"""<section class="tool-background">
<div class="background-part">
<h2 class="sample-title">{esc(data['title'])}</h2>
<p>{esc(data['introduction'])}</p>
<p class="note">{esc(data['scope_note'])}</p>
{title(workflow['title'])}<p>{esc(workflow['introduction'])}</p>{table(workflow)}
{title(metrics['title'])}{evidence}
</div>
<div class="background-part">
{title(comparison['title'])}<p>{esc(comparison['introduction'])}</p>{table(comparison)}
</div>
<div class="background-part">
{title(profiles['title'])}{''.join(profile_blocks)}
<p class="background-source">{esc(data['source'])}</p>
</div>
</section>"""


def load_report_info(config_path=None):
    if config_path is None:
        return {}
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("YAML configuration requires PyYAML: pip install PyYAML") from exc
    with Path(config_path).open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError("Report configuration must be a YAML mapping of field names to values")
    return data


def load_logo(logo_path=None):
    if logo_path is None:
        return ""
    path = Path(logo_path)
    mime = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}.get(path.suffix.lower())
    if mime is None:
        raise ValueError("Logo must be a PNG or JPEG image")
    return "data:" + mime + ";base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def report_header(info, logo_uri, generated_at, subtitle):
    def rows(fields):
        result = []
        for label, keys in fields:
            value = next((info[k] for k in keys if info.get(k) is not None and str(info[k]).strip()), None)
            if value is not None:
                result.append(f'<p><strong>{html.escape(label)}:</strong> {html.escape(str(value))}</p>')
        return "".join(result)

    project = rows([
        ("PI", ["PI_NAME"]), ("Project", ["PROJECT_NAME"]), ("Task", ["TASK_ID"]),
        ("Project Lead(s)", ["PROJECT_LEAD_NAME", "PPOJECT_LEAD_NAME"]),
        ("Department", ["DEPARTMENT"]), ("SRM", ["SRM"]), ("Reference genome", ["SPECIES"]),
    ])
    team = rows([
        ("Lead Analyst(s)", ["LEAD_ANALYSTS"]), ("Group Lead", ["GROUP_LEAD"]),
        ("Contact E-mail", ["CONTACT_EMAIL"]), ("CAB Pipeline", ["PIPELINE"]),
    ])
    dates = rows([("Date started", ["START_DATE"]), ("Date completed", ["COMPLETION_DATE"])])
    logo = f'<img class="cab-logo" src="{logo_uri}" alt="St. Jude Center for Applied Bioinformatics">' if logo_uri else ""
    team_block = f'<h2>CAB Analysis Team</h2><div class="team-details">{team}</div>' if team else ""
    cover_class = " report-cover" if info or logo_uri else ""
    return f"""<header class="header{cover_class}">
  {logo}
  <h1>RNA-Seq Fusion Analysis Report</h1>
  <p class="report-subtitle">{html.escape(subtitle)}</p>
  <div class="report-info">{project}{team_block}
    <div class="report-dates">{dates}<p><strong>Report generated:</strong> {html.escape(generated_at)}</p></div>
  </div>
</header>"""


# ---------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------

PATTERNS = {
    "arriba": [
        re.compile(r"^(?P<sample>.+?)[_-]arriba(?:\(\d+\))?\.tsv$", re.I),
    ],
    "fusioncatcher": [
        re.compile(r"^(?P<sample>.+?)[_-]fusion[_-]?catcher(?:\(\d+\))?\.(?:txt|tsv)$", re.I),
        re.compile(r"^(?P<sample>.+?)[_-]fusioncatcher(?:\(\d+\))?\.(?:txt|tsv)$", re.I),
    ],
    "starfusion": [
        re.compile(r"^(?P<sample>.+?)[_-]star[_-]?fusion(?:\(\d+\))?\.tsv$", re.I),
        re.compile(r"^(?P<sample>.+?)[_-]star-fusion(?:\(\d+\))?\.tsv$", re.I),
    ],
}


def discover_samples(input_dir: Path):
    samples = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        for tool, regexes in PATTERNS.items():
            for rx in regexes:
                m = rx.match(path.name)
                if m:
                    sample = m.group("sample").strip()
                    samples.setdefault(sample, {})[tool] = path
                    break
            else:
                continue
            break
    return samples


def read_table(path: Path):
    return pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False, low_memory=False
    )


def numeric(series):
    return pd.to_numeric(series, errors="coerce").fillna(0)


def join_unique(values):
    vals = []
    for v in values:
        v = str(v).strip()
        if v and v not in vals:
            vals.append(v)
    return "; ".join(vals)


# ---------------------------------------------------------------------
# Top-N extraction
# ---------------------------------------------------------------------

def top_arriba(path: Path, top_n=5):
    df = read_table(path)
    required = ["#gene1", "gene2", "confidence", "split_reads1", "split_reads2", "discordant_mates"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing Arriba columns: {missing}")

    df["Split_Reads"] = numeric(df["split_reads1"]) + numeric(df["split_reads2"])
    df["Discordant_Mates"] = numeric(df["discordant_mates"])

    rank = {"high": 3, "medium": 2, "low": 1}
    df["_confidence_rank"] = df["confidence"].str.lower().map(rank).fillna(0)
    df["Fusion"] = df["#gene1"].astype(str) + "--" + df["gene2"].astype(str)

    df = (
        df.sort_values(
            ["_confidence_rank", "Split_Reads", "Discordant_Mates"],
            ascending=[False, False, False]
        )
        .drop_duplicates("Fusion")
        .head(top_n)
        .copy()
    )

    cols = {
        "Fusion": df["Fusion"],
        "Confidence": df["confidence"],
        "Split_Reads": df["Split_Reads"].astype(int),
        "Discordant_Mates": df["Discordant_Mates"].astype(int),
    }
    if "reading_frame" in df.columns:
        cols["Reading_Frame"] = df["reading_frame"]
    if "type" in df.columns:
        cols["Type"] = df["type"]

    out = pd.DataFrame(cols)
    out.insert(0, "Tool", "Arriba")
    return out.reset_index(drop=True)


def top_fusioncatcher(path: Path, top_n=5):
    df = read_table(path)
    g1 = "Gene_1_symbol(5end_fusion_partner)"
    g2 = "Gene_2_symbol(3end_fusion_partner)"
    required = [g1, g2, "Spanning_pairs", "Spanning_unique_reads"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing FusionCatcher columns: {missing}")

    df["Fusion"] = df[g1].astype(str) + "--" + df[g2].astype(str)
    df["_pairs"] = numeric(df["Spanning_pairs"])
    df["_unique"] = numeric(df["Spanning_unique_reads"])

    agg = {"_pairs": "max", "_unique": "max"}
    if "Predicted_effect" in df.columns:
        agg["Predicted_effect"] = join_unique
    if "Fusion_description" in df.columns:
        agg["Fusion_description"] = join_unique

    g = (
        df.groupby("Fusion", as_index=False)
        .agg(agg)
        .sort_values(["_pairs", "_unique"], ascending=[False, False])
        .head(top_n)
        .copy()
    )

    cols = {
        "Fusion": g["Fusion"],
        "Spanning_Pairs": g["_pairs"].astype(int),
        "Unique_Reads": g["_unique"].astype(int),
    }
    if "Predicted_effect" in g.columns:
        cols["Predicted_Effect"] = g["Predicted_effect"]
    if "Fusion_description" in g.columns:
        cols["Description"] = g["Fusion_description"]

    out = pd.DataFrame(cols)
    out.insert(0, "Tool", "FusionCatcher")
    return out.reset_index(drop=True)


def top_starfusion(path: Path, top_n=5):
    df = read_table(path)
    required = ["FusionName", "JunctionReadCount", "SpanningFragCount"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing STAR-Fusion columns: {missing}")

    df["_junction"] = numeric(df["JunctionReadCount"])
    df["_spanning"] = numeric(df["SpanningFragCount"])
    df["_total"] = df["_junction"] + df["_spanning"]

    df = (
        df.sort_values(["_total", "_junction", "_spanning"], ascending=[False, False, False])
        .drop_duplicates("FusionName")
        .head(top_n)
        .copy()
    )

    cols = {
        "Fusion": df["FusionName"],
        "Junction_Reads": df["_junction"].astype(int),
        "Spanning_Fragments": df["_spanning"].astype(int),
    }
    if "FFPM" in df.columns:
        cols["FFPM"] = df["FFPM"]
    if "SpliceType" in df.columns:
        cols["Splice_Type"] = df["SpliceType"]

    out = pd.DataFrame(cols)
    out.insert(0, "Tool", "STAR-Fusion")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------
# Consensus
# ---------------------------------------------------------------------

def build_consensus(arriba, fc, star):
    mapping = {}
    for tool, df in [("Arriba", arriba), ("FusionCatcher", fc), ("STAR-Fusion", star)]:
        for fusion in df["Fusion"].astype(str):
            mapping.setdefault(fusion, []).append(tool)

    rows = []
    for fusion, tools in mapping.items():
        rows.append({
            "Fusion": fusion,
            "Number_of_Tools": len(tools),
            "Supporting_Tools": ", ".join(tools),
        })

    return (
        pd.DataFrame(rows)
        .sort_values(["Number_of_Tools", "Fusion"], ascending=[False, True])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def plot_arriba(df, top_n=5):
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    ax.bar(df["Fusion"], pd.to_numeric(df["Split_Reads"], errors="coerce").fillna(0))
    ax.set_ylabel("Split reads")
    ax.set_title(f"Arriba Top {top_n} Fusion Candidates")
    ax.tick_params(axis="x", rotation=30)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_fusioncatcher(df, top_n=5):
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    x = list(range(len(df)))
    width = 0.38
    pairs = pd.to_numeric(df["Spanning_Pairs"], errors="coerce").fillna(0)
    unique = pd.to_numeric(df["Unique_Reads"], errors="coerce").fillna(0)
    ax.bar([i - width/2 for i in x], pairs, width=width, label="Spanning pairs")
    ax.bar([i + width/2 for i in x], unique, width=width, label="Unique reads")
    ax.set_xticks(x)
    ax.set_xticklabels(df["Fusion"], rotation=30, ha="right")
    ax.set_ylabel("Read count")
    ax.set_title(f"FusionCatcher Top {top_n} Fusion Candidates")
    ax.legend()
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_starfusion(df, top_n=5):
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    x = list(range(len(df)))
    width = 0.38
    junction = pd.to_numeric(df["Junction_Reads"], errors="coerce").fillna(0)
    spanning = pd.to_numeric(df["Spanning_Fragments"], errors="coerce").fillna(0)
    ax.bar([i - width/2 for i in x], junction, width=width, label="Junction reads")
    ax.bar([i + width/2 for i in x], spanning, width=width, label="Spanning fragments")
    ax.set_xticks(x)
    ax.set_xticklabels(df["Fusion"], rotation=30, ha="right")
    ax.set_ylabel("Read count")
    ax.set_title(f"STAR-Fusion Top {top_n} Fusion Candidates")
    ax.legend()
    fig.tight_layout()
    return fig_to_base64(fig)


# ---------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------

CSS = """
:root{
  --cabcolor:#FF7F00;
  --cablight:#fff3e6;
  --dark:#2c2a25;
  --border:#ded8d1;
  --muted:#666666;
  --bg:#faf9f7;
}
*{box-sizing:border-box}
body{
  margin:0;
  font-family:Arial,Helvetica,sans-serif;
  background:var(--bg);
  color:#222;
}
.header{
  background:white;
  padding:28px 34px 18px;
  border-bottom:1px solid var(--border);
}
.header h1{
  margin:0;
  font-size:30px;
  color:var(--cabcolor);
}
.header p{
  margin:7px 0 0;
  color:var(--muted);
}
.tabs{
  display:flex;
  gap:4px;
  overflow-x:auto;
  padding:12px 24px 0;
  background:white;
  border-bottom:1px solid var(--border);
  position:sticky;
  top:0;
  z-index:10;
}
.tab-button{
  border:none;
  background:#f2ece5;
  padding:12px 18px;
  cursor:pointer;
  font-weight:700;
  border-radius:7px 7px 0 0;
  white-space:nowrap;
}
.tab-button.active{
  background:var(--cabcolor);
  color:#222;
}
.sheet{
  display:none;
  max-width:1200px;
  margin:24px auto;
  background:white;
  padding:30px;
  box-shadow:0 2px 10px rgba(0,0,0,.08);
}
.sheet.active{display:block}
.sample-title{
  font-size:25px;
  color:var(--cabcolor);
  margin:0 0 6px;
}
.section-title{
  color:var(--cabcolor);
  border-bottom:2px solid var(--cablight);
  padding-bottom:7px;
  margin-top:28px;
}
.cards{
  display:grid;
  grid-template-columns:repeat(3,minmax(0,1fr));
  gap:14px;
  margin:16px 0;
}
.card{
  background:var(--cablight);
  padding:16px;
  border-radius:8px;
}
.card strong{
  display:block;
  font-size:25px;
  color:#a65000;
}
.table-wrap{overflow-x:auto}
table{
  border-collapse:collapse;
  width:100%;
  margin:12px 0 22px;
  font-size:13px;
}
th{
  background:var(--cabcolor);
  color:#222;
  text-align:left;
  padding:9px;
  border:1px solid #e4b789;
}
td{
  padding:8px 9px;
  border:1px solid var(--border);
  vertical-align:top;
  overflow-wrap:anywhere;
}
tr:nth-child(even) td{background:#fcfaf7}
.chart{
  text-align:center;
  margin:8px 0 24px;
}
.chart img{
  max-width:100%;
  height:auto;
  border:1px solid var(--border);
  border-radius:5px;
}
.note{
  background:#fff8e6;
  border-left:4px solid #e0a800;
  padding:12px 14px;
}
.good{
  background:#fff3e6;
  border-left:4px solid #FF7F00;
  padding:12px 14px;
}
.report-cover{max-width:1200px;margin:24px auto;padding:36px;background:white}
.cab-logo{display:block;width:420px;max-width:100%;height:auto;margin:0 0 30px}
.report-cover h1{font-size:30px;margin-bottom:12px}
.report-info{margin-top:32px;line-height:1.5;color:#222}
.report-info p{margin:3px 0;color:#222;overflow-wrap:anywhere}
.report-info h2{font-size:18px;color:var(--cabcolor);margin:24px 0 10px}
.team-details{margin-left:24px}
.report-dates{margin-top:24px}
thead{display:table-header-group}
.tool-background{max-width:1200px;margin:24px auto;padding:30px;background:white;line-height:1.5}
.tool-background p{margin:10px 0}
.tool-background h4{color:#995000;font-size:17px;margin:18px 0 8px}
.background-source{font-size:11px;color:#666}
.igv-section{margin-top:32px}
.igv-snapshot{margin:20px 0;break-inside:avoid}
.igv-snapshot img{width:100%;height:auto}
.igv-snapshot figcaption{color:#995000;margin-bottom:8px}
.report-appendix{max-width:1200px;margin:24px auto;padding:30px;background:white;line-height:1.5}
.references li{margin-bottom:14px;overflow-wrap:anywhere}
a{color:#995000;text-decoration:underline;overflow-wrap:anywhere}
.sample-attachment{background:#fff3e6;padding:12px;border-left:4px solid #FF7F00}
@media(max-width:800px){
  .cards{grid-template-columns:1fr}
  .sheet{margin:0;padding:18px}
}
@page{
  size:A4;
  margin:18mm;
  @bottom-center{content:counter(page);font-size:9px;color:#666}
}
@media print{
  .tabs{display:none}
  .screen-only{display:none}
  .tool-background{margin:0;padding:0;font-size:10.5px;break-after:page}
  .tool-background .sample-title{font-size:22px}
  .tool-background .section-title{font-size:16px;margin-top:16px}
  .tool-background table{font-size:9px}
  .tool-background h4{font-size:13px;break-after:avoid}
  .background-part + .background-part{break-before:page}
  .background-part > .section-title:first-child{margin-top:0}
  .tool-profile{break-inside:avoid}
  .igv-section{break-before:page}
  .igv-snapshot{break-before:page}
  .igv-snapshot img{max-height:210mm;object-fit:contain}
  .report-appendix{break-before:page;margin:0;padding:0;font-size:11px}
  .report-appendix h3{break-after:avoid}
  .references li{break-inside:avoid}
  .report-cover{break-after:page;margin:0;padding:0;border:0}
  .cab-logo{width:105mm;margin-bottom:12mm}
  .report-cover h1{font-size:23px}
  .report-info{margin-top:16mm;font-size:12px}
  .report-info h2{font-size:14px}
  .sheet + .sheet{break-before:page}
  body{background:white}
  .header{
    padding:0 0 12px;
    border-bottom:1px solid var(--border);
  }
  .header h1{font-size:24px}
  .header.report-cover{border-bottom:0}
  .caller-title{break-before:page;margin-top:0}
  .section-title + p{break-after:avoid}
  .table-wrap{overflow:visible}
  th,td{padding:5px 6px}
  .sheet{
    display:block !important;
    box-shadow:none;
    margin:0;
    max-width:none;
    padding:0;
  }
  .section-title{
    break-after:avoid;
  }
  table{
    font-size:9.5px;
  }
  tr{
    break-inside:avoid;
  }
  .chart{
    break-inside:avoid;
  }
  .chart img{
    max-height:115mm;
  }
}
/* The side outline follows the floating table of contents in the CAB reference. */
.report-sidebar{display:none}
.skip-nav{position:absolute;left:-10000px;top:8px}
.skip-nav:focus{left:12px;z-index:100;background:white;padding:12px}
@media screen{
  body{background:white}
  .report-layout{display:grid;grid-template-columns:250px minmax(0,1fr);gap:30px;max-width:1440px;margin:0 auto;padding:24px}
  .report-main{min-width:0}
  .report-sidebar{display:block;position:sticky;top:24px;align-self:start;max-height:calc(100vh - 48px);overflow:auto;border:1px solid #ddd;border-radius:6px;background:white}
  .contents-disclosure summary{padding:13px 16px;font-weight:bold;color:#333;cursor:pointer;border-bottom:1px solid #eee}
  .report-sidebar ul{list-style:none;padding:5px 0;margin:0}
  .report-sidebar li{margin:0;padding:0}
  .report-sidebar li[hidden]{display:none}
  .report-sidebar a{display:block;padding:8px 14px;color:#333;text-decoration:none;font-size:13px;line-height:1.4;border-left:3px solid transparent}
  .report-sidebar .toc-depth-1 a{padding-left:25px;font-size:12px}
  .report-sidebar .toc-depth-2 a{padding-left:38px;font-size:12px}
  .report-sidebar a:hover{background:#fff3e6}
  .report-sidebar a:focus-visible{outline:2px solid #995000;outline-offset:-2px}
  .report-sidebar a.active{background:#ff7f00;color:#222;border-left-color:#995000}
  .report-main .header,.report-main .tool-background,.report-main .sheet,.report-main .report-appendix{max-width:none;margin:0;padding:26px 12px;box-shadow:none}
  .report-main .header{padding-top:10px}
  .report-main h2,.report-main h3,.report-main h4{scroll-margin-top:88px}
  #report-information{scroll-margin-top:24px}
  .report-main .tabs{padding-left:12px}
}
@media screen and (max-width:900px){
  .report-layout{display:block;padding:12px}
  .report-sidebar{position:relative;top:auto;max-height:none;width:100%;margin-bottom:18px}
  .report-sidebar nav{max-height:55vh;overflow:auto}
  .report-main .report-cover,.report-main .tool-background,.report-main .sheet,.report-main .report-appendix{padding-left:6px;padding-right:6px}
}
@media print{
  @page{size:255mm 297mm;margin:18mm 18mm 18mm 63mm}
  .skip-nav{display:none!important}
  .report-layout,.report-main{display:block;margin:0;padding:0;width:auto;max-width:none}
  .report-sidebar{display:block!important;position:fixed;top:0;left:-49mm;width:39mm;
    max-height:none;overflow:visible;border:1px solid #ddd;border-radius:5px;background:white}
  .contents-disclosure summary{display:block;padding:9px 10px;background:#ff7f00;color:#222;
    font-weight:bold;font-size:11px;border-bottom:1px solid #ddd}
  .report-sidebar nav{display:block;max-height:none;overflow:visible}
  .report-sidebar ul{list-style:none;padding:4px 0;margin:0}
  .report-sidebar li{margin:0;padding:0;break-inside:avoid}
  .report-sidebar li[hidden]{display:block}
  .report-sidebar a{display:block;color:#333;text-decoration:none;font-size:9px;line-height:1.3;padding:5px 8px;overflow-wrap:anywhere}
  .report-sidebar .toc-depth-1 a{padding-left:14px;font-size:8px}
  .report-sidebar .toc-depth-2 a{padding-left:21px;font-size:8px}
  .interactive-link{padding:9px 12px;background:#fff3e6;border-left:3px solid #ff7f00}
}

"""


JS = """
function openSample(sampleId, btn){
  const sheet=document.getElementById(sampleId);
  if(!sheet || !sheet.classList.contains('sheet')) return;
  document.querySelectorAll('.sheet').forEach(x=>x.classList.toggle('active',x===sheet));
  document.querySelectorAll('.tab-button').forEach(x=>{
    const active=x.dataset.sample===sampleId;
    x.classList.toggle('active',active);x.setAttribute('aria-pressed',String(active));
  });
  updateContents();
}
function updateContents(){
  const current=document.querySelector('.sheet.active');
  document.querySelectorAll('.report-sidebar li[data-sample]').forEach(item=>{
    item.hidden=!!(current && item.dataset.sample && item.dataset.sample!==current.id && item.dataset.level!=='2');
  });
}
(function(){
  const links=Array.from(document.querySelectorAll('.report-sidebar a[data-target]'));
  const outline=document.querySelector('.contents-disclosure');
  const media=window.matchMedia('(max-width: 900px)');
  function responsive(){if(outline) outline.open=!media.matches;}
  responsive();
  if(media.addEventListener) media.addEventListener('change',responsive);
  function showTarget(id,scroll){
    const target=document.getElementById(id);
    if(!target) return;
    const sheet=target.closest('.sheet');
    if(sheet) openSample(sheet.id);
    if(scroll) target.scrollIntoView({block:'start'});
  }
  function fromHash(){
    let id;try{id=decodeURIComponent(location.hash.slice(1));}catch(e){return;}
    if(id) showTarget(id,true);
  }
  links.forEach(link=>link.addEventListener('click',()=>{
    showTarget(link.dataset.target,false);
    if(media.matches && outline) outline.open=false;
  }));
  window.addEventListener('hashchange',fromHash);
  let scheduled=false;
  function highlight(){
    scheduled=false;
    const visible=links.filter(link=>!link.closest('li').hidden && document.getElementById(link.dataset.target)?.getClientRects().length);
    let active=visible[0];
    for(const link of visible){if(document.getElementById(link.dataset.target).getBoundingClientRect().top<=130)active=link;}
    links.forEach(link=>{link.classList.toggle('active',link===active);
      if(link===active)link.setAttribute('aria-current','location');else link.removeAttribute('aria-current');});
  }
  window.addEventListener('scroll',()=>{if(!scheduled){scheduled=true;requestAnimationFrame(highlight);}},{passive:true});
  window.addEventListener('resize',highlight);
  updateContents();fromHash();highlight();
})();
"""



def copy_report_attachment(source, output_dir):
    source = Path(source)
    destination = Path(output_dir) / "attachments" / source.name
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)
    return "attachments/" + quote(source.name)


def bib_text(value):
    # Decode common BibTeX accents without interpreting any commands.
    marks = {'"': '\u0308', "'": '\u0301', '`': '\u0300', '^': '\u0302', '~': '\u0303', 'c': '\u0327'}
    for command, mark in marks.items():
        value = re.sub(r"\\" + re.escape(command) + r"\s*\{?([A-Za-z])\}?",
                       lambda m: unicodedata.normalize("NFC", m.group(1) + mark), value)
    return value.replace("{", "").replace("}", "").replace("--", "–")


def parse_bib_entries(text):
    """Read braced/quoted BibTeX fields, including nested author accents."""
    entries = []
    start_pattern = re.compile(r"@\w+\s*\{\s*([^,]+),")
    pos = 0
    while match := start_pattern.search(text, pos):
        i = match.end()
        fields = {}
        while i < len(text):
            while i < len(text) and (text[i].isspace() or text[i] == ','):
                i += 1
            if i >= len(text) or text[i] == '}':
                i += 1
                break
            field = re.match(r"(\w+)\s*=\s*", text[i:])
            if not field:
                raise ValueError(f"Invalid BibTeX field near {text[i:i+50]!r}")
            key = field.group(1).lower()
            i += field.end()
            if i >= len(text):
                raise ValueError("Incomplete BibTeX entry")
            if text[i] == '{':
                depth = 1
                i += 1
                begin = i
                while i < len(text) and depth:
                    if text[i] == '\\':
                        i += 2
                        continue
                    if text[i] == '{': depth += 1
                    elif text[i] == '}': depth -= 1
                    i += 1
                if depth:
                    raise ValueError("Unclosed BibTeX field")
                value = text[begin:i-1]
            elif text[i] == '"':
                i += 1
                begin = i
                while i < len(text) and text[i] != '"':
                    i += 2 if text[i] == '\\' else 1
                if i >= len(text):
                    raise ValueError("Unclosed quoted BibTeX field")
                value = text[begin:i]
                i += 1
            else:
                begin = i
                while i < len(text) and text[i] not in ',}': i += 1
                value = text[begin:i].strip()
            fields[key] = bib_text(value)
        entries.append(fields)
        pos = i
    return entries


def resource_markdown(text):
    """Render the supplied resource prose with escaped text and HTTP(S) links."""
    def inline(value):
        parts = []
        pos = 0
        for match in re.finditer(r"\[([^\]]+)\]\(([^)]+)\)", value):
            parts.append(html.escape(value[pos:match.start()]))
            label, url = match.groups()
            if url.startswith("ttps://"):
                url = "h" + url
            if urlparse(url).scheme in ("http", "https"):
                parts.append(f'<a href="{html.escape(url, quote=True)}">{html.escape(label)}</a>')
            else:
                parts.append(html.escape(label))
            pos = match.end()
        parts.append(html.escape(value[pos:]))
        return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", "".join(parts))
    text = text.replace(r"\pagebreak", "").replace("<br/>", "\n\n")
    blocks, paragraph = [], []
    def flush():
        if paragraph:
            blocks.append("<p>" + inline(" ".join(paragraph)) + "</p>")
            paragraph.clear()
    for line in text.splitlines():
        line = line.strip().rstrip('\\').strip()
        if not line:
            flush()
        elif line.startswith('#'):
            flush()
            title = line.lstrip('#').strip()
            if title.lower() != "references":
                blocks.append(f'<h3 class="section-title">{html.escape(title)}</h3>')
        else:
            paragraph.append(line)
    flush()
    return "\n".join(blocks)


def load_references(path, output_dir):
    if path is None:
        return ""
    text = Path(path).read_text(encoding="utf-8-sig")
    entries = parse_bib_entries(text)
    if not entries:
        raise ValueError(f"No BibTeX entries found in {path}")
    preamble = text[:re.search(r'@\w+\s*\{', text).start()]
    items = []
    for entry in entries:
        authors = entry.get('author', '').replace(' and ', '; ')
        citation = '. '.join(v for v in [authors, entry.get('year'), entry.get('title'),
                            entry.get('journal'), entry.get('volume'), entry.get('pages')] if v)
        doi = entry.get('doi', '')
        link = f' <a href="https://doi.org/{quote(doi, safe="/")}">doi:{html.escape(doi)}</a>' if doi else ''
        items.append(f'<li>{html.escape(citation)}.{link}</li>')
    source_link = copy_report_attachment(path, output_dir)
    return f"""<section class="report-appendix">
{resource_markdown(preamble)}
<h3 class="section-title">References</h3>
<ol class="references">{''.join(items)}</ol>
<p class="screen-only"><a href="{source_link}" download>Download {html.escape(Path(path).name)}</a></p>
</section>"""


def portable_pdf_links(pdf_path):
    """Keep sample PDF and interactive HTML links relative to the report folder."""
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import DictionaryObject, NameObject, TextStringObject, ArrayObject, NumberObject
    reader = PdfReader(str(pdf_path))
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    for page in writer.pages:
        for annotation in page.get('/Annots', []):
            annotation = annotation.get_object()
            action = annotation.get('/A')
            if not action or action.get('/S') != '/URI':
                continue
            url = urlparse(str(action.get('/URI', '')))
            if url.scheme != 'file':
                continue
            target = Path(unquote(url.path)).resolve()
            attachments = (pdf_path.parent / 'attachments').resolve()
            if target.parent == attachments and target.suffix.lower() == '.pdf':
                annotation[NameObject('/A')] = DictionaryObject({
                    NameObject('/S'): NameObject('/GoToR'),
                    NameObject('/F'): TextStringObject('attachments/' + target.name),
                    NameObject('/D'): ArrayObject([NumberObject(0), NameObject('/Fit')]),
                })
            elif target.name == 'interactive.html' and (pdf_path.parent / 'igv').resolve() in target.parents:
                relative = target.relative_to(pdf_path.parent.resolve()).as_posix()
                annotation[NameObject('/A')] = DictionaryObject({
                    NameObject('/S'): NameObject('/Launch'),
                    NameObject('/F'): TextStringObject(relative),
                })
    with pdf_path.open('wb') as handle:
        writer.write(handle)


def dataframe_html(df):
    return df.to_html(index=False, border=0, classes="data-table", escape=True)


def safe_id(text):
    return "sample_" + re.sub(r"[^A-Za-z0-9_-]+", "_", text)


def sample_sheet_html(sample, arriba, fc, star, consensus, imgs, top_n=5, sample_pdf_link="", igv_html=""):
    attachment = (f'<p class="sample-attachment"><strong>Supporting PDF:</strong> '
                  f'<a href="{html.escape(sample_pdf_link, quote=True)}">Open {html.escape(sample)}.pdf</a></p>') if sample_pdf_link else ""
    multi = consensus[consensus["Number_of_Tools"] >= 2]
    multi_count = len(multi)
    three_count = int((consensus["Number_of_Tools"] == 3).sum()) if not consensus.empty else 0

    consensus_msg = (
        f"{multi_count} fusion(s) appear in the Top {top_n} of at least two callers; "
        f"{three_count} appear in all three."
        if not consensus.empty else
        f"No cross-tool Top-{top_n} consensus candidates were identified."
    )

    return f"""
<section id="{safe_id(sample)}" class="sheet">
  <h2 class="sample-title">{html.escape(sample)} — RNA-Seq Fusion Analysis</h2>
  <p>Top {top_n} fusion candidates from Arriba, FusionCatcher, and STAR-Fusion.</p>

  {attachment}

  <div class="cards">
    <div class="card"><strong>{len(arriba)}</strong>Arriba Top candidates</div>
    <div class="card"><strong>{len(fc)}</strong>FusionCatcher Top candidates</div>
    <div class="card"><strong>{len(star)}</strong>STAR-Fusion Top candidates</div>
  </div>

  <h3 class="section-title">Cross-Tool Consensus</h3>
  <div class="good">{html.escape(consensus_msg)}</div>
  <div class="table-wrap">{dataframe_html(consensus)}</div>

  <h3 class="section-title caller-title">Arriba — Top {top_n}</h3>
  <p>Ranked by confidence, followed by split-read and discordant-mate support.</p>
  <div class="table-wrap">{dataframe_html(arriba.drop(columns=["Tool"]))}</div>
  <div class="chart"><img src="data:image/png;base64,{imgs['arriba']}" alt="Arriba top {top_n} chart"></div>

  <h3 class="section-title caller-title">FusionCatcher — Top {top_n}</h3>
  <p>Ranked by spanning-pair support followed by unique spanning reads.</p>
  <div class="table-wrap">{dataframe_html(fc.drop(columns=["Tool"]))}</div>
  <div class="chart"><img src="data:image/png;base64,{imgs['fusioncatcher']}" alt="FusionCatcher top {top_n} chart"></div>

  <h3 class="section-title caller-title">STAR-Fusion — Top {top_n}</h3>
  <p>Ranked by total junction-read and spanning-fragment support.</p>
  <div class="table-wrap">{dataframe_html(star.drop(columns=["Tool"]))}</div>
  <div class="chart"><img src="data:image/png;base64,{imgs['starfusion']}" alt="STAR-Fusion top {top_n} chart"></div>

  <h3 class="section-title">Interpretation</h3>
  <p>
    Multi-caller candidates should generally be reviewed first, followed by strongly supported
    single-caller findings. Breakpoints, original alignments, reading frame, recurrent-normal
    annotations, read-through status, and other caller-specific artifact flags should be reviewed
    before biological interpretation.
  </p>
  <div class="note">
    This report summarizes computational fusion-caller output only. It does not establish
    pathogenicity, diagnosis, treatment eligibility, or clinical significance.
  </div>
{igv_html}
</section>
"""


def add_side_navigation(document):
    """Add static, escaped anchors and a shared TOC without external libraries."""
    class OutlineParser(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=False)
            self.parts, self.entries, self.sections = [], [], []
            self.in_header = False
            self.heading = None

        def handle_starttag(self, tag, attrs):
            values = dict(attrs)
            raw = self.get_starttag_text()
            if tag == 'header':
                self.in_header = True
                raw = raw[:-1] + ' id="report-information">'
            if tag == 'section':
                parent = self.sections[-1] if self.sections else ''
                self.sections.append(values.get('id', '') if 'sheet' in values.get('class', '').split() else parent)
            if tag in ('h2', 'h3', 'h4') and not self.in_header:
                anchor = f"report-section-{len(self.entries) + 1}"
                sample = self.sections[-1] if self.sections else ''
                self.heading = dict(id=anchor, sample=sample, level=int(tag[1]), text=[])
                self.entries.append(self.heading)
                raw = raw[:-1] + f' id="{anchor}" tabindex="-1">'
            self.parts.append(raw)

        def handle_endtag(self, tag):
            self.parts.append(f"</{tag}>")
            if tag == 'header': self.in_header = False
            if tag == 'section' and self.sections: self.sections.pop()
            if tag in ('h2', 'h3', 'h4'): self.heading = None

        def handle_data(self, data):
            self.parts.append(data)
            if self.heading is not None: self.heading['text'].append(data)

        def handle_entityref(self, name):
            token = '&' + name + ';'
            self.parts.append(token)
            if self.heading is not None: self.heading['text'].append(html.unescape(token))

        def handle_charref(self, name):
            token = '&#' + name + ';'
            self.parts.append(token)
            if self.heading is not None: self.heading['text'].append(html.unescape(token))

        def handle_comment(self, data): self.parts.append('<!--' + data + '-->')
        def handle_decl(self, decl): self.parts.append('<!' + decl + '>')
        def handle_startendtag(self, tag, attrs): self.parts.append(self.get_starttag_text())

    parser = OutlineParser()
    parser.feed(document)
    links = ['<li><a href="#report-information" data-target="report-information">Report information</a></li>']
    for item in parser.entries:
        title = ''.join(item['text']).strip()
        if title == 'Fusion analysis overview': title = 'Overview'
        if item['sample'] and item['level'] == 2:
            title = title.split(' — ')[0] + ' results'
        depth = 0 if item['level'] == 2 or not item['sample'] and title in ('Additional Resources', 'Acknowledgements', 'References') else item['level'] - 2
        links.append(f'<li class="toc-depth-{depth}" data-sample="{html.escape(item["sample"], quote=True)}" data-level="{item["level"]}">'
                     f'<a href="#{item["id"]}" data-target="{item["id"]}">{html.escape(title)}</a></li>')
    sidebar = ('<a class="skip-nav" href="#report-main">Skip to report</a><div class="report-layout">'
               '<aside class="report-sidebar"><details class="contents-disclosure" open><summary>Contents</summary>'
               '<nav aria-label="Report contents"><ul>' + ''.join(links) + '</ul></nav></details></aside>'
               '<main class="report-main" id="report-main">')
    rendered = ''.join(parser.parts)
    rendered = rendered.replace('<body>', '<body>' + sidebar, 1)
    rendered = rendered.replace('</body>', '</main></div><script>' + JS + '</script></body>', 1)
    return rendered


def build_full_html(sample_reports, report_info=None, logo_uri="", generated_at="", references_html="", tool_background_html=""):
    header = report_header(report_info or {}, logo_uri, generated_at, "All samples | RNA-Seq fusion calling")
    buttons = []
    sheets = []

    for i, r in enumerate(sample_reports):
        active = " active" if i == 0 else ""
        sid = safe_id(r["sample"])
        buttons.append(
            f'<button class="tab-button{active}" data-sample="{sid}" aria-pressed="{str(i == 0).lower()}" onclick="openSample(\'{sid}\', this)">'
            f'{html.escape(r["sample"])}</button>'
        )

        sheet = sample_sheet_html(
            r["sample"], r["arriba"], r["fc"], r["star"], r["consensus"], r["imgs"], r.get("top_n", 5), r.get("sample_pdf_link", ""), r.get("igv_html", "")
        )
        if i == 0:
            sheet = sheet.replace('class="sheet"', 'class="sheet active"', 1)
        sheets.append(sheet)

    return add_side_navigation(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RNA-Seq Fusion Analysis — All Samples</title>
<style>{CSS}</style>
</head>
<body>
{header}
{tool_background_html}
<div class="tabs">
  {''.join(buttons)}
</div>
{''.join(sheets)}
{references_html}
</body>
</html>
""")


def build_individual_html(sample_report):
    r = sample_report
    header = report_header(r.get("report_info", {}), r.get("logo_uri", ""), r.get("generated_at", ""),
                           f"Sample: {r['sample']} | Top {r.get('top_n', 5)} fusion pairs per tool")
    sheet = sample_sheet_html(
        r["sample"], r["arriba"], r["fc"], r["star"], r["consensus"], r["imgs"], r.get("top_n", 5), r.get("sample_pdf_link", ""), r.get("igv_html", "")
    ).replace('class="sheet"', 'class="sheet active"', 1)

    return add_side_navigation(f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(r["sample"])} Fusion Report</title>
<style>{CSS}</style>
</head>
<body>
{header}
{r.get("tool_background_html", "")}
{sheet}
{r.get("references_html", "")}
</body>
</html>
""")



# ---------------------------------------------------------------------
# Individual HTML -> PDF
# ---------------------------------------------------------------------

def html_to_pdf(html_path: Path, pdf_path: Path):
    """
    Convert an individual sample HTML report to PDF using WeasyPrint.
    Returns (True, "") on success or (False, error_message) on failure.
    """
    if not WEASYPRINT_AVAILABLE:
        return False, (
            "WeasyPrint is unavailable. Install weasyprint and its system libraries. "
            f"Import error: {WEASYPRINT_IMPORT_ERROR}"
        )

    try:
        WeasyHTML(
            filename=str(html_path),
            base_url=str(html_path.parent)
        ).write_pdf(str(pdf_path))
        portable_pdf_links(pdf_path)
        return True, ""
    except Exception as exc:
        return False, str(exc)


# ---------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------

def process_sample(sample, files, output_dir, top_n=5, report_info=None, logo_uri="", generated_at="", sample_pdf_link="", references_html="", igv_options=None, tool_background_html=""):
    required = {"arriba", "fusioncatcher", "starfusion"}
    missing = required - set(files)
    if missing:
        raise RuntimeError(f"Missing caller file(s): {', '.join(sorted(missing))}")

    arriba = top_arriba(files["arriba"], top_n)
    fc = top_fusioncatcher(files["fusioncatcher"], top_n)
    star = top_starfusion(files["starfusion"], top_n)
    consensus = build_consensus(arriba, fc, star)

    imgs = {
        "arriba": plot_arriba(arriba, top_n),
        "fusioncatcher": plot_fusioncatcher(fc, top_n),
        "starfusion": plot_starfusion(star, top_n),
    }

    # All final per-sample files are written directly into outdir/
    combined = pd.concat([arriba, fc, star], ignore_index=True, sort=False)
    combined.insert(0, "Sample", sample)

    top_tsv = output_dir / f"{sample}_top{top_n}_summary.tsv"
    combined.to_csv(top_tsv, sep="\t", index=False)

    report = {
        "sample": sample,
        "top_n": top_n,
        "report_info": report_info or {},
        "logo_uri": logo_uri,
        "generated_at": generated_at,
        "sample_pdf_link": sample_pdf_link,
        "references_html": references_html,
        "tool_background_html": tool_background_html,
        "arriba": arriba,
        "fc": fc,
        "star": star,
        "consensus": consensus,
        "imgs": imgs,
        "combined": combined,
    }

    if igv_options is not None:
        igv_result = build_igv(sample, files, output_dir,
                              selectors={"arriba": top_arriba, "fusioncatcher": top_fusioncatcher, "starfusion": top_starfusion},
                              read_table=read_table, numeric=numeric, **igv_options)
        report["igv_html"] = igv_report_html(igv_result)
    individual_html = build_individual_html(report)
    html_path = output_dir / f"{sample}_fusion_report.html"
    html_path.write_text(individual_html, encoding="utf-8")

    pdf_path = output_dir / f"{sample}_fusion_report.pdf"
    pdf_ok, pdf_error = html_to_pdf(html_path, pdf_path)

    if not pdf_ok and pdf_path.exists():
        pdf_path.unlink()

    report["html_path"] = html_path
    report["pdf_path"] = pdf_path if pdf_ok else None
    report["pdf_error"] = pdf_error
    report["top_tsv"] = top_tsv

    return report


def run_batch(input_dir, output_dir, top_n=5, config_path=None, logo_path=None, references_path=None, igv_top_n=5, bam_dir="bam", igv_genome=None, igv_flank=250, bam_pattern=None, igv_command=None, tool_background_path=None):
    if top_n not in (5, 10):
        raise ValueError("top_n must be 5 or 10")
    tool_background_html = load_tool_background(tool_background_path)
    report_info = load_report_info(config_path)
    logo_uri = load_logo(logo_path)
    igv_options = dict(top_n=igv_top_n, genome=igv_genome or str(report_info.get("SPECIES", "hg38")),
                       bam_dir=bam_dir, flank=igv_flank, bam_pattern=bam_pattern, command=igv_command)
    if igv_top_n not in (5, 10) or igv_flank < 1:
        raise ValueError("--igv-top-n must be 5 or 10 and --igv-flank must be positive")
    generated_at = datetime.now().astimezone().strftime("%H:%M:%S %Z %m/%d/%Y")
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if references_path is None:
        candidate = input_dir / "Fusion_reference.bib"
        references_path = candidate if candidate.is_file() else None
    references_html = load_references(references_path, output_dir)
    samples = discover_samples(input_dir)
    if not samples:
        raise RuntimeError(f"No recognizable fusion caller files found in {input_dir}")

    successful = []
    status = []

    print(f"Found {len(samples)} sample(s): {', '.join(sorted(samples))}")

    for sample in sorted(samples):
        try:
            sample_pdf = input_dir / f"{sample}.pdf"
            sample_pdf_link = copy_report_attachment(sample_pdf, output_dir) if sample_pdf.is_file() else ""
            report = process_sample(sample, samples[sample], output_dir, top_n, report_info, logo_uri, generated_at,
                                    sample_pdf_link, references_html, igv_options, tool_background_html)
            successful.append(report)

            pdf_status = "OK" if report["pdf_path"] else "FAILED"

            status.append({
                "Sample": sample,
                "Status": "OK",
                "HTML_Report": str(report["html_path"]),
                "PDF_Report": str(report["pdf_path"]) if report["pdf_path"] else "",
                "PDF_Status": pdf_status,
                f"Top{top_n}_Summary": str(report["top_tsv"]),
                "Error": report["pdf_error"],
            })

            print(f"[OK] {sample}")
            print(f"     HTML: {report['html_path']}")
            if report["pdf_path"]:
                print(f"     PDF:  {report['pdf_path']}")
            else:
                print(f"     PDF:  not generated - {report['pdf_error']}")

        except Exception as exc:
            status.append({
                "Sample": sample,
                "Status": "FAILED",
                "HTML_Report": "",
                "PDF_Report": "",
                "PDF_Status": "FAILED",
                f"Top{top_n}_Summary": "",
                "Error": str(exc),
            })
            print(f"[FAILED] {sample}: {exc}")

    if successful:
        # Combined tabbed HTML report.
        combined_html = build_full_html(successful, report_info, logo_uri, generated_at, references_html, tool_background_html)
        main_html = output_dir / "fusion_report_all_samples.html"
        main_html.write_text(combined_html, encoding="utf-8")

        # Combined Top-N data table.
        all_summary = pd.concat(
            [r["combined"] for r in successful],
            ignore_index=True,
            sort=False
        )
        all_summary.to_csv(
            output_dir / f"all_samples_top{top_n}_summary.tsv",
            sep="\t",
            index=False
        )

        print(f"\nCombined HTML report: {main_html}")

    pd.DataFrame(status).to_csv(
        output_dir / "batch_run_summary.tsv",
        sep="\t",
        index=False
    )

    print(f"Batch summary: {output_dir / 'batch_run_summary.tsv'}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate individual HTML + PDF RNA-Seq fusion reports and a combined tabbed HTML report for multiple samples."
    )
    parser.add_argument("-i", "--inputs", default="inputs", help="Input folder")
    parser.add_argument("-o", "--outdir", default="results", help="Output folder")
    parser.add_argument(
        "--top-n", type=int, choices=[5, 10], default=5,
        help="Maximum fusion pairs per tool per sample (5 or 10; default: 5)",
    )
    parser.add_argument("--config", type=Path, help="YAML file containing project and analysis-team information")
    parser.add_argument("--logo", type=Path, help="CAB logo image (PNG or JPEG), embedded in reports")
    parser.add_argument("--references", type=Path, help="Reference file (default: Fusion_reference.bib in the input folder)")
    parser.add_argument("--igv-top-n", type=int, choices=[5, 10], default=5, help="IGV events per caller (independent of --top-n; default: 5)")
    parser.add_argument("--bam-dir", type=Path, default=Path("bam"), help="Folder of sample BAMs and indexes (default: bam)")
    parser.add_argument("--bam-pattern", help="Optional BAM glob relative to --bam-dir; use {sample}, e.g. '{sample}/*.sorted.bam'")
    parser.add_argument("--igv-genome", help="IGV genome ID (default: SPECIES from YAML, otherwise hg38)")
    parser.add_argument("--igv-flank", type=int, default=250, help="Bases flanking each breakpoint in IGV (default: 250)")
    parser.add_argument("--igv-command", help="IGV Desktop launcher command; runs exported batch files and embeds screenshots. Example: 'xvfb-run -a /opt/IGV/igv.sh'")
    parser.add_argument("--tool-background", type=Path, help="Shared background JSON (default: fusion_tool_background.json beside the script)")
    args = parser.parse_args()

    run_batch(args.inputs, args.outdir, args.top_n, args.config, args.logo, args.references,
              args.igv_top_n, args.bam_dir, args.igv_genome, args.igv_flank, args.bam_pattern, args.igv_command, args.tool_background)


if __name__ == "__main__":
    main()
