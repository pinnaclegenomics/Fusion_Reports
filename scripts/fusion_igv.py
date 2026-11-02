"""IGV exports and report views for the fusion report generator.

No read images are synthesized. Screenshots must come from IGV Desktop.
BAMs stay in their original location and are never copied into the report.
"""
import base64
import csv
import gzip
import hashlib
import html
import json
import re
import shlex
import struct
import subprocess
from pathlib import Path
from urllib.parse import quote

CALLERS = {'arriba': 'Arriba', 'fusioncatcher': 'FusionCatcher', 'starfusion': 'STAR-Fusion'}


def parse_breakpoint(value):
    match = re.fullmatch(r'([^:\s]+):(\d+)(?::([+-]))?', str(value).strip())
    if not match or int(match[2]) < 1:
        raise ValueError(f'Invalid breakpoint: {value!r}')
    return match[1], int(match[2]), match[3] or '.'


def select_events(files, top_n, read_table, selectors, numeric):
    events, errors = [], []
    for tool, label in CALLERS.items():
        df = read_table(files[tool])
        if tool == 'arriba':
            df['_fusion'] = df['#gene1'] + '--' + df['gene2']
            df['_rank'] = df['confidence'].str.lower().map({'high': 3, 'medium': 2, 'low': 1}).fillna(0)
            df['_a'] = numeric(df['split_reads1']) + numeric(df['split_reads2'])
            df['_b'] = numeric(df['discordant_mates'])
            df = df.sort_values(['_rank', '_a', '_b'], ascending=False)
            columns = ('breakpoint1', 'breakpoint2')
        elif tool == 'fusioncatcher':
            df['_fusion'] = df['Gene_1_symbol(5end_fusion_partner)'] + '--' + df['Gene_2_symbol(3end_fusion_partner)']
            df['_a'] = numeric(df['Spanning_pairs'])
            df['_b'] = numeric(df['Spanning_unique_reads'])
            df = df.sort_values(['_a', '_b'], ascending=False)
            columns = ('Fusion_point_for_gene_1(5end_fusion_partner)', 'Fusion_point_for_gene_2(3end_fusion_partner)')
        else:
            df['_fusion'] = df['FusionName']
            df['_a'] = numeric(df['JunctionReadCount'])
            df['_b'] = numeric(df['SpanningFragCount'])
            df['_total'] = df['_a'] + df['_b']
            df = df.sort_values(['_total', '_a', '_b'], ascending=False)
            columns = ('LeftBreakpoint', 'RightBreakpoint')
        selected = selectors[tool](files[tool], top_n)
        for rank, fusion in enumerate(selected['Fusion'], 1):
            try:
                row = df[df['_fusion'] == fusion].iloc[0]
                left = parse_breakpoint(row[columns[0]])
                right = parse_breakpoint(row[columns[1]])
                if tool == 'arriba':
                    strands = []
                    for col in ('strand1(gene/fusion)', 'strand2(gene/fusion)'):
                        value = str(row.get(col, '.')).split('/')[-1]
                        strands.append(value if value in ('+', '-') else '.')
                    left = (*left[:2], strands[0]); right = (*right[:2], strands[1])
                events.append(dict(tool=label, tool_key=tool, rank=rank, fusion=fusion,
                                   chrom1=left[0], pos1=left[1], strand1=left[2],
                                   chrom2=right[0], pos2=right[1], strand2=right[2]))
            except (KeyError, ValueError, IndexError) as exc:
                errors.append(f'{label} rank {rank} ({fusion}): {exc}')
    return events, errors


def find_bam(sample, bam_dir, pattern=None):
    root = Path(bam_dir)
    if pattern:
        candidates = sorted(root.glob(pattern.format(sample=sample)))
    else:
        candidates = []
        exact = root / f'{sample}.bam'
        if exact.is_file():
            candidates = [exact]
        else:
            candidates.extend((root / sample).glob('*.bam'))
            candidates.extend(root.glob(f'{sample}.*.bam'))
            candidates.extend(root.glob(f'{sample}_*.bam'))
    candidates = sorted({p.resolve() for p in candidates if p.is_file() and p.suffix.lower() == '.bam'})
    if len(candidates) > 1:
        return None, None, 'Multiple BAMs match this sample. Set --bam-pattern to choose one.'
    if not candidates:
        return None, None, f'BAM not found for {sample} in {root}.'
    bam = candidates[0]
    indexes = [Path(str(bam)+'.bai'), bam.with_suffix('.bai'), Path(str(bam)+'.csi'), bam.with_suffix('.csi')]
    index = next((p for p in indexes if p.is_file()), None)
    if index is None:
        return bam, None, f'BAM index (.bai or .csi) not found for {bam.name}.'
    return bam, index, ''


def bam_references(path):
    """Read only the BAM header, including contig lengths, via its BGZF stream."""
    def integer(handle):
        return struct.unpack('<i', handle.read(4))[0]
    with gzip.open(path, 'rb') as handle:
        if handle.read(4) != b'BAM\x01':
            raise ValueError('File does not have a BAM header')
        length = integer(handle)
        if not 0 <= length <= 100_000_000:
            raise ValueError('Invalid BAM header length')
        header = handle.read(length).decode('utf-8', errors='replace')
        if 'SO:coordinate' not in header.split('\n')[0]:
            raise ValueError('BAM header must declare coordinate sorting (SO:coordinate)')
        n_ref = integer(handle)
        if not 0 <= n_ref <= 1_000_000:
            raise ValueError('Invalid BAM reference count')
        refs = {}
        for _ in range(n_ref):
            n = integer(handle)
            if not 1 <= n <= 1_000_000: raise ValueError('Invalid BAM reference name')
            name = handle.read(n).rstrip(b'\0').decode()
            refs[name] = integer(handle)
        return refs


def match_contig(chrom, refs):
    if chrom in refs: return chrom
    alternatives = [chrom[3:] if chrom.startswith('chr') else 'chr'+chrom]
    if chrom in ('M', 'MT', 'chrM', 'chrMT'):
        alternatives += ['M', 'MT', 'chrM', 'chrMT']
    matches = list(dict.fromkeys(c for c in alternatives if c in refs))
    if len(matches) != 1: raise ValueError(f'Cannot match contig {chrom} to BAM header')
    return matches[0]


def quoted_command(value):
    text = str(value)
    if any(c in text for c in ('\n', '\r', '"')):
        raise ValueError('IGV command values cannot contain quotes or newlines')
    return '"'+text+'"'


def write_viewer(path, sample, genome, events):
    # JSON is escaped for embedding in a script, not evaluated as source from input.
    payload = json.dumps(dict(sample=sample, genome=genome, events=events)).replace('<', '\\u003c')
    template = '''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Interactive fusion IGV</title><style>
body{font:15px Arial;margin:20px;color:#222}h1{color:#ff7f00;font-size:24px}label{display:block;margin:12px 0}
button,select{padding:9px;max-width:100%}button{background:#ff7f00;border:1px solid #995000;cursor:pointer}
#status{padding:12px;background:#fff3e6;white-space:pre-wrap}#igv{margin-top:16px}p{line-height:1.5}
</style><script src="https://cdn.jsdelivr.net/npm/igv@3.8.7/dist/igv.min.js"></script></head><body>
<h1 id="title"></h1><p>Select this sample's coordinate-sorted BAM and its matching BAI/CSI index.
BAM data are read locally by your browser. Internet is needed for IGV and the reference genome.</p>
<label>BAM file <input id="bam" type="file" accept=".bam"></label>
<label>BAM index <input id="bai" type="file" accept=".bai,.csi"></label>
<button id="load">Load sample alignments</button>
<label>Fusion event <select id="event"></select></label>
<p id="status">Select files to load read alignments. The event list contains caller-reported breakpoint pairs.</p>
<div id="igv"></div><script>
const data=__DATA__;let browser=null;let loading=false;
document.getElementById('title').textContent=data.sample+' — Interactive IGV';
const selector=document.getElementById('event'),status=document.getElementById('status');
data.events.forEach((e,i)=>{const option=document.createElement('option');option.value=i;
option.textContent=e.tool+' #'+e.rank+': '+e.fusion+' | '+e.locus;selector.appendChild(option)});
async function selected(){if(browser){await browser.search(data.events[Number(selector.value)].locus)}}
selector.addEventListener('change',()=>selected().catch(e=>{status.textContent=e.message}));
document.getElementById('load').addEventListener('click',async()=>{
if(loading)return;const bam=document.getElementById('bam').files[0],bai=document.getElementById('bai').files[0];
if(!bam||!bai){status.textContent='Select both the BAM and its matching index.';return}
if(!data.events.length){status.textContent='No valid breakpoint pairs were exported.';return}
loading=true;status.textContent='Loading '+bam.name+' …';
try{if(typeof igv==='undefined')throw new Error('IGV could not load. Check the internet connection and reload.');
if(browser){igv.removeBrowser(browser);browser=null}
browser=await igv.createBrowser(document.getElementById('igv'),{genome:data.genome,
locus:data.events[Number(selector.value)].locus,tracks:[{name:data.sample+' reads',type:'alignment',format:'bam',
url:bam,indexURL:bai,height:500,viewAsPairs:true,colorBy:'pairOrientation'}]});
status.textContent='Loaded '+bam.name+'. Choose an event to inspect both breakpoints. Genome: '+data.genome;
}catch(e){status.textContent='Unable to load alignments: '+e.message}finally{loading=false}
});
</script></body></html>'''
    path.write_text(template.replace('__DATA__', payload), encoding='utf-8')


def build_igv(sample, files, output_dir, top_n, genome, bam_dir, flank, selectors, read_table,
              numeric, bam_pattern=None, command=None):
    if top_n not in (5,10): raise ValueError('IGV top_n must be 5 or 10')
    if flank < 1: raise ValueError('IGV flank must be positive')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', genome):
        raise ValueError('--igv-genome must be an IGV genome ID, for example hg38 or hg19')
    events, errors = select_events(files, top_n, read_table, selectors, numeric)
    bam, index, bam_error = find_bam(sample, bam_dir, bam_pattern)
    refs = {}
    if bam and index:
        try: refs = bam_references(bam)
        except (OSError, ValueError, struct.error) as exc: bam_error = f'BAM validation failed: {exc}'
    accepted = []
    for event in events:
        try:
            loci = []
            for n in (1,2):
                chrom, pos = event[f'chrom{n}'], event[f'pos{n}']
                if refs:
                    chrom = match_contig(chrom, refs)
                    if pos > refs[chrom]: raise ValueError(f'Breakpoint exceeds {chrom} length')
                    event[f'chrom{n}'] = chrom
                end = min(pos+flank, refs[chrom]) if refs else pos+flank
                loci.append(f'{chrom}:{max(1,pos-flank)}-{end}')
            event['locus'] = ' '.join(loci)
            accepted.append(event)
        except ValueError as exc:
            errors.append(f"{event['tool']} #{event['rank']} {event['fusion']}: {exc}")
    events = accepted
    signature = {'events':events,'genome':genome,'flank':flank,'top_n':top_n,
                 'bam':str(bam) if bam else None,
                 'bam_stat':(bam.stat().st_size,bam.stat().st_mtime_ns) if bam else None,
                 'index_stat':(index.stat().st_size,index.stat().st_mtime_ns) if index else None}
    digest = hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    safe_sample = re.sub(r'[^A-Za-z0-9_.-]', '_', sample) + '_' + hashlib.sha256(sample.encode()).hexdigest()[:6]
    folder = Path(output_dir).resolve() / 'igv' / safe_sample / f'top{top_n}_{digest}'
    folder.mkdir(parents=True,exist_ok=True)
    shots = folder/'snapshots';shots.mkdir(exist_ok=True)
    for event in events:
        event['snapshot'] = f"{event['tool_key']}_{event['rank']:02d}.png"
    (folder/'events.json').write_text(json.dumps({'sample':sample,**signature,'events':events,'errors':errors},indent=2))
    with (folder/'events.tsv').open('w',newline='') as handle:
        fields=['tool','rank','fusion','chrom1','pos1','strand1','chrom2','pos2','strand2','locus','snapshot']
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(events)
    with (folder/'breakpoints.bed').open('w') as handle:
        handle.write('track name="Fusion breakpoints" description="Top fusion breakpoint positions" itemRgb="On"\n')
        for event in events:
            for n in (1,2):
                name=re.sub(r'\s+', '_', f"{event['tool_key']}_{event['rank']}_{event['fusion']}_partner{n}")
                handle.write(f"{event[f'chrom{n}']}\t{event[f'pos{n}']-1}\t{event[f'pos{n}']}\t{name}\t0\t{event[f'strand{n}']}\t{event[f'pos{n}']-1}\t{event[f'pos{n}']}\t255,127,0\n")
    # BED is 0-based, half-open; IGV goto loci and input breakpoints are 1-based.
    expected = bam or (Path(bam_dir)/f'{sample}.bam').resolve()
    load = 'load '+quoted_command(expected)
    if index:load += ' '+quoted_command('index='+str(index))
    lines=['new','genome '+quoted_command(genome),load,
           'load '+quoted_command(folder/'breakpoints.bed'),'snapshotDirectory '+quoted_command(shots),
           'maxPanelHeight 600','viewaspairs','colorBy PAIR_ORIENTATION']
    for event in events:lines += ['goto '+event['locus'],'snapshot '+event['snapshot']]
    lines += ['exit']
    batch=folder/'igv_batch.txt';batch.write_text('\n'.join(lines)+'\n')
    write_viewer(folder/'interactive.html',sample,genome,events)
    if command and not bam_error and events:
        # Explicit CLI option authorizes launching IGV; no shell interpretation.
        with (folder/'igv.log').open('w') as log:
            try:
                result=subprocess.run(shlex.split(command)+['-b',str(batch)],stdout=log,stderr=subprocess.STDOUT,timeout=900)
                if result.returncode:errors.append(f'IGV exited with code {result.returncode}; see igv.log.')
            except (OSError,subprocess.TimeoutExpired) as exc:errors.append(f'IGV screenshot run failed: {exc}')
    elif command and bam_error:
        errors.append('IGV screenshot run was skipped: '+bam_error)
    return dict(sample=sample,folder=folder,events=events,errors=errors,bam_error=bam_error,
                top_n=top_n,genome=genome,output_dir=Path(output_dir).resolve())


def report_html(result):
    if not result:return ''
    folder=result['folder'];base=folder.relative_to(result['output_dir']).as_posix()
    link=lambda name:quote(base+'/'+name,safe='/')
    esc=html.escape
    count=sum((folder/'snapshots'/e['snapshot']).is_file() for e in result['events'])
    messages=[]
    if result['bam_error']:messages.append(result['bam_error'])
    messages.extend(result['errors'])
    pending='' if count==len(result['events']) else '<p>IGV screenshots are pending. Run the exported batch file in IGV Desktop, then regenerate this report with the same options.</p>'
    table_rows=[];figures=[]
    for event in result['events']:
        table_rows.append(f"<tr><td>{esc(event['tool'])}</td><td>{event['rank']}</td><td>{esc(event['fusion'])}</td><td>{esc(event['locus'])}</td></tr>")
        path=folder/'snapshots'/event['snapshot']
        if path.is_file():
            data=path.read_bytes()
            if not data.startswith(b'\x89PNG\r\n\x1a\n'):continue
            image=base64.b64encode(data).decode()
            figures.append(f"<figure class='igv-snapshot'><figcaption>{esc(event['tool'])} #{event['rank']}: {esc(event['fusion'])}<br>{esc(event['locus'])}</figcaption><img src='data:image/png;base64,{image}' alt='IGV read alignments for {esc(event['fusion'],quote=True)}'></figure>")
    notices=''.join('<p>'+esc(m)+'</p>' for m in messages)
    return f'''<section class="igv-section">
<h3 class="section-title">IGV Visualisation</h3>
<p class="interactive-link"><a href="{link('interactive.html')}" target="_blank" rel="noopener">Open {esc(result['sample'])} interactive.html</a></p>
<p>Selected events: top {result['top_n']} per caller.</p>
<p>Genome: {esc(result['genome'])}. {len(result['events'])} breakpoint pairs selected; {count} IGV screenshots available.
The caller's original ranking is retained. For repeated rows of a fusion pair, the highest-ranked row supplies its coordinates.
FusionCatcher ranks pairs using aggregated support; coordinates use the row with the highest spanning-pair then unique-read support.</p>
{notices}{pending}
<div class="screen-only"><p><a href="{link('igv_batch.txt')}" download>IGV screenshot batch</a> · <a href="{link('breakpoints.bed')}" download>Breakpoint BED</a> ·
<a href="{link('events.tsv')}" download>Selected events TSV</a></p>
<details><summary>Show interactive IGV in this report</summary>
<iframe loading="lazy" title="Interactive IGV for {esc(result['sample'],quote=True)}" src="{link('interactive.html')}" style="width:100%;height:900px;border:1px solid #ddd"></iframe></details></div>
<table><thead><tr><th>Caller</th><th>Rank</th><th>Fusion</th><th>IGV loci (1-based)</th></tr></thead><tbody>{''.join(table_rows)}</tbody></table>
{''.join(figures)}</section>'''
