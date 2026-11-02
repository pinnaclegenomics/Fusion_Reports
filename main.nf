nextflow.enable.dsl=2

// Quote values as single shell arguments, including spaces and apostrophes.
def quote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def stagedOne(value) { value instanceof List ? value[0] : value }

def optionalFile(value) {
    value ? [file(value, checkIfExists: true)] : []
}

process GENERATE_FUSION_REPORTS {
    tag 'all samples'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // Include nested file contents in caching, including BAM/index changes.
    cache 'deep'

    input:
    path input_dir, stageAs: 'inputs'
    path bam_dirs, stageAs: 'bam'
    path metadata, stageAs: 'metadata/*'
    path logos, stageAs: 'logo/*'
    path references, stageAs: 'references/*'
    path background, stageAs: 'background/*'
    path scripts, stageAs: 'scripts'
    val options

    output:
    path 'reports', emit: reports

    script:
    def args = ['--inputs', input_dir, '--outdir', 'reports',
                '--top-n', options.top_n, '--igv-top-n', options.igv_top_n,
                '--igv-flank', options.igv_flank]
    if (bam_dirs) args += ['--bam-dir', stagedOne(bam_dirs)]
    if (metadata) args += ['--config', stagedOne(metadata)]
    if (logos) args += ['--logo', stagedOne(logos)]
    if (references) args += ['--references', stagedOne(references)]
    if (background) args += ['--tool-background', stagedOne(background)]
    if (options.bam_pattern) args += ['--bam-pattern', options.bam_pattern]
    if (options.igv_genome) args += ['--igv-genome', options.igv_genome]
    if (options.igv_command) args += ['--igv-command', options.igv_command]
    def commandArgs = args.collect { quote(it) }.join(' ')
    """
    export PYTHONDONTWRITEBYTECODE=1
    export MPLBACKEND=Agg
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    export XDG_CACHE_HOME="\$PWD/.cache"
    ${quote(options.python)} scripts/Generate_Fusion_reports_html_PDF.py ${commandArgs}
    ${quote(options.python)} scripts/validate_reports.py --inputs inputs --reports reports
    """
}

workflow {
    if (params.help) {
        log.info '''
CAB fusion reporting pipeline
  nextflow run main.nf -profile conda -params-file params.example.yaml

Required: --inputs DIRECTORY
Optional: --outdir DIRECTORY --top_n 5|10 --igv_top_n 5|10
          --report_config FILE --logo FILE --references FILE
          --bam_dir DIRECTORY --bam_pattern '{sample}/*.sorted.bam'
          --igv_genome hg38 --igv_flank 250 --igv_command 'IGV launcher'
          --tool_background FILE --python python3
See README.md for local, Conda, Slurm, and screenshot instructions.
'''
    } else {
        if (!params.inputs) error 'Provide --inputs with your fusion caller results directory.'
        if (!(params.top_n.toString() in ['5', '10'])) error '--top_n must be 5 or 10.'
        if (!(params.igv_top_n.toString() in ['5', '10'])) error '--igv_top_n must be 5 or 10.'
        if (!(params.igv_flank.toString() ==~ /[1-9][0-9]*/)) error '--igv_flank must be a positive integer.'
        def inputs = file(params.inputs, checkIfExists: true)
        if (!inputs.isDirectory()) error '--inputs must be a directory.'
        def bams = optionalFile(params.bam_dir)
        if (bams && !bams[0].isDirectory()) error '--bam_dir must be a directory.'
        if (params.igv_command && !params.bam_dir) error '--igv_command requires --bam_dir.'
        GENERATE_FUSION_REPORTS(
            inputs, bams, optionalFile(params.report_config), optionalFile(params.logo),
            optionalFile(params.references), optionalFile(params.tool_background),
            file("${projectDir}/scripts", checkIfExists: true),
            [top_n: params.top_n, igv_top_n: params.igv_top_n, igv_flank: params.igv_flank,
             bam_pattern: params.bam_pattern, igv_genome: params.igv_genome,
             igv_command: params.igv_command, python: params.python]
        )
    }
}
