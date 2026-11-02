#!/usr/bin/env bash
# Optional wrapper around the IGV Desktop launcher supplied by IGV.
# Example: export IGV_HOME=/opt/IGV_2.xx.x
# Then: bash run_igv.sh -b /absolute/path/to/igv_batch.txt
set -euo pipefail

: "${IGV_HOME:?Set IGV_HOME to the directory containing the installed IGV igv.sh}"
launcher="${IGV_HOME%/}/igv.sh"
if [[ ! -f "$launcher" ]]; then
    printf 'IGV launcher not found: %s\n' "$launcher" >&2
    exit 1
fi

# Preserve all arguments, including the -b batch-file argument added by the report.
# IGV's own launcher handles its Java runtime, classpath and memory settings.
exec bash "$launcher" "$@"
