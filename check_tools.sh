#!/bin/bash

# check_tools.sh — verify all pipeline tool dependencies before submitting jobs.
#
# For each required executable this script:
#   1. Checks whether the tool is already on PATH.
#   2. If not, attempts to load the module name configured in pipeline.conf.
#   3. If neither works, reports the tool as MISSING with installation guidance.
#
# JAR files, the adapter FASTA, and bowtie2 index prefixes are also validated.
#
# Usage (run interactively on a login / interactive node):
#   bash check_tools.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="${SCRIPT_DIR}/pipeline.conf"

if [[ ! -f "$CONF" ]]; then
    echo "ERROR: pipeline.conf not found at $CONF"
    exit 1
fi
source "$CONF"

# Resolve tilde in JAR/adapter paths
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR/#\~/$HOME}"
TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS/#\~/$HOME}"
PICARD_JAR="${PICARD_JAR/#\~/$HOME}"

# ── Ensure the 'module' command is available ──────────────────────────────────
# It is a shell function sourced from the cluster init scripts; it is not
# present in non-login shells.  Try common init locations if needed.
if ! command -v module &>/dev/null; then
    for _init in \
        /etc/profile.d/modules.sh \
        /usr/share/Modules/init/bash \
        /usr/local/Modules/init/bash \
        "${MODULESHOME:-/dev/null}/init/bash" \
        /opt/modules/init/bash \
        /opt/lmod/lmod/init/bash; do
        if [[ -f "$_init" ]]; then
            # shellcheck source=/dev/null
            source "$_init"
            break
        fi
    done
fi

MODULES_AVAILABLE=false
command -v module &>/dev/null && MODULES_AVAILABLE=true

# ── Counters ──────────────────────────────────────────────────────────────────
_PASS=0
_FAIL=0
_LOADED_MODULES=()   # track modules loaded during this check

# ── Helpers ───────────────────────────────────────────────────────────────────

# check_exe NAME CMD MODULE INSTALL_HINT
#   Verifies that CMD is executable.  Tries PATH first, then module load.
check_exe() {
    local name=$1 cmd=$2 module=$3 hint=$4
    printf "  %-14s : " "$name"

    if command -v "$cmd" &>/dev/null; then
        echo "OK        $(command -v "$cmd")"
        (( _PASS++ ))
        return 0
    fi

    # Not on PATH — try the configured module
    if [[ -n "$module" ]]; then
        if [[ "$MODULES_AVAILABLE" == false ]]; then
            echo "MISSING   (not on PATH; 'module' command unavailable — cannot try '$module')"
        else
            # Run module load in the CURRENT shell (not a subshell) so that PATH
            # changes take effect here.  Redirect stderr to a temp file to avoid
            # polluting output; module systems like Lmod print status messages
            # ("Loading X") to stderr even on success, so exit code is the only
            # reliable indicator of whether the load itself failed.
            local tmpfile mod_exit
            tmpfile=$(mktemp)
            module load "$module" 2>"$tmpfile"
            mod_exit=$?
            rm -f "$tmpfile"

            if command -v "$cmd" &>/dev/null; then
                echo "OK        $(command -v "$cmd")  [loaded module: $module]"
                _LOADED_MODULES+=("$module")
                (( _PASS++ ))
                return 0
            elif [[ $mod_exit -ne 0 ]]; then
                echo "MISSING   (not on PATH; 'module load $module' failed — check module name in pipeline.conf)"
            else
                echo "MISSING   (not on PATH; module '$module' loaded but '$cmd' not found — wrong module name?)"
            fi
        fi
    else
        echo "MISSING   (not on PATH; no module configured in pipeline.conf)"
    fi

    echo "               --> $hint"
    (( _FAIL++ ))
    return 1
}

# check_file LABEL PATH
#   Verifies that a file exists and is readable.
check_file() {
    local label=$1 path=$2
    printf "  %-22s : " "$label"
    if [[ -f "$path" && -r "$path" ]]; then
        echo "OK        $path"
        (( _PASS++ ))
    else
        echo "MISSING   $path"
        (( _FAIL++ ))
    fi
}

# check_index LABEL PREFIX
#   Verifies that at least one bowtie2 index file (.1.bt2 or .1.bt2l) exists.
check_index() {
    local label=$1 prefix=$2
    printf "  %-22s : " "$label"
    if [[ -f "${prefix}.1.bt2" || -f "${prefix}.1.bt2l" ]]; then
        echo "OK        $prefix"
        (( _PASS++ ))
    else
        echo "MISSING   $prefix"
        echo "               --> Build with: bowtie2-build <genome.fa> $prefix"
        (( _FAIL++ ))
    fi
}

# ── Run checks ────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo " Pipeline tool check"
echo " Config : $CONF"
echo "============================================================"

echo ""
echo "[ Executables ]"
check_exe "java"       java       ""                  \
    "Install JDK 8+: https://adoptium.net  (or load a java module)"
check_exe "bowtie2"    bowtie2    "$MODULE_BOWTIE2"   \
    "Install bowtie2: https://bowtie-bio.sourceforge.net/bowtie2"
check_exe "samtools"   samtools   "$MODULE_SAMTOOLS"  \
    "Install samtools: https://www.htslib.org"
check_exe "bamCoverage" bamCoverage "$MODULE_DEEPTOOLS" \
    "Install deeptools: pip install deeptools  or  conda install -c bioconda deeptools"

echo ""
echo "[ JAR files ]"
check_file "TRIMMOMATIC_JAR"  "$TRIMMOMATIC_JAR"
check_file "PICARD_JAR"       "$PICARD_JAR"

echo ""
echo "[ Adapter file ]"
check_file "TRIMMOMATIC_ADAPTERS" "$TRIMMOMATIC_ADAPTERS"

echo ""
echo "[ Bowtie2 indices ]"
check_index "BOWTIE2_INDEX_MAIN"  "$BOWTIE2_INDEX_MAIN"
check_index "BOWTIE2_INDEX_SPIKE" "$BOWTIE2_INDEX_SPIKE"

# ── Summary ───────────────────────────────────────────────────────────────────

echo ""
echo "============================================================"
echo " Results: ${_PASS} passed, ${_FAIL} failed"
if [[ ${#_LOADED_MODULES[@]} -gt 0 ]]; then
    echo " Note: the following modules were loaded during this check:"
    for m in "${_LOADED_MODULES[@]}"; do echo "   module load $m"; done
    echo " Add these to your ~/.bashrc or jobscript if not already present."
fi
echo "============================================================"
echo ""

if [[ $_FAIL -gt 0 ]]; then
    echo "Action required: install or configure the MISSING items above,"
    echo "then edit pipeline.conf and re-run this script."
    exit 1
else
    echo "All checks passed. You are ready to run the pipeline."
    exit 0
fi
