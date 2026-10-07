#!/usr/bin/env bash
# Checks that the options used in Table 3 of the paper exist in the current tool versions,
# and records the CPU features needed to explain the Parabricks SIGILL.
# Order of use for each tool: binary on PATH -> conda/mamba env in $HOME/tool_check_envs
# (Bioconda, works on aarch64) -> Apptainer BioContainers image. Nothing is installed system-wide.
set -u
LOG=check_tools.log
: > "$LOG"
ENVROOT=$HOME/tool_check_envs
export CONDA_PKGS_DIRS=$ENVROOT/pkgs
mkdir -p "$ENVROOT"

echo "### host"
uname -m
grep -m1 -E "^(model name|CPU part)" /proc/cpuinfo
lscpu | grep -E "^(Architecture|Model name|Vendor ID|CPU\(s\)|Flags)" | cut -c1-400
echo "SVE: $(grep -m1 -ow -E 'sve' /proc/cpuinfo || echo no)  SVE2: $(grep -m1 -ow -E 'sve2' /proc/cpuinfo || echo no)"
command -v nvidia-smi >/dev/null && nvidia-smi -L
free -g | head -2

RT=$(command -v apptainer || command -v singularity || true)
CONDA=$(command -v mamba || command -v micromamba || command -v conda ||
        ls /data/galaxy/*/_conda/bin/conda /data/galaxy/*/*/_conda/bin/conda 2>/dev/null | head -1 || true)
echo "apptainer: ${RT:-none}   conda: ${CONDA:-none}"
echo "### tools"

# run <bioconda package> <command...>
run() {
  local pkg=$1 cmd=$2; shift 2
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "[PATH $(command -v "$cmd")]" >> "$LOG"; "$cmd" "$@" 2>&1
  elif [ -n "$CONDA" ]; then
    local p=$ENVROOT/$pkg
    if [ ! -x "$p/bin/$cmd" ]; then
      "$CONDA" create -y -q -p "$p" -c conda-forge -c bioconda "$pkg" >> "$LOG" 2>&1
    fi
    echo "[conda $pkg $(ls "$p/conda-meta" 2>/dev/null | grep -m1 "^$pkg-[0-9]")]" >> "$LOG"
    "$p/bin/$cmd" "$@" 2>&1
  elif [ -n "$RT" ]; then
    local tag
    tag=$(curl -s "https://quay.io/api/v1/repository/biocontainers/$pkg/tag/?limit=1&onlyActiveTags=true" |
          python3 -c 'import json,sys; print(json.load(sys.stdin)["tags"][0]["name"])' 2>/dev/null)
    echo "[image $pkg:$tag]" >> "$LOG"
    "$RT" exec "docker://quay.io/biocontainers/$pkg:$tag" "$cmd" "$@" 2>&1
  else
    echo "NO_RUNTIME"
  fi
}

# check <label> <package> "<flag1>|<flag2>|..." <command...>
check() {
  local label=$1 pkg=$2 flags=$3; shift 3
  local out; out=$(run "$pkg" "$@")
  { echo "=== $label: $*"; echo "$out" | head -n 400; } >> "$LOG"
  local miss=""
  IFS='|' read -ra F <<< "$flags"
  for f in "${F[@]}"; do grep -qF -- "$f" <<< "$out" || miss="$miss $f"; done
  if [ -z "$miss" ]; then echo "OK      $label"; else echo "MISSING $label:$miss"; fi
}

check "GATK HaplotypeCaller"  gatk4 "--reference|--input|--output|-ERC" gatk HaplotypeCaller --help
check "GATK Mutect2"          gatk4 "--normal-sample|-normal|--input" gatk Mutect2 --help
check "GATK Funcotator"       gatk4 "--variant|--reference|--ref-version|--data-sources-path|--output-file-format" gatk Funcotator --help
check "SnpEff ann"            snpeff "ann" snpEff -h
check "PLINK 2"               plink2 "--vcf|--make-bed|--bfile|--pheno|--glm|--out" plink2 --help
check "VCFtools"              vcftools "--gzvcf|--minQ|--max-missing|--recode|--out" vcftools --help
check "bcftools view -i"      bcftools "--include" bcftools view
check "bcftools mpileup"      bcftools "--fasta-ref" bcftools mpileup
check "featureCounts"         subread "--countReadPairs|-p|-T|-a|-o" featureCounts
check "HISAT2"                hisat2 "-x|-1|-2|-S|--dta|-p" hisat2 --help
check "Picard MarkDuplicates" picard "INPUT|OUTPUT|METRICS_FILE" picard MarkDuplicates --help
check "Trimmomatic"           trimmomatic "-threads|PE" trimmomatic
check "Bismark extractor"     bismark "--bedGraph|--gzip" bismark_methylation_extractor --help
check "miRDeep2 quantifier"   mirdeep2 "-p|-m|-r|-t" quantifier.pl
check "clusterProfiler"       bioconductor-clusterprofiler "OrgDb|ont|organism" Rscript -e 'suppressMessages(library(clusterProfiler)); print(names(formals(enrichGO))); print(names(formals(enrichKEGG)))'
check "DESeq2"                bioconductor-deseq2 "countData|colData|design" Rscript -e 'suppressMessages(library(DESeq2)); print(names(formals(DESeqDataSetFromMatrix)))'
echo "full output: $(pwd)/$LOG   (remove the test environments with: rm -rf $ENVROOT)"
