#!/bin/bash
# Submit a Snakemake job to LSF (Snakemake 7 pipes job script to stdin).
# Prints ONLY the LSF job ID to stdout so Snakemake tracks cluster jobs.

set -euo pipefail

RULE="${SMK_RULE:-snakemake}"
LABEL="${SMK_LABEL:-}"
MEM_MB="${SMK_MEM_MB:-64000}"
QUEUE="${SMK_QUEUE:-none}"
GPU="${SMK_GPU:-0}"

mkdir -p logs/cluster

JOBNAME="smk-${RULE}"
if [[ -n "${LABEL}" ]]; then
  JOBNAME="smk-${RULE}-${LABEL}"
fi
JOBNAME="${JOBNAME:0:120}"

BSUB=(
  bsub
  -J "${JOBNAME}"
  -oo "logs/cluster/${RULE}.out"
  -eo "logs/cluster/${RULE}.err"
  -R "rusage[mem=${MEM_MB}]"
)

if [[ -n "${QUEUE}" && "${QUEUE}" != "none" && "${QUEUE}" != "0" ]]; then
  BSUB+=(-q "${QUEUE}")
fi

if [[ "${GPU}" != "0" ]]; then
  BSUB+=(-gpu "num=1/host")
  BSUB+=(-env "PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True")
fi

output="$("${BSUB[@]}" 2>&1)" || {
  echo "${output}" >&2
  exit 1
}

jobid="$(echo "${output}" | sed -n 's/.*<\([0-9][0-9]*\)>.*/\1/p' | head -1)"
if [[ -z "${jobid}" ]]; then
  echo "Failed to parse LSF job id from: ${output}" >&2
  exit 1
fi

echo "${output}" >&2
echo "${jobid}"
