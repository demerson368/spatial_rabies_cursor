#!/bin/bash
# Snakemake 7 cluster submit wrapper for LSF.
# Snakemake appends the compiled job script path as the LAST positional argument.
# Upstream LSF job IDs arrive via {dependencies}; prints ONLY the LSF job ID to stdout.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"
mkdir -p logs/cluster

JOBSCRIPT="${@: -1}"
DEP_IDS=()
if [[ $# -gt 1 ]]; then
  for arg in "${@:1:$#-1}"; do
    if [[ "${arg}" =~ ^[0-9]+$ ]]; then
      DEP_IDS+=("${arg}")
    fi
  done
fi

if [[ ! -f "${JOBSCRIPT}" ]]; then
  echo "Expected Snakemake job script as last argument, got: ${JOBSCRIPT}" >&2
  exit 1
fi

read -r RULE LABEL MEM_MB QUEUE GPU <<< "$(
  python3 - "${JOBSCRIPT}" <<'PY'
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
props = read_job_properties(jobscript)
resources = props.get("resources", {})
wildcards = props.get("wildcards", {})
print(
    props["rule"],
    wildcards.get("label", ""),
    resources.get("mem_mb", 64000),
    resources.get("queue", "none"),
    resources.get("gpu", 0),
)
PY
)"

echo "$(date -Iseconds) rule=${RULE} label=${LABEL} mem=${MEM_MB} queue=${QUEUE} gpu=${GPU} script=${JOBSCRIPT}" \
  >> logs/cluster/submit.log

JOBNAME="smk-${RULE}"
[[ -n "${LABEL}" ]] && JOBNAME="smk-${RULE}-${LABEL}"
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
  BSUB+=(-env "all,PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True")
fi

if [[ ${#DEP_IDS[@]} -gt 0 ]]; then
  DEP_EXPR=""
  for dep in "${DEP_IDS[@]}"; do
    DEP_EXPR+="done(${dep}) && "
  done
  DEP_EXPR="${DEP_EXPR% && }"
  BSUB+=(-w "${DEP_EXPR}")
fi

BSUB+=(/bin/bash "${JOBSCRIPT}")

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
