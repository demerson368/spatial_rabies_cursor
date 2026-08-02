#!/bin/bash
# Run Snakemake workflow for one sample (LSF cluster submission).
#
# Usage:
#   sh snakemake.sh                          # default SAMPLE + ref_tag from config
#   sh snakemake.sh <loupe_target.csv>       # explicit target path
#   SAMPLE=J74_3_A1_GFP_R1_008um sh snakemake.sh
#   FORCE=label_transfer sh snakemake.sh     # force-rerun rule(s) + downstream
#   FORCE=all sh snakemake.sh                # force all rules (-F)
#
# Reference variants: set ref_out/ref_file/ref_rds/ref_label and ref_tag in
# config/config.yaml. Outputs are namespaced by ref_tag, e.g.
#   ref_tag: compv1_mmv2
#   -> for_loupe/spatial_{label}_compv1_mmv2_for_loupe.csv
# Legacy merged reference runs use ref_tag: v2 (untagged loupe/transfer paths).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "${ROOT}"

CONFIG="${ROOT}/config/config.yaml"
read -r REF_TAG <<< "$(
  python3 - "${CONFIG}" <<'PY'
import sys
from pathlib import Path

cfg_path = Path(sys.argv[1])
ref_tag = "v2"
for line in cfg_path.read_text().splitlines():
    line = line.split("#", 1)[0].strip()
    if line.startswith("ref_tag:"):
        ref_tag = line.split(":", 1)[1].strip().strip('"').strip("'")
print(ref_tag)
PY
)"

#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um}"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um}"
SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um}"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um}"
SUFFIX=""
if [[ "${REF_TAG}" != "v2" ]]; then
  SUFFIX="_${REF_TAG}"
fi
DEFAULT_TARGET="for_loupe/spatial_${SAMPLE}${SUFFIX}_for_loupe.csv"
TARGET="${1:-${DEFAULT_TARGET}}"
FORCE="${FORCE:-none}"

chmod +x profiles/lsf/cluster-submit.sh profiles/lsf/status_lsf.sh profiles/lsf/submit_lsf.sh 2>/dev/null || true

snakemake -s workflow/Snakefile --configfile config/config.yaml --unlock >/dev/null 2>&1 || true
# Optional: CLEAN=1 sh snakemake.sh clears Snakemake metadata cache (can trigger extra reruns).
if [[ "${CLEAN:-0}" == "1" ]]; then
  rm -rf .snakemake/metadata 2>/dev/null || true
fi

FORCE_ARGS=()
case "${FORCE}" in
  all|"-F")
    FORCE_ARGS+=(-F)
    ;;
  none|"")
    ;;
  *)
    # e.g. FORCE=label_transfer or FORCE="label_transfer assemble_loupe"
    # shellcheck disable=SC2206
    FORCE_ARGS+=(--forcerun ${FORCE})
    ;;
esac




echo "REF_TAG=${REF_TAG}"
echo "TARGET=${TARGET}"
echo "FORCE=${FORCE}"

snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profiles/lsf \
  "${FORCE_ARGS[@]}" \
  --cluster "${ROOT}/profiles/lsf/cluster-submit.sh {dependencies}" \
  --cluster-status "${ROOT}/profiles/lsf/status_lsf.sh" \
  --cores 10 \
  --jobs 50 \
  --latency-wait 120 \
  "${TARGET}"
