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
# With run_nontumor_filter=true in samples.tsv for that label:
#   -> for_loupe/spatial_{label}_compv1_mmv2_nontumor_for_loupe.csv
# Legacy merged reference runs use ref_tag: v2 (untagged loupe/transfer paths).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "${ROOT}"

CONFIG="${ROOT}/config/config.yaml"
SAMPLES_TSV="${ROOT}/config/samples.tsv"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J74_3_A1_GFP_R1_008um}"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_J75_3_D1_GFP_R1_008um}"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_7_D1_GFP_R1_008um}"
#SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_008um}"
SAMPLE="${SAMPLE:-Tuyu_Apr2026_846531_HD_BA_WM_C_DNB_SM_D90_5_A1_GFP_R1_016um}"

read -r REF_TAG DEFAULT_TARGET <<< "$(
  python3 - "${CONFIG}" "${SAMPLES_TSV}" "${SAMPLE}" <<'PY'
import sys
from pathlib import Path

import pandas as pd

cfg_path = Path(sys.argv[1])
samples_path = Path(sys.argv[2])
label = sys.argv[3]

ref_tag = "v2"
for line in cfg_path.read_text().splitlines():
    line = line.split("#", 1)[0].strip()
    if line.startswith("ref_tag:"):
        ref_tag = line.split(":", 1)[1].strip().strip('"').strip("'")

suffix = "" if ref_tag == "v2" else f"_{ref_tag}"

df = pd.read_csv(samples_path, sep="\t", comment="#")
df.columns = df.columns.str.strip()
matches = df.loc[df["label"].astype(str) == label]
if len(matches) == 1 and "run_nontumor_filter" in df.columns:
    flag = str(matches["run_nontumor_filter"].iloc[0]).lower()
    if flag in {"1", "true", "yes", "y"}:
        suffix = f"{suffix}_nontumor"

target = f"for_loupe/spatial_{label}{suffix}_for_loupe.csv"
print(ref_tag, target)
PY
)"

TARGET="${1:-${DEFAULT_TARGET}}"
FORCE="${FORCE:-none}"
#FORCE="${FORCE:-run_cellestimate_nontumor}"

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
