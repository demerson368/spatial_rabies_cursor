#!/bin/bash
# Poll LSF job status for Snakemake 7 cluster-status.
# Must ALWAYS print exactly one of: running | success | failed
# and exit 0 — a non-zero exit makes Snakemake abort with
# "Failed to obtain job status".
#
# Snakemake appends the external job id as argv[1].
# Do NOT use set -e: bjobs often returns non-zero when a job is not yet
# visible or has aged out of the active queue.

JOBID="${1:-}"
JOBID="${JOBID#\{}"
JOBID="${JOBID%\}}"

if [[ -z "${JOBID}" || "${JOBID}" == "jobid" ]]; then
  echo running
  exit 0
fi

stat="$(bjobs -noheader -o 'stat' "${JOBID}" 2>/dev/null | head -1 | tr -d '[:space:]' || true)"

if [[ -n "${stat}" ]]; then
  case "${stat}" in
    RUN|PEND|PSUSP|USUSP|SSUSP|WAIT|PROV)
      echo running
      exit 0
      ;;
    DONE)
      echo success
      exit 0
      ;;
    EXIT|ZOMBI|UNKWN)
      echo failed
      exit 0
      ;;
    *)
      echo running
      exit 0
      ;;
  esac
fi

hist="$(bjobs -a -noheader -o 'stat' "${JOBID}" 2>/dev/null | head -1 | tr -d '[:space:]' || true)"

case "${hist}" in
  DONE)
    echo success
    ;;
  EXIT|ZOMBI|UNKWN)
    echo failed
    ;;
  RUN|PEND|PSUSP|USUSP|SSUSP|WAIT|PROV)
    echo running
    ;;
  "")
    # Job not visible yet (LSF lag) or aged out — keep polling
    echo running
    ;;
  *)
    echo running
    ;;
esac
exit 0
