#!/bin/bash
# Production campaign controller (robust).
# Submits one batch (<=100 array tasks) at a time and waits for THAT job id to
# fully leave the queue before submitting the next, so concurrent jobs never
# exceed 100. Robust to transient squeue failures: a squeue error is treated as
# "keep waiting", never as "done" (the bug that overlapped batches earlier).
# Resumable: completed decks are skipped (guard in run_batch.sbatch); fully
# finished batches are not resubmitted.
set -u
cd "$(dirname "$(readlink -f "$0")")"
USER=$(whoami)
LOG=campaign.log

# count my batch jobs in queue; echo -1 on squeue failure (NEVER 0 on error)
count_mine() {
  local out rc
  out=$(squeue -u "$USER" -h -o "%j" 2>/dev/null); rc=$?
  [ $rc -ne 0 ] && { echo -1; return; }
  printf '%s\n' "$out" | grep -c '^mol_b'
}

# block until a specific job id has no tasks left (two consecutive empty reads)
wait_job_gone() {
  local jid=$1 empties=0 out rc
  while true; do
    out=$(squeue -h -j "$jid" -o "%i" 2>/dev/null); rc=$?
    if [ "$rc" -ne 0 ]; then sleep 60; continue; fi   # glitch: keep waiting
    if [ -z "$out" ]; then
      empties=$((empties+1)); [ "$empties" -ge 2 ] && return 0
    else
      empties=0
    fi
    sleep 30
  done
}

echo "[$(date)] campaign start (robust controller)" | tee -a "$LOG"
for d in $(ls -d batch*/ | sort); do
  d=${d%/}
  n_in=$(wc -l < "$d/inputs.list")
  n_done=$(grep -l "This run used" "$d"/*.out 2>/dev/null | wc -l)
  if [ "$n_done" -ge "$n_in" ]; then
    echo "[$(date)] $d complete ($n_done/$n_in) - skip" | tee -a "$LOG"; continue
  fi
  # ensure NONE of my batch jobs are active before submitting (robust zero-check)
  while true; do c=$(count_mine); [ "$c" = "0" ] && break; sleep 30; done
  jid=$(cd "$d" && sbatch --parsable run_batch.sbatch)
  echo "[$(date)] submitted $d -> $jid ($n_done/$n_in already done)" | tee -a "$LOG"
  # confirm it registered before waiting on it
  for i in $(seq 1 12); do
    [ -n "$(squeue -h -j "$jid" -o %i 2>/dev/null)" ] && break; sleep 5
  done
  wait_job_gone "$jid"
  echo "[$(date)] $d drained" | tee -a "$LOG"
done
echo "[$(date)] campaign COMPLETE: $(grep -l 'This run used' batch*/*.out 2>/dev/null | wc -l)/2000 decks" | tee -a "$LOG"
