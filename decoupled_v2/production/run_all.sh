#!/bin/bash
# Controller: submit the decoupled_v2 cross-section campaign in chunks that
# respect MaxSubmit=500, each array throttled to <=100 concurrent, draining one
# chunk before submitting the next. Resume-safe (sbatch skips finished decks).
set -u
cd "$(dirname "$0")"
CHUNK=180             # <= cluster submit cap of 200
THROTTLE=180          # concurrent tasks the user permits
TOTAL=$(wc -l < inputs.list)
echo "[$(date '+%F %T')] controller start: $TOTAL decks, chunk=$CHUNK, throttle=%$THROTTLE"

my_jobs() {               # total running+pending jobs for me; '' (empty) on error
  squeue -u pzuch007 -h -t PENDING,RUNNING,CONFIGURING,COMPLETING 2>/dev/null | grep -c .
}

wait_gone() {             # block until I have 0 jobs for 2 consecutive reads
  local empties=0 n
  while :; do
    n=$(my_jobs)
    if ! [[ "$n" =~ ^[0-9]+$ ]]; then empties=0; sleep 20; continue; fi  # glitch
    if [ "$n" -eq 0 ]; then
      empties=$((empties+1)); [ $empties -ge 2 ] && return 0
    else empties=0; fi
    sleep 20
  done
}

start=1
while [ "$start" -le "$TOTAL" ]; do
  end=$(( start + CHUNK - 1 )); [ "$end" -gt "$TOTAL" ] && end=$TOTAL
  jid=$(sbatch --parsable --array="${start}-${end}%${THROTTLE}" run_campaign.sbatch) || {
    echo "[$(date '+%F %T')] sbatch failed for ${start}-${end}; retry in 60s"; sleep 60; continue; }
  echo "[$(date '+%F %T')] submitted chunk ${start}-${end} as job $jid"
  wait_gone "$jid"
  echo "[$(date '+%F %T')] chunk ${start}-${end} drained"
  start=$(( end + 1 ))
done

done_n=$(grep -l "This run used" decks/*.out 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] ALL CHUNKS DONE. finished outputs: $done_n / $TOTAL"
