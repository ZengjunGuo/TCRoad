#!/bin/bash
# Advance future Lin windows: prepared inputs -> environment -> GLx5000 tracks.
# Sequential per window. Never starts a second track job. Does not touch the
# 192-way historical hazard workers. prepare_future_lin_windows.sh may still
# be filling later prepared-input windows in parallel.
set -euo pipefail

ROOT=/mnt/sdb_test/tang/zengjun/TC_Road_Risk
PY=$ROOT/software/venvs/lin-v1.1-py310/bin/python
RUNTIME=$ROOT/scratch/lin_runtime/v1.1-global-periodic-seed20260809-stream0
LOGDIR=$ROOT/runs/lin_environment/MPI-ESM1-2-LR/future_pipeline_logs
YEAR_WORKERS=${YEAR_WORKERS:-8}
TRACK_TIMEOUT=${TRACK_TIMEOUT:-129600}
mkdir -p "$LOGDIR"
cd "$ROOT"

windows=(
  ssp126:2041:2060
  ssp126:2081:2100
  ssp245:2041:2060
  ssp245:2081:2100
  ssp370:2041:2060
  ssp370:2081:2100
  ssp585:2041:2060
  ssp585:2081:2100
)

log() {
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee -a "$LOGDIR/pipeline.log"
}

prepared_root() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/data/cmip6/derived/lin/MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}"
}

env_record() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/runs/lin_environment/MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}"
}

env_output() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/data/lin/environment/MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}"
}

env_work() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/scratch/lin_environment_work/MPI-ESM1-2-LR_${exp}_r1i1p1f1_${start}-${end}"
}

track_publish() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/data/lin/tracks/MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}/GLx5000peryear_stream0"
}

track_work() {
  local exp=$1 start=$2 end=$3
  echo "$ROOT/scratch/lin_track_window_work/MPI-ESM1-2-LR_${exp}_r1i1p1f1_${start}-${end}_GLx5000peryear_stream0"
}

manifest_pass() {
  local path=$1
  [[ -f "$path" ]] || return 1
  python3 - "$path" << 'PY'
import json, sys
data = json.load(open(sys.argv[1], encoding="utf-8"))
sys.exit(0 if data.get("status") in {"pass", "completed"} else 1)
PY
}

run_env() {
  local exp=$1 start=$2 end=$3
  local rec out work prep
  rec=$(env_record "$exp" "$start" "$end")
  out=$(env_output "$exp" "$start" "$end")
  work=$(env_work "$exp" "$start" "$end")
  prep=$(prepared_root "$exp" "$start" "$end")
  if manifest_pass "$rec/manifest.json"; then
    log "skip env $exp $start-$end (published)"
    return 0
  fi
  if [[ -e "$work" || -e "$rec" || -e "$out" ]]; then
    log "env $exp $start-$end has leftover paths; waiting rather than relaunching"
    return 1
  fi
  log "start env $exp $start-$end"
  nice -n 10 env PYTHONUNBUFFERED=1 "$PY" "$ROOT/code/run_lin_environment.py" \
    --prepared-root "$prep" \
    --formal-runtime "$RUNTIME" \
    --work-root "$work" \
    --output-root "$out" \
    --record-root "$rec" \
    --python "$PY" \
    > "$LOGDIR/${exp}_${start}-${end}_env.log" 2>&1
  log "done env $exp $start-$end"
}

run_tracks() {
  local exp=$1 start=$2 end=$3
  local pub work prep rec
  pub=$(track_publish "$exp" "$start" "$end")
  work=$(track_work "$exp" "$start" "$end")
  prep=$(prepared_root "$exp" "$start" "$end")
  rec=$(env_record "$exp" "$start" "$end")
  if manifest_pass "$pub/record/manifest.json"; then
    log "skip tracks $exp $start-$end (published)"
    return 0
  fi
  if [[ -e "$work" || -e "$pub" ]]; then
    log "tracks $exp $start-$end has leftover paths; not relaunching"
    return 1
  fi
  log "preflight tracks $exp $start-$end"
  env PYTHONUNBUFFERED=1 "$PY" "$ROOT/code/run_lin_tracks_climate_window.py" \
    --project "$ROOT" \
    --experiment "$exp" \
    --start-year "$start" \
    --end-year "$end" \
    --prepared-root "$prep" \
    --environment-manifest "$rec/manifest.json" \
    --formal-runtime "$RUNTIME" \
    --python "$PY" \
    --year-workers "$YEAR_WORKERS" \
    --track-timeout-seconds "$TRACK_TIMEOUT" \
    --preflight-only \
    > "$LOGDIR/${exp}_${start}-${end}_tracks_preflight.log" 2>&1
  log "start tracks $exp $start-$end workers=$YEAR_WORKERS"
  nice -n 5 env PYTHONUNBUFFERED=1 "$PY" "$ROOT/code/run_lin_tracks_climate_window.py" \
    --project "$ROOT" \
    --experiment "$exp" \
    --start-year "$start" \
    --end-year "$end" \
    --prepared-root "$prep" \
    --environment-manifest "$rec/manifest.json" \
    --formal-runtime "$RUNTIME" \
    --python "$PY" \
    --year-workers "$YEAR_WORKERS" \
    --track-timeout-seconds "$TRACK_TIMEOUT" \
    > "$LOGDIR/${exp}_${start}-${end}_tracks.log" 2>&1
  log "done tracks $exp $start-$end"
}

SLEEP_SECONDS=${SLEEP_SECONDS:-120}

while true; do
  all_done=1
  progressed=0
  for spec in "${windows[@]}"; do
    exp=${spec%%:*}
    rest=${spec#*:}
    start=${rest%%:*}
    end=${rest##*:}
    prep=$(prepared_root "$exp" "$start" "$end")
    rec=$(env_record "$exp" "$start" "$end")
    pub=$(track_publish "$exp" "$start" "$end")
    ework=$(env_work "$exp" "$start" "$end")
    twork=$(track_work "$exp" "$start" "$end")

    if manifest_pass "$pub/record/manifest.json"; then
      continue
    fi
    all_done=0

    if ! manifest_pass "$prep/manifest.json"; then
      log "waiting on prepared inputs $exp $start-$end"
      continue
    fi

    if ! manifest_pass "$rec/manifest.json"; then
      if pgrep -f "run_lin_environment.py .*MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}" >/dev/null; then
        log "env $exp $start-$end already running"
        continue
      fi
      if [[ -e "$ework" || -e "$rec" || -e "$(env_output "$exp" "$start" "$end")" ]]; then
        log "env $exp $start-$end leftover without a live process; inspect before relaunch"
        continue
      fi
      run_env "$exp" "$start" "$end"
      progressed=1
    fi

    if ! manifest_pass "$rec/manifest.json"; then
      continue
    fi

    if pgrep -f "run_lin_tracks_climate_window.py .*--experiment ${exp} .*--start-year ${start} .*--end-year ${end}" >/dev/null; then
      log "tracks $exp $start-$end already running"
      continue
    fi
    if [[ -e "$twork" || -e "$pub" ]]; then
      log "tracks $exp $start-$end leftover without a live process; inspect before relaunch"
      continue
    fi
    run_tracks "$exp" "$start" "$end"
    progressed=1
  done

  if [[ $all_done -eq 1 ]]; then
    log "all eight future windows published"
    exit 0
  fi
  if [[ $progressed -eq 0 ]]; then
    log "idle; sleeping ${SLEEP_SECONDS}s for prepare/env to advance"
    sleep "$SLEEP_SECONDS"
  fi
done
