#!/bin/bash
# Prepare Lin inputs for MPI-ESM1-2-LR future windows. Sequential, nice'd,
# so it does not starve the 192-way hazard workers.
set -euo pipefail
ROOT=/mnt/sdb_test/tang/zengjun/TC_Road_Risk
PY=$ROOT/software/venvs/lin-v1.1-py310/bin/python
CDO=$ROOT/software/cdo-2.6.1-env/bin/cdo
ACCEPT=$ROOT/runs/cmip6_acceptance/MPI-ESM1-2-LR_r1i1p1f1.json
LOGDIR=$ROOT/runs/lin_environment/MPI-ESM1-2-LR/prepare_future_logs
mkdir -p "$LOGDIR"
cd "$ROOT"
for exp in ssp126 ssp245 ssp370 ssp585; do
  for pair in 2041:2060 2081:2100; do
    start=${pair%%:*}
    end=${pair##*:}
    out=$ROOT/data/cmip6/derived/lin/MPI-ESM1-2-LR/${exp}/r1i1p1f1/${start}-${end}
    if [[ -f $out/manifest.json ]]; then
      echo "skip existing $exp $start-$end"
      continue
    fi
    echo "prepare $exp $start-$end"
    nice -n 10 "$PY" "$ROOT/code/prepare_lin_inputs.py" \
      --acceptance-manifest "$ACCEPT" \
      --experiment "$exp" \
      --start-year "$start" \
      --end-year "$end" \
      --output "$out" \
      --cdo "$CDO" \
      > "$LOGDIR/${exp}_${start}-${end}.log" 2>&1
    echo "done $exp $start-$end"
  done
done
echo "all future windows prepared or skipped"
