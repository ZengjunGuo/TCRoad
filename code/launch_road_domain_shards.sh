#!/bin/bash
set -euo pipefail
ROOT=/mnt/sdb_test/tang/zengjun/TC_Road_Risk
OUT=$ROOT/runs/hazard_production/lin_road_domain_300km_v1
source "$ROOT/software/activate-climada-core6.1-petals6.2-py310.sh"
cd "$ROOT"
export MPLCONFIGDIR
MPLCONFIGDIR=$(mktemp -d /tmp/tcroad-mpl.XXXXXX)
IDX="$1"
SHARDS="$2"
python code/run_lin_10k_batch.py \
  --sample "$ROOT/data/lin/samples/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/road_domain_300km/events_road_domain_300km.nc" \
  --fixed-r0-catalogue "$ROOT/data/lin/samples/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/road_domain_300km/fixed_r0_catalogue_v1.nc" \
  --fixed-r0-manifest "$ROOT/data/lin/samples/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/road_domain_300km/fixed_r0_catalogue_v1.manifest.json" \
  --event-worker "$ROOT/code/run_lin_event_worker.py" \
  --output-root "$OUT" \
  --scratch-root "$OUT/scratch" \
  --track-native-last-index 360 \
  --shard-index "$IDX" \
  --shard-count "$SHARDS" \
  --allow-non-10k-sample \
  --worker-arg=--tracks \
  --worker-arg="$ROOT/data/lin/tracks/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/GLx5000peryear_stream0/output/tracks_GL_MPI-ESM1-2-LR_historical_r1i1p1f1_199501_201412.nc" \
  --worker-arg=--cmip6-ta \
  --worker-arg="$ROOT/data/cmip6/derived/lin/MPI-ESM1-2-LR/historical/r1i1p1f1/1995-2014/input/MPI-ESM1-2-LR_historical_r1i1p1f1_ta_1995-2014.nc" \
  --worker-arg=--roads-nc \
  --worker-arg="$ROOT/data/osm/derived/planet-260803/motor_road_length_density_0p1deg.nc" \
  --worker-arg=--elevation-tif \
  --worker-arg="$ROOT/data/static/tcr_public/climada/topography_land_360as/v1/topography_land_360as.tif" \
  --worker-arg=--c-drag-tif \
  --worker-arg="$ROOT/data/static/tcr_public/climada/c_drag_500/v1/c_drag_500.tif"
