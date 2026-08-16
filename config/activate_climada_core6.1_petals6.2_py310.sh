#!/bin/sh

CLIMADA_ENV_PREFIX=/mnt/sdb_test/tang/zengjun/TC_Road_Risk/software/conda-envs/climada-core6.1-petals6.2-py310
export CLIMADA_ENV_PREFIX
export PATH="$CLIMADA_ENV_PREFIX/bin:$PATH"
export PROJ_DATA="$CLIMADA_ENV_PREFIX/share/proj"
export GDAL_DATA="$CLIMADA_ENV_PREFIX/share/gdal"
export GDAL_DRIVER_PATH="$CLIMADA_ENV_PREFIX/lib/gdalplugins"
export CPL_ZIP_ENCODING=UTF-8
