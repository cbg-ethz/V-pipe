#!/usr/bin/env bash

source /miniconda3/bin/activate base

if [ "$1" = "vpipe" ]; then
    shift
    exec -a "$0" snakemake -s "/V-pipe/vpipe.snake" --use-conda "$@"
fi

exec "$@"
