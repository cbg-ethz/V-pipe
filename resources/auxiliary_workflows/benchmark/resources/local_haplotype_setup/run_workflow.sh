#!/usr/bin/env bash

bsub \
  -N \
  -R 'rusage[mem=5000]' \
  -W 120:00 \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  --rerun-incomplete \
  -pr \
  --cores 200 \
  --use-conda \
  --latency-wait 30 \
  --show-failed-logs \
  "$@"
