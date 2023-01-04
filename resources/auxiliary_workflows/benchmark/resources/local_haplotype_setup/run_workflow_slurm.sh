#!/usr/bin/env bash

sbatch \
  --mail-type=END \
  --mem-per-cpu=5000 \
  --time=120:00:00 \
  -o snake.out -e snake.err \
snakemake \
  --profile slurm \
  --rerun-incomplete \
  -pr \
  --cores 200 \
  --use-conda \
  --latency-wait 30 \
  --show-failed-logs \
  "$@"
