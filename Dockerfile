FROM snakemake/snakemake:latest

# NOTE rsync only used with local scratch
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# TODO: only move workflow files
COPY . /V-pipe

VOLUME /work
WORKDIR /work

RUN echo '{"output": {"snv": true, "local": true, "global": true, "visualization": true, "QA": true}}' > config.yaml \
 && rm -f samples samples.tsv && ln -sf /V-pipe/tests/data/hiv/ ./samples  \
 && snakemake -s /V-pipe/vpipe.snake -j 1 --conda-create-envs-only --use-conda --conda-prefix /conda_prefix --config "general={virus_base_config: hiv}" \
 && rm -f samples samples.tsv && ln -sf /V-pipe/tests/data/sars-cov-2 samples && cp -f samples/samples.tsv ./samples.tsv \
 && snakemake -s /V-pipe/vpipe.snake -j 1 --conda-create-envs-only --use-conda --conda-prefix /conda_prefix --config "general={virus_base_config: sars-cov-2}" \
 && rm -f samples samples.tsv config.yaml \
 && mamba clean --yes --all && rm -rf /opt/conda/pkgs

ENTRYPOINT [ \
    "snakemake", \
    "-s", "/V-pipe/vpipe.snake", \
    "--use-conda", \
    "--conda-prefix", "/conda_prefix" \
]
