FROM snakemake/snakemake:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
 && rm -rf /var/lib/apt/lists/*

# TODO: only move workflow files
COPY . /V-pipe

VOLUME /work
WORKDIR /work

RUN echo '{"output": {"snv": true, "local": true, "global": true, "visualization": true, "QA": true}}' > config.yaml \
 && ln -sf /V-pipe/tests/data/hiv samples \
 && snakemake -s /V-pipe/vpipe.snake -j 1 --conda-create-envs-only --use-conda --conda-prefix /conda_prefix --config "general={virus_base_config: hiv}" \
 && ln -sf /V-pipe/tests/data/sars-cov-2 samples \
 && snakemake -s /V-pipe/vpipe.snake -j 1 --conda-create-envs-only --use-conda --conda-prefix /conda_prefix --config "general={virus_base_config: sars-cov-2}" \
 && rm config.yaml

ENTRYPOINT [ \
    "snakemake", \
    "-s", "/V-pipe/vpipe.snake", \
    "--use-conda", \
    "--conda-prefix", "/conda_prefix" \
]
