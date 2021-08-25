# syntax=docker/dockerfile:1

###
### Stage 0: download conda environments
###
FROM snakemake/snakemake:latest AS create-envs

RUN apt-get update && apt-get install -y --no-install-recommends \
    jdupes

# TODO: only move workflow files
#COPY . /V-dock/V-pipe
WORKDIR /V-dock/V-pipe/
COPY README.md 	./README.md
COPY LICENSE 	./LICENSE
COPY vpipe.snake 	./vpipe.snake
COPY rules 	./rules
COPY envs 	./envs
COPY scripts 	./scripts
COPY resources 	./resources
COPY config 	./config
COPY functions.sh 	./functions.sh
COPY utils 	./utils
COPY init_project.sh 	./init_project.sh

COPY tests/data /test-data

WORKDIR /work

# configuration: activate all steps
RUN echo 'output:\n  snv: true\n  local: true\n  global: true\n  visualization: true\n  QA: true' > config.yaml

# TODO harmonize list with CI tests and Docker tests
RUN for virus in $(ls /test-data/); do echo "\n\n\e[36;1mvirus: ${virus}\e[0m\n" \
 &&   ln -sf "/test-data/${virus}/" ./samples \
 &&   if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./samples.tsv; fi \
 &&   snakemake -s /V-dock/V-pipe/vpipe.snake -j 1 --conda-create-envs-only --use-conda --conda-prefix /V-dock/conda_prefix --config "general={virus_base_config: ${virus}}" \
 &&   rm -f samples samples.tsv \
  ; done

RUN jdupes -Lr /V-dock/conda_prefix/


###
### Stage 1: base layer with V-pipe and environments
###
FROM snakemake/snakemake:latest AS vpipe-base

# NOTE rsync only used with local scratch
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# NOTE V-pipe/envs/*.yaml and conda_prefix/* must be in sync so that env checksums match
COPY --from=create-envs /V-dock /V-dock



###
### Test 1: test the base layer with hiv
###
FROM vpipe-base AS test_hiv
ENV virus=hiv
WORKDIR /work
RUN echo 'output:\n  snv: true\n  local: true\n  global: false\n  visualization: true\n  QA: true' > config.yaml
COPY --from=create-envs /test-data/${virus} ./samples
RUN if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./samples.tsv; fi
RUN snakemake -s /V-dock/V-pipe/vpipe.snake -j 4 --use-conda --conda-prefix /V-dock/conda_prefix --config "general={virus_base_config: ${virus}}"



###
### Test 2: test the base layer with sars-cov-2
###
FROM vpipe-base AS test_sars-cov-2
ENV virus=sars-cov-2
WORKDIR /work
RUN echo 'output:\n  snv: true\n  local: true\n  global: false\n  visualization: true\n  QA: true' > config.yaml
COPY --from=create-envs /test-data/${virus} ./samples
RUN if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./samples.tsv; fi
RUN snakemake -s /V-dock/V-pipe/vpipe.snake -j 4 --use-conda --conda-prefix /V-dock/conda_prefix --config "general={virus_base_config: ${virus}}"



###
### Final stage: setup image ready to run
###
FROM vpipe-base

VOLUME /work
WORKDIR /work

ENTRYPOINT [ \
    "snakemake", \
    "-s", "/V-dock/V-pipe/vpipe.snake", \
    "--use-conda", \
    "--conda-prefix", "/V-dock/conda_prefix" \
]
