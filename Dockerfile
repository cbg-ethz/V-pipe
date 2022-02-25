# syntax=docker/dockerfile:1.3

# ====================================================================================
# NOTE this docker uses `--network=none` option for RUN directives, you need to either
#  - activate BuildKit with `export DOCKER_BUILDKIT=1`
#  - run this using buildx
# ====================================================================================

# NOTE the following must match the ENTRYPOINT
ARG install_path=/opt/V-dock
ARG vpipe_path=${install_path}/V-pipe
ARG envs_path=${install_path}/conda_envs
ARG test_data=/test-data
ARG virus_download_list
ARG snaketag=stable

###
### Stage 0: download conda environments
###
FROM snakemake/snakemake:${snaketag} AS create-envs

ARG install_path
ARG vpipe_path
ARG envs_path
ARG test_data

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
    jdupes

# TODO: only move workflow files
#COPY . ${vpipe_path}
WORKDIR ${vpipe_path}/
COPY LICENSE.md ./LICENSE.md
COPY workflow ./workflow
COPY resources ./resources
COPY config ./config
COPY utils ./utils
COPY init_project.sh ./init_project.sh

COPY tests/data ${test_data}

WORKDIR /work

# configuration: activate all steps
RUN mkdir config \
 && printf 'output:\n  snv: true\n  local: true\n  global: true\n  visualization: true\n  diversity: true\n  QA: true\n  upload: true\nupload:\n  orig_cram: true' > config/config.yaml

# TODO harmonize list with CI tests and Docker tests
RUN for virus in ${virus_download_list:-$(ls ${test_data}/)}; do printf '\n\n\e[36;1mvirus: %s\e[0m\n' "${virus}" \
 &&   ln -sf "${test_data}/${virus}/" ./samples \
 &&   if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./config/samples.tsv; fi \
 &&   PYTHONUNBUFFERED=1 snakemake -s ${vpipe_path}/workflow/Snakefile -j 1 --conda-create-envs-only --use-conda --conda-prefix ${envs_path} --config "general={virus_base_config: ${virus}}" \
 &&   rm -f samples config/samples.tsv \
  ; done \
 && jdupes -Lr ${envs_path}/


###
### Stage 1: base layer with V-pipe and environments
###
FROM snakemake/snakemake:${snaketag} AS vpipe-tests-base
ARG install_path

# NOTE rsync only used with local scratch
# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# NOTE V-pipe/envs/*.yaml and conda_prefix/* must be in sync so that env checksums match
COPY --from=create-envs ${install_path} ${install_path}



###
### Test 1: test the base layer with hiv
###
FROM vpipe-tests-base AS test_hiv
ARG install_path
ARG vpipe_path
ARG envs_path
ARG test_data
ENV virus=hiv

WORKDIR /work
RUN mkdir config \
 && printf 'output:\n  snv: true\n  local: true\n  global: false\n  visualization: true\n  diversity: true\n  QA: true\n  upload: true\nupload:\n  orig_cram: true' > config/config.yaml
COPY --from=create-envs ${test_data}/${virus} ./samples
RUN if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./config/samples.tsv; fi
# NOTE see top comment if `--network=none` breaks build process
RUN --network=none \
    PYTHONUNBUFFERED=1 snakemake -s ${vpipe_path}/workflow/Snakefile -j 4 --use-conda --conda-prefix ${envs_path} --config "general={virus_base_config: ${virus}}" \
 && echo "$(date --iso-8601=sec ; grep -E 'failed|for error' .snakemake/log/*.snakemake.log)" > ${install_path}/${virus}.teststamp



###
### Test 2: test the base layer with sars-cov-2
###
FROM vpipe-tests-base AS test_sars-cov-2
ARG install_path
ARG vpipe_path
ARG envs_path
ARG test_data
ENV virus=sars-cov-2

WORKDIR /work
RUN mkdir config \
 && printf 'output:\n  snv: true\n  local: true\n  global: false\n  visualization: true\n  diversity: true\n  QA: true\n  upload: true\nupload:\n  orig_cram: true' > config/config.yaml
COPY --from=create-envs ${test_data}/${virus} ./samples
RUN if test -e samples/samples.tsv; then cp -f samples/samples.tsv ./config/samples.tsv; fi
# NOTE see top comment if `--network=none` breaks build process
RUN --network=none \
    PYTHONUNBUFFERED=1 snakemake -s ${vpipe_path}/workflow/Snakefile -j 4 --use-conda --conda-prefix ${envs_path} --config "general={virus_base_config: ${virus}}" \
 && echo "$(date --iso-8601=sec ; grep -E 'failed|for error' .snakemake/log/*.snakemake.log)" > ${install_path}/${virus}.teststamp



###
### Final base: gather tests
###
FROM vpipe-tests-base as vpipe-final-base
ARG install_path

# NOTE individual test can be forced using the following on github actions
COPY --from=test_hiv	${install_path}/hiv.teststamp	${install_path}
COPY --from=test_sars-cov-2	${install_path}/sars-cov-2.teststamp	${install_path}



###
### Final stage: setup image ready to run
###

# =============================================
# NOTE this will *skip* tests on GitHub actions
#FROM vpipe-test-base
#ARG install_path
#ARG vpipe_path
#ARG envs_path
# ---------------------------------------------
# HACK this will *force* tests on GitHub actions
FROM snakemake/snakemake:${snaketag}
ARG install_path
ARG vpipe_path
ARG envs_path

# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY --from=vpipe-final-base ${install_path} ${install_path}
# =============================================

LABEL maintainer="V-pipe Dev Team <v-pipe@bsse.ethz.ch>"
VOLUME /work
WORKDIR /work

# NOTE current docker versions do not offer a way to bake the content of an ARG into an ENTRYPOINT
ENTRYPOINT [ \
    "snakemake", \
    "-s", "/opt/V-dock/V-pipe/workflow/Snakefile", \
    "--use-conda", \
    "--conda-prefix", "/opt/V-dock/conda_envs" \
]
