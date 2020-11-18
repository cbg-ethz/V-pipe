FROM continuumio/miniconda:4.7.12

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
 && rm -rf /var/lib/apt/lists/*

RUN curl -O "https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh" \
    && bash quick_install.sh -b master -p /V-pipe_source -w /V-pipe_workdir

WORKDIR /V-pipe_workdir
RUN \
    ./vpipe -j 1 --conda-create-envs-only

VOLUME /V-pipe
WORKDIR /V-pipe

ENTRYPOINT ["/V-pipe_workdir/vpipe"]
