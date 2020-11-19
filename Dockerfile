FROM debian:stable

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
 && rm -rf /var/lib/apt/lists/*

RUN curl -O "https://raw.githubusercontent.com/cbg-ethz/V-pipe/master/utils/quick_install.sh" \
    && bash quick_install.sh -b master -p /V-pipe_source -w /V-pipe_workdir

WORKDIR /V-pipe_workdir
RUN ln -s /V-pipe_source/V-pipe/testdata/ samples \
    && ./vpipe -j 1 --conda-create-envs-only --conda-prefix /conda_prefix

VOLUME /V-pipe
WORKDIR /V-pipe

ENTRYPOINT ["/V-pipe_workdir/vpipe", "--conda-prefix", "/conda_prefix"]
