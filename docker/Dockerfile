FROM debian:stable

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    curl \
 && rm -rf /var/lib/apt/lists/*


RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && sh Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3 \
    && /miniconda3/bin/conda config --add channels defaults \
    && /miniconda3/bin/conda config --add channels bioconda \
    && /miniconda3/bin/conda config --add channels conda-forge \
    && /miniconda3/bin/conda install --yes snakemake-minimal

COPY . /V-pipe
RUN mkdir /work

WORKDIR /work

COPY docker/setup_conda.sh /V-pipe
RUN /V-pipe/setup_conda.sh

COPY docker/entry_point.sh /V-pipe

ENTRYPOINT ["/V-pipe/entry_point.sh"]
