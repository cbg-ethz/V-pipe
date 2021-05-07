#!/bin/bash
#
# Author: Uwe Schmitt <uwe.schmitt@id.ethz.ch>


# to enable offline execution of v-pipe we must make sure that all needed
# conda packages are stored in the container.
# to achieve this we setup our own local conda channel at /local_conda_channel
#
# why using channels?
#
# 1. snakemake uses hashes for the names of the created conda environments which
#    is an snakemake implementation detail. So we cannot pre-create these
#    environments with correct names reliably.
#
# 2. snakemake offers the cli flag --conda-create-envs-only for creating these
#    conda environments without running the actual workflow, but this feature
#    appears to be broken.

set -e

CHANNEL=/local_conda_channel
PACKAGES=${CHANNEL}/noarch

mkdir -p ${PACKAGES}

. /miniconda3/bin/activate base

# make sure that all downloaded packages
# are cached in $PACKAGES
conda config --remove-key pkg_dirs || true
conda config --append pkgs_dirs ${PACKAGES}

# create all environments, this will populate
# the $PACKAGES
for envfile in /V-pipe/envs/*.yaml; do
    echo ${envfile}
    echo
    filename=$(basename ${envfile})
    envname=${filename%.yaml}
    conda env create -n ${envname} -f ${envfile}
    echo
done

# split source and binary packages as demanded
# by conda to implement our own local file system
# based chanel at $CHANNEL:
mkdir -p ${CHANNEL}/linux-64
mv ${PACKAGES}/*linux-64* ${CHANNEL}/linux-64

# create index for $CHANNEL
conda install conda-build
echo
echo
echo conda index ${CHANNEL}
conda index ${CHANNEL}
echo
echo

# remove venvs since we have our own store now and
# environments are not needed in the container
echo rm -rf /miniconda3/envs
rm -rf /miniconda3/envs

# configur conda  to use our local channel only:
conda config --remove-key channels || true
conda config --append channels file://${CHANNEL}
conda config --set offline True

# a bit more log output when running the pipeline
# with --verbose:
conda config --set verbosity 1
