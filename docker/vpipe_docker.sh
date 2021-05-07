#!/bin/bash
#
# vpipe_docker.sh
# Copyright (C) 2021 Uwe Schmitt <uwe.schmitt@id.ethz.ch>
#
# Distributed under terms of the MIT license.
#

# https://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

BRANCH=$(git -C ${SCRIPT_DIR} branch --show-current)

docker run -it --init -v $(pwd):/work/ vpipe-${BRANCH} vpipe "$@" 2>&1 | sed "s/\/work\/.snakemake/.snakemake/g"
