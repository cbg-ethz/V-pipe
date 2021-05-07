#!/bin/bash
#
# build_docker_container.sh
# Copyright (C) 2021 Uwe Schmitt <uwe.schmitt@id.ethz.ch>
#
# Distributed under terms of the MIT license.
#

# https://stackoverflow.com/questions/59895/
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

BRANCH=$(git -C ${SCRIPT_DIR} branch --show-current)

docker build --progress plain -t vpipe-${BRANCH} -f ${SCRIPT_DIR}/Dockerfile ..
