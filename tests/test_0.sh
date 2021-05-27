#!/bin/bash
set -e

# https://stackoverflow.com/questions/59895/
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
VPIPEROOT=${HERE}/..

CWD=$(pwd)
function restore_wd {
    cd ${CWD}
}
trap restore_wd EXIT

PROJECT_DIR=/tmp/project

function setup_project {
    PROJECT_DIR=$(mktemp -d)
    pushd ${PROJECT_DIR}
    ${VPIPEROOT}/init_project.sh
    mkdir samples
    cp -R ${VPIPEROOT}/testdata/pos_M* samples
    popd
}

# setup project files when not run on via github actions
[ x$CI == x ] && setup_project


function run_workflow {

    pushd ${PROJECT_DIR}
    PYTHONUNBUFFERED=1 snakemake -s ${VPIPEROOT}/vpipe.snake --use-conda --dry-run
    echo
    cat samples.tsv
    echo
    PYTHONUNBUFFERED=1 snakemake -s ${VPIPEROOT}/vpipe.snake --use-conda -p -j 2
    popd
}


TEST_NAME=$(basename ${0%.*})
EXIT_CODE=0

DIFF_FILE=/tmp/diffs_${TEST_NAME}.txt
LOG_FILE=/tmp/log_${TEST_NAME}.txt
rm -f ${DIFF_FILE}
rm -f ${LOG_FILE}

function compare_to_recorded_results {

    cd ${CWD}/expected_outputs/${TEST_NAME}

    for RECORDED_OUTPUT in $(find . -type f); do
        CURRENT_OUTPUT=${PROJECT_DIR}/${RECORDED_OUTPUT}
        echo COMPARE ${RECORDED_OUTPUT} AND ${CURRENT_OUTPUT}
        if diff -I '^#' ${RECORDED_OUTPUT} ${CURRENT_OUTPUT} >> ${DIFF_FILE}; then
            :
        else
            echo
            echo RESULTS ${RECORDED_OUTPUT} AND ${CURRENT_OUTPUT} DIFFER
            echo
            EXIT_CODE=1;
        fi
    done
}

# setup_project
run_workflow 2>&1 | tee ${LOG_FILE}
echo
echo
compare_to_recorded_results

echo

if [ ${EXIT_CODE} = 1 ]; then
    echo TESTS FAILED, CHECK ${DIFF_FILE} and ${LOG_FILE} FOR FUTHER INFORMATION
else
    echo TESTS SUCEEDED
fi

exit ${EXIT_CODE}
