#!/bin/bash
set -e

THREADS=5

# https://stackoverflow.com/questions/59895/
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
VPIPEROOT=${HERE}/..

if [[ -z "${1}" ]]; then
	echo "Usage: $0 <virus>" >&2
	exit 2
elif [[ ! "${1}" =~ ^[[:alnum:]_-]+$ ]]; then
	echo "Bad virus name $1" >&2
	exit 2
elif [[ ! -f ${VPIPEROOT}/config/${1}.yaml ]]; then
	echo "Missing virus base config for ${1}" >&2
	exit 2
elif [[ ! -d ${VPIPEROOT}/tests/data/${1}/ ]]; then
	echo "Missing virus test data for ${1}" >&2
	exit 2
fi
VIRUS=$1

CWD=$(pwd)
function restore_wd {
    cd ${CWD}
}
trap restore_wd EXIT

PROJECT_DIR=/tmp/project/${VIRUS}

function setup_project {
    PROJECT_DIR=$(mktemp -d)
    pushd ${PROJECT_DIR}
    ${VPIPEROOT}/init_project.sh
    mkdir samples
    cp -R ${VPIPEROOT}/tests/data/${VIRUS}/* samples/
    [ -e samples/samples.tsv ] && mv samples/samples.tsv ./
    popd
}

# setup project files when not run on via github actions
[ x$CI == x ] && setup_project


function run_workflow {

    pushd ${PROJECT_DIR}
    cat > config.yaml <<CONFIG
general:
    virus_base_config: "${VIRUS}"

output:
    snv: true
    local: true
    global: false
    visualization: true
    QA: true

snv:
    threads: ${THREADS}
CONFIG

    PYTHONUNBUFFERED=1 snakemake -s ${VPIPEROOT}/vpipe.snake --configfile config.yaml --use-conda --cores ${THREADS} --dry-run
    echo
    cat samples.tsv
    echo
    PYTHONUNBUFFERED=1 snakemake -s ${VPIPEROOT}/vpipe.snake --configfile config.yaml --use-conda --cores ${THREADS} -p --keep-going
    popd
}


TEST_NAME=$(basename ${0%.*})_${VIRUS}
EXIT_CODE=0

function check_logs {
    grep -E 'failed|for error' ${PROJECT_DIR}/.snakemake/log/*.snakemake.log && EXIT_CODE=1 || echo "snakemake execution successful"
}

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
check_logs
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
