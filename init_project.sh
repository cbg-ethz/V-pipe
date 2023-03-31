#!/usr/bin/env bash

VPIPE_DIR="$(dirname "$0")"

if [[ -r "${VPIPE_DIR}/init_project.conf" ]]; then
   # shellcheck disable=SC1091
   . "${VPIPE_DIR}/init_project.conf" ''
fi

USAGE="
usage: $0 [options]

-m         bootstrap only a minimal set of files (vpipe_config.yaml and vpipe wrapper)
-n         disable auto-detection and management of conda environments
-c <PATH>  path where conda is installed
-e <ENV>   conda environment to use
-d         debug by dumping the conda search paths
-b         include wrappers for benchmark runs
-o         overwrite config.yaml even if there is one existing
-h         print this help message and exit
"

while getopts 'mnc:e:bdoh' opt; do
    case "${opt}" in
        m)
            MINIMAL=1;
            ;;
        n)
            NOCONDAAUTODETECT=1;
            ;;
        c)
            vp_conda_path="${OPTARG}"
            ;;
        e)
            vp_conda_env="${OPTARG}"
            ;;
        b)
            BENCHMARK=1;
            ;;
        d)
            debug=1
            ;;
        o)
            overwrite=1
            ;;
        h)
            printf "%s\\n" "$USAGE"
            exit 0
            ;;
        *)
            printf "%s\\n" "$USAGE" > /dev/stderr
            exit 2
            ;;
    esac
done

PROJECT_DIR="$(pwd)"
if [[ ! -e "$PROJECT_DIR/config.yaml" || -n "${overwrite}" ]]; then
    if [[ -e "$PROJECT_DIR/config.yaml" ]]; then
        echo "Warning: there is already a $PROJECT_DIR/config.yaml file. Backing up:" > /dev/stderr
        mv -v "$PROJECT_DIR/config.yaml"{,.backup}
    fi
    sed $'s@^output:@input:\\\n    samples_file: samples.tsv\\\n\\\noutput:\\\n    datadir: samples/\\\n@' "$VPIPE_DIR/config/config.yaml" > "$PROJECT_DIR/config.yaml"
else
    echo "Warning: there is already a $PROJECT_DIR/config.yaml file. Use option '-o' to overwrite." > /dev/stderr
fi

# guess activation command
ACTIVATE=
EXTRA_VPIPE_OPTS=
if [ -z "$NOCONDAAUTODETECT" ]; then
    conda_search=( )
    if [[ -n "${vp_conda_path}" ]]; then
        # only search where explicitely asked to
        conda_search=( "${vp_conda_path}" )
    else
        # search order:
        # - 'V-pipe' conda environment
        # - base conda environment
        conda_search=( "$VPIPE_DIR/../"{mambaforge,miniconda3}{/envs/V-pipe,} )
        # - also search any currently active environment
        if [[ -x $(which conda) ]]; then
            conda_search+=( "$(conda info --base)" )
        elif [[ -n "${CONDA_EXE}" ]]; then
            conda_search+=( "${CONDA_EXE%/bin/conda}" )
        fi
        if [[ -n "${CONDA_PREFIX}" &&  "${CONDA_PREFIX}" =~ /envs/(^[^/]+)/ ]]; then
            conda_search+=( "${CONDA_PREFIX}" )
            if [[ -z "${vp_conda_env}" ]]; then
                vp_conda_env="${BASH_REMATCH[1]}"
            fi
        fi
    fi

    # debug the search
    if [[ -n "${debug}" ]]; then
        printf -- " - %s\n" "${conda_search[@]}"
    fi

    # Search !
    for VPIPEENV in "${conda_search[@]}"; do
        if [ -d "${VPIPEENV}" ] && [ -x "${VPIPEENV}/bin/activate" ] && [ -x "${VPIPEENV}/bin/conda" ]; then
            # Search the specified environment (if provided)
            if [[ -n "${vp_conda_env}" && -x "${VPIPEENV}/envs/${vp_conda_env}/bin/snakemake" ]]; then
                echo "Conda environment found in ${VPIPEENV}, environment <${vp_conda_env}>"
                ACTIVATE=". ${VPIPEENV}/bin/activate '${vp_conda_env}'"
                break
            fi

            # Explicit specific place? => skip autodetection
            if [[ -n "${vp_conda_env}" &&  -n "${vp_conda_path}" ]]; then
                break
            fi

            # Continue autodetection
            if [ -x "${VPIPEENV}/envs/V-pipe/bin/snakemake" ]; then
                echo "Conda environment found in ${VPIPEENV}, environment V-pipe"
                ACTIVATE=". ${VPIPEENV}/bin/activate 'V-pipe'"
                break
            elif [ -x "${VPIPEENV}/bin/snakemake" ]; then
                echo "Conda environment found in ${VPIPEENV}"
                ACTIVATE=". ${VPIPEENV}/bin/activate 'base'"
                break
            fi
        fi
    done
    if [ -z "${ACTIVATE}" ]; then
        echo "Warning: cannot detect conda environment" 1>&2
    fi
    EXTRA_VPIPE_OPTS="--use-conda"
else
    echo "Activate the appropriate conda environment before use"
fi


cat > "$PROJECT_DIR/vpipe" <<EOF
#!/usr/bin/env bash
${ACTIVATE:+"# shellcheck disable=SC1091
${ACTIVATE}"}
exec -a "\$0" snakemake -s "$VPIPE_DIR/workflow/Snakefile" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
chmod +x "$PROJECT_DIR/vpipe"

if [[ -n "$BENCHMARK" ]]; then
    cat > "$PROJECT_DIR/vpipeBench" <<EOF
#!/usr/bin/env bash
${ACTIVATE:+"# shellcheck disable=SC1091
${ACTIVATE}"}
exec -a "\$0" snakemake -s "$VPIPE_DIR/resources/auxiliary_workflows/benchmark/vpipeBench.snake" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
    chmod +x "$PROJECT_DIR/vpipeBench"

    cat > "$PROJECT_DIR/vpipeBenchRunner" <<EOF
#!/usr/bin/env bash
${ACTIVATE:+"# shellcheck disable=SC1091
${ACTIVATE}"}
exec -a "\$0" snakemake -s "$VPIPE_DIR/resources/auxiliary_workflows/benchmark/vpipeBenchRunner.snake" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
    chmod +x "$PROJECT_DIR/vpipeBenchRunner"
fi

cat <<EOF
V-pipe project initialized!

Create and populate 'samples' directory and/or adjust config.yaml.
Then, use ./vpipe to run V-pipe.
EOF


if [ -z "$MINIMAL" ]; then
    # NOTE currently all the necessary files are packaged as resources. Here is an example how to provide extra files to a project
    : #mkdir -p "$PROJECT_DIR/references"
    : #cp -riv "$VPIPE_DIR/references"/* "$PROJECT_DIR/references/"
fi
