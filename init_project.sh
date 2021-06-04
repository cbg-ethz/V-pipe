#!/usr/bin/env bash

USAGE="
usage: $0 [options]

-m           bootstrap only a minimal set of files (vpipe_config.yaml and vpipe wrapper)
-n           disable auto-detection and management of conda environments
-b           include wrappers for benchmark runs
-h           print this help message and exit
"

args=$(getopt 'mnbh' "$*")
if [ "$?" != 0 ]; then
    printf "%s\\n" "$USAGE"
    exit 2
fi
eval set -- "$args"

for i; do
    case "$i" in
        -m)
            MINIMAL=1;
            shift
            ;;
        -n)
            NOCONDAAUTODETECT=1;
            shift
            ;;
        -b)
            BENCHMARK=1;
            shift
            ;;
        -h)
            printf "%s\\n" "$USAGE"
            exit 2
            shift
            ;;
        --)
            shift
            break
            ;;
    esac
done

PROJECT_DIR=$(pwd)

# https://stackoverflow.com/a/242550
VPIPE_DIR=$(dirname "$0")

cp -iv "$VPIPE_DIR/vpipe_config.yaml" "$PROJECT_DIR/"

# guess activation command
ACTIVATE=
EXTRA_VPIPE_OPTS=
if [ -z "$NOCONDAAUTODETECT" ]; then
    # search order:
    # - 'V-pipe' conda environment
    # - base conda environment
    for VPIPEENV in "$VPIPE_DIR/../miniconda3"{/envs/V-pipe,}; do
        if [ -d "${VPIPEENV}" ] && [ -x "${VPIPEENV}/bin/activate" ] && [ -x "${VPIPEENV}/bin/conda" ] && [ -x "${VPIPEENV}/bin/snakemake" ]; then
            echo "Conda environment found in ${VPIPEENV}"
            if [ -x "${VPIPEENV}/envs/V-pipe/bin/snakemake" ]; then
                ACTIVATE=". ${VPIPEENV}/bin/activate 'V-pipe'"
            else
                ACTIVATE=". ${VPIPEENV}/bin/activate 'base'"
            fi
            break
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
${ACTIVATE}
exec -a "\$0" snakemake -s "$VPIPE_DIR/vpipe.snake" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
chmod +x "$PROJECT_DIR/vpipe"

if [[ -n "$BENCHMARK" ]]; then
    cat > "$PROJECT_DIR/vpipeBench" <<EOF
    #!/usr/bin/env bash
    ${ACTIVATE}
    exec -a "\$0" snakemake -s "$VPIPE_DIR/vpipeBench.snake" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
    chmod +x "$PROJECT_DIR/vpipeBench"

    cat > "$PROJECT_DIR/vpipeBenchRunner" <<EOF
    #!/usr/bin/env bash
    ${ACTIVATE}
    exec -a "\$0" snakemake -s "$VPIPE_DIR/vpipeBenchRunner.snake" ${EXTRA_VPIPE_OPTS} "\$@"
EOF
    chmod +x "$PROJECT_DIR/vpipeBenchRunner"
fi

cat <<EOF
V-pipe project initialized!

Create and populate 'references' and 'samples' directories or adjust vpipe_config.yaml.
Then, use ./vpipe to run V-pipe.
EOF


if [ -z "$MINIMAL" ]; then
    mkdir -p "$PROJECT_DIR/references"
    cp -riv "$VPIPE_DIR/references"/* "$PROJECT_DIR/references/"
fi
