#!/bin/sh

USAGE="
usage: $0 [options]

-m           bootstrap only a minimal set of files (vpipe.config and vpipe wrapper)
-h           print this help message and exit
"

args=$(getopt 'mh' "$*")
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

cp -iv "$VPIPE_DIR/vpipe.config" "$PROJECT_DIR/"

cat > "$PROJECT_DIR/vpipe" <<EOF
#!/bin/sh
snakemake -s "$VPIPE_DIR/vpipe.snake" "\$@"
EOF
chmod +x "$PROJECT_DIR/vpipe"

cat <<EOF
V-pipe project initialized!

Create and populate 'references' and 'samples' directories or adjust vpipe.config.
Then, use ./vpipe to run V-pipe.
EOF


if [ -z "$MINIMAL" ]; then
    mkdir -p "$PROJECT_DIR/references"
    cp -iv "$VPIPE_DIR/references"/* "$PROJECT_DIR/references/"
fi
