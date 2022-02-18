#!/usr/bin/env bash

#
# Example script to prepare and assist uploads.
#  - it will get called whenever the consensus FASTA files or the CRAM-compressed raw-reads are updated by snakemake
#  - it will create a per-sample directory '.../uploads/' with symlinks to all uploadable files
#  - it will crate a global directory 'uploads' with symlinks to samples that were updated
#    - these file aren't tracked by snakemake's DAG
#    - thus they can be deleted without triggering a re-build by snakemake
#    - but they will be re-created whenever an input file dependency changes
#    - it is possible to iteratively scan this global directory between runs of V-pipe to determine
#      which are new/updated samples to consider for upload
#  - meanwhile the OUTPUT file is tracked by the snakemake DAG
#    - it must be provided
#    - its destruction will retrigger the prepare_upload rule and this script
#    - it can optionally be used by the script to store arbitrary data (e.g. json)
#

usage() { echo "Usage: $0 [ -h ] [ -n ] [ -- ] <OUTPUT> <SAMPLE_ID> <SAMPLE_DIR> [ <UPLOAD_FILES> ... ]

options:
	-n : no random nonce at the end of the global symlinks,
	     they will be not unique
	-- : end of options, start of positional parameters

positional parameters:
	<OUTPUT>    : the output file that must be created by the rule
	<SAMPLE_ID> : a string (with no path separator slashes) that uniquely
	             identifies the sample and the date
	<SAMPLE_DIR>: the base directory of the sample
	<UPLOAD_FILES>: a list of files to consider for upload

Generates symlinks that help tracking new and updated samples to consider
for upload. Serves also as a demo for the upload parameters of V-pipe." 1>&2; exit "$1"; }

# NOTE it is possible to have named options (e.g. with getops) before the named options begin
do_random_nonce=1
while getopts "nh" o; do
	case "${o}" in
		n)	do_random_nonce=0	;;
		h)	usage 0	;;
		*)	usage 1	;;
	esac
done
shift $((OPTIND-1))




## get the positional parameters

output="${1}"
sample_id="${2//\//_}"
sample_dir="${3}"
shift 3
to_upload=( "${@}" )



## process

# create and populate per-sample directory
mkdir -p "${sample_dir}/uploads/"

for p in "${to_upload[@]}"; do
	test -e "${p}" || continue
	fixed_p="$(realpath --relative-to "${sample_dir}/uploads/" "${p}")"
	( set -x; ln -f -s "$fixed_p" "${sample_dir}/uploads/" )
done


# create and populare global directory
mkdir -p uploads/

fixed_uploads="$(realpath --relative-to "uploads" "${sample_dir}/uploads/")"

# make unique symbolic link:
if (( do_random_nonce )); then
	read -r random o < <(dd if=/dev/urandom bs=30 count=1 2>/dev/null | sha1sum -b)
	unique_id="${sample_id}__${random}"
	force=""
else
	unique_id="${sample_id}"
	force="f"
fi

( set -x; ln "-s${force}" "$fixed_uploads" "uploads/$unique_id" )


# create the mandatory output file
touch "${output}"
