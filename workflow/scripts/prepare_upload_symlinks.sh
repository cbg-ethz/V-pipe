#!/usr/bin/env bash

#
# Example script to prepare and assist uploads.
#  - it will get called whenever the consensus FASTA files or the CRAM-compressed raw-reads are updated by snakemake
#  - it will create a per-sample directory '.../uploads/' with symlinks to all uploadable files
#  - it will create a global directory 'uploads' with symlinks to samples that were updated
#    - these files aren't tracked by snakemake's DAG
#    - thus they can be deleted without triggering a re-build by snakemake
#    - but they will be re-created whenever an input file dependency changes
#    - it is possible to iteratively scan this global directory between runs of V-pipe to determine
#      which are new/updated samples to consider for upload
#  - meanwhile the OUTPUT file is tracked by the snakemake DAG
#    - it must be provided
#    - its destruction will retrigger the prepare_upload rule and this script
#    - it can optionally be used by the script to store arbitrary data (e.g. json)
#

do_random_nonce=0
exec_cmd=

usage() { echo "Usage: $0 [ -h ] [ -n ] [ -e <CMD> ] [ -- ] <OUTPUT> <SAMPLE_ID> <SAMPLE_DIR> [ <UPLOAD_FILES> ... ]

options:
	-r : append a random nonce at the end of the global symlinks,
	     multiple updates of the same sample and date will generate
	     multiple unique symlinks, one for each update.
	-R : no random nonce at the end of the global symlinks,
	     they will not be unique, a sample and date will always have
	     a single symlink, no matter how many time it was updated.
	[ Default: random nonce are $( if (( do_random_nonce )); then echo "enabled"; else echo "disabled"; fi ) ]

	-e : run <CMD> at the end of this script with the same positionals;
	       <OUTPUT> <SAMPLE_ID> <SAMPLE_DIR> [ <UPLOAD_FILES> ... ]
	     it becomes that script's job to create the output file,
	     if successful.

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
while getopts "rRe:h" o; do
	case "${o}" in
		r)	do_random_nonce=1	;;
		R)	do_random_nonce=0	;;
		e)	exec_cmd="${OPTARG}"	;;
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

( set -x; rm -f "uploads/$unique_id"; ln "-s${force}" "$fixed_uploads" "uploads/$unique_id" )


# run command if asked to

if [[ -n "${exec_cmd}" ]]; then
	exec ${exec_cmd} "${output}" "${sample_id}" "${sample_dir}" "${to_upload[@]}"
	echo "Failed to exec ${exec_cmd}" > /dev/stderr
	exit 2
fi

# create the mandatory output file
touch "${output}"
