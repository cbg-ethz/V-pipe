#!/usr/bin/env bash

# defaults
PREFIX=$(pwd)
FORCE=
BRANCH=master
RELEASE=
WORKDIR=
MINIMAL=

# Helper
fail() {
	printf '\e[31;1mArgh: %s\e[0m\n'	"$1"	1>&2
	[[ -n "$2" ]] && echo "$2" 1>&2
	exit 1
}

oops() {
	printf '\e[33;1mOops:%s\e[0m\n'	"$1"	 1>&2
}

title() {
	printf '\e[34;1m======================\n%s\n======================\e[0m\n\n'	"$1"
}

message() {
	printf '\e[37;1m%s\t%s\e[0m\n' "$1" "$2"
}

status() {
	printf '\e[36;1m%s\e[0m\n'	"$1"
}

check_directory() {
	if [[ -d "$1" ]]; then
		if (( FORCE )); then
			rm -rf "$1"
		else
			fail "${2:-Directory} ${1} already exists" 'check and use -f if needed'
		fi
	fi
}

usage() {
echo "usage: $0 [options]
options:
-f           force overwriting directories
-p PREFIX    prefix directory under which to install V-pipe
             [default: current directory]
-b BRANCH    install specified branch of V-pipe${BRANCH:+
             [default: $BRANCH]}
-r RELEASE   install specified release package of V-pipe${RELEASE:+
             [default: $RELEASE]}
-w WORKDIR   create and populate working directory
-m           only minimal working directory
-h           print this help message and exit"
}


# parameters
while getopts ":fp:b:r:w:mh" opt; do
	case ${opt} in
		f)
			FORCE=1
		;;
		p)
			PREFIX="${OPTARG}"
		;;
		b)
			BRANCH="${OPTARG}"
		;;
		r)
			RELEASE="${OPTARG}"
		;;
		w)
			WORKDIR="${OPTARG}"
		;;
		m)
			MINIMAL='-m'
		;;
		h)
			usage
			exit 0
		;;
		\?)
			fail "Invalid option: ${OPTARG}"
		;;
		:)
			fail "Invalid option: ${OPTARG} requires an argument"
		;;
	esac
done
shift $((OPTIND -1))

[[ -z "${BRANCH}" && -z "${RELEASE}" ]] && fail "Please specify either a branch or a release package"


###################
#                 #
#   Basic stuff   #
#                 #
###################

DOWNLOAD=
if [[ -x $(which wget) ]] || wget --version; then # most Linux distros (Gentoo)
	DOWNLOAD='wget'
elif [[ -x $(which curl) ]] || curl --version; then # a few Linux distros (CentOS, Archlinux) and Mac OS X
	DOWNLOAD='curl -O'
else
	fail 'Please install either wget or curl'
fi;

# HACK if not available locally, we will install git from conda later on (together with the other conda packages: snakemake, etc.)
GIT=
if [[ -x $(which git) ]] && git --version 2>/dev/null >/dev/null; then # most computers with development tools installed, clusters, etc
	GIT=
else # out-of-the box Mac OS X and some Linux dockers
	oops 'git is missing, I will download it from conda'
	GIT=git
fi;

# CHECK having environment modifiers (such as conda or modules) is a *horrendously bad idea* that can cause hard to understand errors much later
ENVIRONMENTWARNING=
for PROFILE in $HOME/.bash_profile $HOME/.bashrc $HOME/.profile; do
	if [[ -e $PROFILE ]] &&  grep -H 'conda initialize\|CONDA\|module \(add\|load\)' $PROFILE; then
		ENVIRONMENTWARNING=1
	fi
done
if [[ -n "$CONDA_PREFIX" ]]; then
	echo 'CONDA_PREFIX environment variable set'
	ENVIRONMENTWARNING=1
fi
if (( ENVIRONMENTWARNING )); then
	oops 'You have conda or modules automatically loaded in your profile. This is can lead to potential conflicts and errors.'
	echo 'consider always loading such environment modifiers either manually or at the beginning of your jobs, never in your profile.'
	sleep 3
fi


##################
#                #
#   Miniconda3   #
#                #
##################

MINICONDA=
MINICONDAPATH="${PREFIX}/miniconda3"

title 'installing Miniconda3'

# Check if directory is free
check_directory "${MINICONDAPATH}" 'Miniconda3 installation path'

# Check OS for OS-Spefic Miniconda3 installer
if [[ "$OSTYPE" == linux* ]]; then 
	MINICONDA=Miniconda3-latest-Linux-x86_64.sh 
elif [[ "$OSTYPE" == darwin* ]]; then
	MINICONDA=Miniconda3-latest-MacOSX-x86_64.sh
else
	fail 'I cannot detect OS. Only Linux and Mac OS X are supported' 'manually override OSTYPE environment variable if needed'
fi

message 'Using installer:' ${MINICONDA}

# Get and install Miniconda3
mkdir -p ${PREFIX}
cd ${PREFIX}
[[ -f ${MINICONDA} ]] && rm ${MINICONDA}
${DOWNLOAD} https://repo.anaconda.com/miniconda/${MINICONDA}
sh ${MINICONDA} -b -p miniconda3
# -b for batch (no question asked)


. miniconda3/bin/activate

# set the channel precedence (lowest to highest)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge 
# NOTE conda-forge *HAS TO* be higher than bioconda

VPIPEENV=
if conda install --yes snakemake-minimal $GIT; then	# NOTE Mac OS X and some Linux dockers don't have git out of the box
	: # success!
else
	oops 'I cannot install snakemake in base environment. Conflicts ?'

	VPIPEENV=V-pipe
	# HACK Alternate to consider if we have have version conflicts
	conda create --yes -n ${VPIPEENV} -c conda-forge -c bioconda snakemake-minimal conda git || fail "I cannot install snakemake in environment ${VPIPEENV}."
	conda activate ${VPIPEENV}
fi

# NOTE No need to download and install V-pipe dependencies, snakemake --use-conda handles this.

echo $'\n'



##############
#            #
#   V-pipe   #
#            #
##############

title 'installing V-pipe'

if [[ -z "${RELEASE}" ]]; then
	message 'Using branch:' ${BRANCH}

	check_directory 'V-pipe' 'V-pipe installation directory'
	git clone --branch ${BRANCH} https://github.com/cbg-ethz/V-pipe.git || fail "I cannot install branch ${BRANCH}."
else
	message 'Using release:' "${RELEASE}"
	check_directory "V-pipe-${RELEASE}" 'V-pipe installation directory'
	${DOWNLOAD} "https://github.com/cbg-ethz/V-pipe/archive/${RELEASE}.tar.gz" || fail "I cannot download package ${BRANCH}."
	tar xvzf "${RELEASE}.tar.gz" || fail "I cannot install package ${RELEASE}."
fi

echo $'\n'



###############
#             #
#   Working   #
#             #
###############

INIT="$(pwd)/V-pipe${RELEASE:+-${RELEASE}}/init_project.sh"

if [[ -z "${WORKDIR}" ]]; then
	status 'Installation of V-pipe completed'
	echo "
To setup working directory:
	mkdir -p working
	cd working
	${INIT}"

	exit 0
fi

title 'Working directory'
message 'Working directory:' ${WORKDIR}

check_directory "${WORKDIR}" 'Working directory'
mkdir -p "${WORKDIR}" || fail "I cannot create directory"
cd "${WORKDIR}"
"${INIT}" ${MINIMAL} || fail "Populating working directory failed"

echo $'\n'
status 'Installation of V-pipe completed'

exit 0
