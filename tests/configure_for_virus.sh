#!/usr/bin/env bash

if [[ -z "$1" || "$1" == '-h' || "$1" == '--help' ]]; then
	echo "usage:	$0 <virus> [threads]

Description:
	creates configuraiton and samples directory for virus in current working directory

Arguments:
	virus	- name of the virus to test [ 'hiv', 'sars-cov-2' ]
	threads	- overides envrionment variable THREADS [ default: 4 ]"
	exit 0
fi

virus="$1"
threads="${2:-${THREADS:-4}}"
trim_primers="false"

echo "Configuring for ${virus}, using ${threads} threads"

# samples data
data_root="tests/data/${virus}"
mkdir -p samples
cp --link -vrf "${data_root}"/*/ samples/

if [ -e "${data_root}/samples.tsv" ]; then
	cp -vrf "${data_root}/samples.tsv" samples/
	# automatically turn trimming on if 4-columns format in TSV
	# shellcheck disable=SC2162
	while read s b l p o; do
		# no proto?
		if [[ -z "${p}" ]]; then
			continue
		fi

		# proto => trim!
		trim_primers="true"
		echo "with trimming"
		break

		: "${s} ${b} ${l} ${o}" are unused
	done < "samples/samples.tsv"
fi

# configuration file
cat > config/config.yaml <<CONFIG
general:
  virus_base_config: "${virus}"

output:
  snv: true
  local: true
  global: false
  visualization: true
  diversity: true
  QA: true
  upload: true
  trim_primers: ${trim_primers}

upload:
  orig_cram: true

snv:
  threads: ${threads}
CONFIG

# does this test data provides extra configuration options
if [ -e "${data_root}/config-extra.yaml" ]; then
	echo "with extra config:"
	cat "${data_root}/config-extra.yaml"
	echo
	# recursively merge using go-yq
	# shellcheck disable=SC2016
	yq eval-all --inplace '. as $item ireduce ({}; . * $item )' config/config.yaml "${data_root}/config-extra.yaml"
fi

# Display config
echo
echo "Configuration:"
echo
yq config/config.yaml
