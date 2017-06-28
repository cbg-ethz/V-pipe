
sortbam() {
	STEM=${1%.bam}
	RAND_NUM="${RANDOM}"

	OLD="$(dirname ${1})/old_${RAND_NUM}.bam"
	TEMP="$(dirname ${1})/temp_${RAND_NUM}"

	mv "${1}" "${OLD}"
	echo "Sorting BAM file"
	samtools sort -o "${1}" -O bam -T "${TEMP}" "${OLD}"

	echo "Indexing BAM file"
	samtools index "${1}"
	rm "${OLD}"
}

sam2bam() {
	STEM=${1%.sam}
	echo "Writing BAM file"
	samtools view -bS "${STEM}.sam" > "${STEM}.bam"
	sortbam "${STEM}.bam"
}

bam2sam() {
	STEM=${1%.bam}
	samtools view -h "${STEM}.bam" > "${STEM}.sam"
}

InDelFixer() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx10G -jar /cluster/work/bewi/modules/InDelFixer/current/InDelFixer.jar $@
}

ConsensusFixer() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx10G -jar /cluster/work/bewi/modules/ConsensusFixer/current/ConsensusFixer.jar $@
}

QuasiRecomb() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx55G -jar /cluster/work/bewi/modules/QuasiRecomb/QuasiRecomb.jar $@
}

SamToFastq() {
	java -jar /cluster/work/bewi/modules/picard/current/picard.jar SamToFastq $@
}
