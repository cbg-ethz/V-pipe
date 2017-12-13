
sortbam() {
	STEM=${2%.bam}
	RAND_NUM="${RANDOM}"

	OLD="$(dirname ${2})/old_${RAND_NUM}.bam"
	TEMP="$(dirname ${2})/temp_${RAND_NUM}"

	mv "${2}" "${OLD}"
	echo "Sorting BAM file"
	$1 sort -o "${2}" -O bam -T "${TEMP}" "${OLD}"

	echo "Indexing BAM file"
	samtools index "${2}"
	rm "${OLD}"
}

sam2bam() {
	STEM=${2%.sam}
	echo "Writing BAM file"
	$1 view -bS "${STEM}.sam" > "${STEM}.bam"
	sortbam $1 "${STEM}.bam"
}

bam2sam() {
	STEM=${2%.bam}
	$1 view -h "${STEM}.bam" > "${STEM}.sam"
}

InDelFixer() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx10G -jar $@
}

ConsensusFixer() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx10G -jar $@
}

QuasiRecomb() {
	java -XX:+UseParallelGC -XX:NewRatio=9 -Xms2G -Xmx55G -jar $@
}

SamToFastq() {
	java -jar $1 SamToFastq ${@:2}
}
