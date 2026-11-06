__author__ = "Ivan Topolsky"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


rule process_interactive:
    """Processes the sample to ndjson, skip upload to loculus.

    This rule demonstrates the flexible parameter approach:
    - Parameters can be provided via CLI arguments (as shown)
    - Or via environment variables (TIMELINE_FILE, PRIMER_FILE, NEXTCLADE_REFERENCE)
    - CLI arguments take precedence over environment variables

    Note: On ETH cluster, internet access is provided by eth_proxy module loaded at job level.
    """
    input:
        sample_fp=rules.paired_end_read_merger.output.fname_bam_merged,
        timeline_file=ancient(rules.tallymut.input.times),
    output:
        result_fp="{dataset}/interactive/{ref}.ndjson.zst",
    log:
        outfile="{dataset}/interactive/{ref}.ndjson.out.log",
        errfile="{dataset}/interactive/{ref}.ndjson.err.log",
    benchmark:
        "{dataset}/interactive/{ref}.ndjson.benchmark"
    conda:
        config.interactive["conda"]
    threads: 1
    resources:
        disk_mb=1250,
        mem_mb=config.interactive["mem"],
        runtime=config.interactive["time"],
    params:
        sample_id=ID,
        id_sep=config.general["id_separator"],
        organism=lambda wildcards: wildcards.ref,
        cfg_file="interactive_cfg/ref.{ref}.yaml",
        cache_path="interactive_cfg/cache/",
    shell:
        """
        XDG_CACHE_HOME='{params.cache_path}'

        lapis_url="$(yq -r .LAPIS_URL '{params.cfg_file}')"
        reference_accession="$(yq -r .REFERENCE_ACCESSION '{params.cfg_file}')"
        nuc_ref="$(yq -r .NUC_REF '{params.cfg_file}')"
        aa_ref="$(yq -r .AA_REF '{params.cfg_file}')"
        full_id='{params.sample_id}'

        nuc_ref="$(
    cd "$(dirname '{params.cfg_file}')"
    realpath "${{nuc_ref}}"
)"
        aa_ref="$(
    cd "$(dirname '{params.cfg_file}')"
    realpath "${{aa_ref}}"
)"

        sr2silo process-from-vpipe \
            --skip-merge \
            --input-file {input.sample_fp} \
            --sample-id "${{full_id%{params.id_sep}*}}" \
            --timeline-file {input.timeline_file} \
            --organism {params.organism} \
            --output-fp {output.result_fp} \
            --lapis-url "${{lapis_url}}" \
            --nuc-ref "${{nuc_ref}}" \
            --aa-ref "${{aa_ref}}" \
            --reference-accession "${{reference_accession}}" \
            >{log.outfile} 2>{log.errfile}
        """


def get_all_inter_refs(dir):
    if not dir:
        return []
    (REFS,) = glob_wildcards(os.path.join(dir, "ref.{ref}.yaml"))
    # return expand(os.path.join(dir, "{ref}.yaml"), ref=REFS)
    return REFS


all_inter_refs = get_all_inter_refs("interactive_cfg/")


rule allInteractive:
    input:
        expand(
            "{dataset}/interactive/{ref}.ndjson.zst",
            dataset=datasets,
            ref=all_inter_refs,
        ),
