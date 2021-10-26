rule write_summary_json:
    input:
        lambda wildcard: [f for f in all_files if f.endswith("ref_ambig_dels.fasta")]
    output:
        "summary.json"
    conda:
        config.report_sequences["conda"]
    shell:
        # needed {CONDA_PREFIX} on my dev setup on Mac + pyenv:
        """
        ${{CONDA_PREFIX}}/bin/python {VPIPE_BASEDIR}/scripts/report_sequences.py {input}
        """
