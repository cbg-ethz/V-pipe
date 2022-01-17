
rule write_summary_json:
    input:
        consensus = lambda wildcard: [f for f in all_files if f.endswith("ref_majority_dels.fasta")],
        dehumanized = lambda wildcard: [f for f in all_files if f.endswith("dehuman.cram")],
        all_files = all_files
    output:
        "summary.zip",
    conda:
        config.report_sequences["conda"]
    shell:
        # needed {CONDA_PREFIX} on my dev setup on Mac + pyenv:
        """
        ${{CONDA_PREFIX}}/bin/python {VPIPE_BASEDIR}/scripts/report_sequences.py {input.consensus}
        """
