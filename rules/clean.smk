__author__ = "Susana Posada-Cespedes"
__author__ = "David Seifert"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


rule extractclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/extracted_data
        """


rule trimmingclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/preprocessed_data
        """


rule vicunaclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/initial_consensus
        rm -rf {params.DIR}/*/*/references/vicuna_consensus.fasta
        rm -rf {params.DIR}/*/*/references/initial_consensus.fasta
        rm -rf references/initial_aln.fasta
        rm -rf references/initial_aln_gap_removed.fasta
        rm -rf references/MAFFT_initial_aln.*
        """


rule msaclean:
    shell:
        """
        rm -rf references/ALL_aln_*.fasta
        rm -rf references/MAFFT_*_cohort.*
        """


rule alignclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/alignments
        rm -rf {params.DIR}/*/*/QA_alignments
        rm -rf {params.DIR}/*/*/references/ref_ambig.fasta
        rm -rf {params.DIR}/*/*/references/ref_majority.fasta
        rm -rf {params.DIR}/*/*/references/initial_consensus.fasta
        """


rule bwaclean:
    input:
        "{}.bwt".format(reference_file),
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -f {input}
        rm -rf {params.DIR}/*/*/alignments
        rm -rf {params.DIR}/*/*/references/ref_ambig*.fasta
        rm -rf {params.DIR}/*/*/references/ref_majority*.fasta
        """


rule bowtieclean:
    input:
        INDEX1="{}.1.bt2".format(reference_file),
        INDEX2="{}.2.bt2".format(reference_file),
        INDEX3="{}.3.bt2".format(reference_file),
        INDEX4="{}.4.bt2".format(reference_file),
        INDEX5="{}.rev.1.bt2".format(reference_file),
        INDEX6="{}.rev.2.bt2".format(reference_file),
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -f {input}
        rm -rf {params.DIR}/*/*/alignments
        rm -rf {params.DIR}/*/*/references/ref_ambig*.fasta
        rm -rf {params.DIR}/*/*/references/ref_majority*.fasta
        """


rule snvclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/variants/SNVs
        """


rule savageclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/variants/global/contigs_stage_?.fasta
        rm -rf {params.DIR}/*/*/variants/global/stage_?
        """


rule haplocliqueclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm {params.DIR}/*/*/variants/global/quasispecies.*
        """


rule visualizationclean:
    params:
        DIR=config.input["datadir"],
    shell:
        """
        rm -rf {params.DIR}/*/*/visualization
        """
