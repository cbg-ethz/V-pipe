import os


__author__ = "Kim"
__author__ = "Lara Fuhrmann"
__license__ = "Apache2.0"
__maintainer__ = "Ivan Topolsky"
__email__ = "v-pipe@bsse.ethz.ch"


rule compute_diversity_measures:
    input:
        fnames_samples_snvs_vcf="{dataset}/variants/SNVs/snvs.vcf",
        fname_ref=reference_file,
    output:
        diversity_csv="{dataset}/variants/SNVs/diversity_measures.csv",
        shannon_csv="{dataset}/variants/SNVs/position_shannon_entropy.csv",
    conda:
        config.diversity["conda"]
    benchmark:
        "{dataset}/variants/SNVs/diversity_measures.benchmark"
    script:
        "../scripts/compute_diversity_measures.py"


rule aggregate_diversity:
    input:
        fnames_diversity=expand(
            "{dataset}/variants/SNVs/diversity_measures.csv", dataset=datasets
        ),
        fnames_shannon=expand(
            "{dataset}/variants/SNVs/position_shannon_entropy.csv", dataset=datasets
        ),
    benchmark:
        cohortdir("diversity_measures.benchmark")
    output:
        diversity_csv=cohortdir("aggregated_diversity.csv"),
        shannon_csv=cohortdir("aggregated_entropy.csv"),
    conda:
        config.diversity["conda"]
    script:
        "../scripts/aggregate_diversity.py"
