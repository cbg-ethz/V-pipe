import os


rule compute_diversity_measures:
    input:
        fnames_samples_snvs_vcf="{dataset}/variants/SNVs/snvs.vcf",
        fname_ref=reference_file,
    output:
        diversity_csv="{dataset}/variants/SNVs/diversity_measures.csv",
        shannon_csv="{dataset}/variants/SNVs/position_shannon_entropy.csv",
    conda:
        config.diversity["conda"]
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
    output:
        diversity_csv=os.path.join(
            config.output["datadir"],
            config.output["cohortdir"],
            "aggregated_diversity.csv",
        ),
        shannon_csv=os.path.join(
            config.output["datadir"],
            config.output["cohortdir"],
            "aggregated_entropy.csv",
        ),
    conda:
        config.diversity["conda"]
    script:
        "../scripts/aggregate_diversity.py"
