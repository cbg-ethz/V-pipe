import os

rule compute_diversity_measures:
    input:
        fnames_samples_snvs_vcf=expand('{dataset}/variants/SNVs/snvs.vcf', dataset=datasets),
        fname_ref = reference_file,
    output:
        diversity_csv=expand('{dataset}/variants/SNVs/diversity_measures.csv', dataset=datasets),
        shannon_csv=expand('{dataset}/variants/SNVs/position_shannon_entropy.csv', dataset=datasets),
    params:
        output_dir=expand('{dataset}/variants/SNVs/', dataset=datasets),
        script_dir=os.path.join(VPIPE_BASEDIR, "scripts"),
    run:
        for i in range(len(params.output_dir)):
            fname_sample_vcf = input.fnames_samples_snvs_vcf[i]
            out_dir = params.output_dir[i]
            shell("python3 {params.script_dir}/compute_diversity_measures.py {fname_sample_vcf} {input.fname_ref} {out_dir}")
