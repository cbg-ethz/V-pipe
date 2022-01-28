localrules:
    prepare_upload


rule prepare_upload:
    input:
        R1=partial(raw_data_file, pair=1),
        R2=partial(raw_data_file, pair=2),
        consensus_aligned="{dataset}/references/ref_majority_dels.fasta",
        consensus_unaligned="{dataset}/references/consensus_ambig.bcftools.fasta",
    output:
        upload_prepared_touch="{dataset}/upload_prepared.touch"
    threads: 1

    shell:
        """
        to_upload=({input})

        # depending on config dehuman might not exist, this is why we
        # did not add it to the input fields:
        to_upload+=({wildcards.dataset}/raw_data/dehuman.cram)

        for p in ${{to_upload[@]}}; do
            test -e $p || continue
            fixed_p=$(realpath --relative-to {wildcards.dataset}/uploads $p)
            ( set -x; ln -f -s $fixed_p {wildcards.dataset}/uploads )
        done

        mkdir -p uploads

        sample_id=$(basename $(dirname {wildcards.dataset}))
        fixed_uploads=$(realpath --relative-to uploads {wildcards.dataset}/uploads)

        # make unique symbolic link:
        random=$(dd if=/dev/urandom bs=30 count=1 2>/dev/null | sha1sum -b | cut -d" " -f 1)
        unique_id=${{sample_id}}__${{random}}

        ( set -x; ln -s $fixed_uploads uploads/$unique_id )

        touch {output.upload_prepared_touch}
        """
