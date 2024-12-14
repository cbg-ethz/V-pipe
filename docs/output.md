> **Warning:** the documentation on this wiki is deprecated and refers to V-pipe versions 1.0 and 2.0 (published in [Bioinformatics](https://doi.org/10.1093/bioinformatics/btab015)).
> You can find the most recent documentation on running V-pipe 3.0 in the [readme on the _master_ branch](https://github.com/cbg-ethz/V-pipe/blob/master/README.md).

## Output files

For each sample, V-pipe produces several output files that are located in the corresponding sample-specific directory. First, the alignment file and consensus sequences are located in the `alignments` and `references` subdirectories, respectively. Second, output files containing SNVs and viral haplotypes are located in the `variants` subdirectories.

Below, we provide an example of relevant output files and their locations, following the same structure as in the [getting-started](https://github.com/cbg-ethz/V-pipe/wiki/getting-started) section. The output files for the two patient samples will be located in the following subdirectories:

```
working_directory
├─references
│  └───HXB2.fasta
└─samples
  ├──patient1
  │   ├──20100113
  │   │   ├──alignments
  │   │   |  └──REF_aln.bam
  │   │   ├──references
  |   |   |  ├──ref_ambig.fasta
  |   |   |  └──ref_majority.fasta
  |   |   └──variants
  |   |      ├──SNVs
  |   |      |  └──snvs.vcf
  |   |      └──global
  |   |         └──contigs_stage_c.fasta
  │   └──20110202
  │      ├──alignments
  │       |  └──REF_aln.bam
  │       ├──references
  |       |  ├──ref_ambig.fasta
  |       |  └──ref_majority.fasta
  |       └──variants
  |          ├──SNVs
  |          |  └──snvs.vcf
  |          └──global
  |             └──contigs_stage_c.fasta         
  └──patient2
          ├──alignments
          |  └──REF_aln.bam
          ├──references
          |  ├──ref_ambig.fasta
          |  └──ref_majority.fasta
          └──variants
             ├──SNVs
             |  └──snvs.vcf
             └──global
                └──contigs_stage_c.fasta
```

In addition, V-pipe generates a csv file containg the frequencies of all minor alleles that differ from the consensus among analysed samples. This ouput file is located in the `variants` subdirectory

```
working_directory
└─variants
  └───minority_variants.tsv
```


