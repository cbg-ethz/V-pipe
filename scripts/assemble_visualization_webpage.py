import os
import re
import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

import vcf
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def convert_vcf(fname):
    """Convert VCF to JSON."""
    output = []

    if os.path.getsize(fname) == 0:
        print(f'Empty VCF: "{fname}"')
        return output

    print(f'Parsing VCF: "{fname}"')
    with open(fname) as fd:
        vcf_reader = vcf.Reader(fd)

        for record in vcf_reader:
            output.append(
                {
                    "position": record.POS,
                    "reference": record.REF,
                    "variant": [v.sequence for v in record.ALT],
                    "frequency": round(np.mean(
                        [v for k, v in record.INFO.items() if k.startswith("Freq")]
                    ), 3),
                    "posterior": round(1 - 10**(-record.QUAL / 10), 3)
                }
            )

    return json.dumps(output)


def parse_gff(fname):
    """Convert GFF to map."""
    features = []

    print(f'Parsing GFF: "{fname}"')
    with open(fname) as fd:
        for record in GFF.parse(fd):
            for feature in record.features:
                features.append(
                    {
                        "id": record.id,
                        "type": feature.type,
                        "name": feature.qualifiers.get('Name', [feature.id])[0],
                        "start": int(feature.location.start),
                        "end": int(feature.location.end),
                    }
                )
    return features


def arrange_gff_data(features):
    """Add row number to each feature."""
    features.sort(key=lambda f: f["start"])

    rows = []
    for feature in features:
        if not rows:
            feature["row_cnt"] = 0
            rows.append([feature])
        else:
            found = False
            for idx, row in enumerate(rows):
                if row[-1]["end"] <= feature["start"]:
                    feature["row_cnt"] = idx
                    row.append(feature)
                    found = True
                    break
            if not found:
                feature["row_cnt"] = len(rows)
                rows.append([feature])

    return [item for row in rows for item in row]


def get_gff_data(gff_dir):
    """Returns a map with filename key and gff json data."""
    # Hardcode metainformation for the GFF file provided in the repository.
    gff_metainfo = {}
    gff_metainfo["Genes_NC_045512.2"] = "Genes (ORFs)"
    gff_metainfo[
        "Sars-Cov2_Mature_products"
    ] = "Mature processed protein products"
    gff_metainfo["Sars-Cov2_Uniprot_domains"] = "Uniprot domains"
    gff_metainfo[
        "Sars-Cov2_TM_domains"
    ] = "Uniprot predicted transmembrane regions"

    gff_map = {}
    for path in os.listdir(gff_dir):
        full_path = os.path.join(gff_dir, path)
        description = os.path.splitext(path)[0]
        if description in gff_metainfo:
            description = gff_metainfo[description]
        gff_map[description] = arrange_gff_data(parse_gff(full_path))
    return gff_map


def get_primers_data(full_path, consensus):
    """Returns a map with filename key and primers json data."""
    # Hardcode metainformation for the filenames provided in the repository.
    primers_metainfo = {}
    primers_metainfo["nCoV-2019"] = "ARTIC bioinformatics primers for nanopore sequencing of nCoV2019"

    primers_map = {}
    consensus_upper = consensus.upper()
    description = os.path.splitext(os.path.basename(full_path))[0]
    if description in primers_metainfo:
        description = primers_metainfo[description]
    csv = pd.read_csv(full_path, sep=',')
    primers = [{'name': row[0], "seq":row[1].upper()}
               for row in csv[['name', 'seq']].values]
    primer_locations = []
    for entry in primers:
        lookup_sequence = entry['seq']
        # If the sequence corresponds to a `right` primer, then reverse and
        # complement.
        if "RIGHT" in entry['name'].upper():
            seq = Seq(lookup_sequence, IUPAC.unambiguous_dna)
            lookup_sequence = str(seq.reverse_complement())
        offsets = [
            m.start() for m in re.finditer(
                '(?=' + lookup_sequence + ')',
                consensus_upper)]
        for offset in offsets:
            primer_locations.append({'name': entry['name'],
                                     'seq': entry['seq'],
                                     'start': offset,
                                     'end': offset + len(entry['seq']) - 1})
    if primer_locations:
        primers_map[description] = arrange_gff_data(primer_locations)
    else:
        print("No primer was mapped from ", path, ".")
    return primers_map


def convert_coverage(fname, sample_name, genome_length):
    """Convert the read coverage to bp coverage."""
    csv = pd.read_csv(fname, sep="\t", index_col=0, header=0)
    return csv[sample_name].values.tolist()


def main(
    consensus_file,
    coverage_file,
    vcf_file,
    gff_directory,
    primers_file,
    html_file_in,
    html_file_out,
    wildcards_dataset,
    reference_file,
):
    # parse the sample name
    path_components = os.path.normpath(wildcards_dataset).split(os.path.sep)
    sample_name = '/'.join(path_components[1:])

    # parse the consensus sequence
    consensus = next(SeqIO.parse(consensus_file, "fasta")).seq.upper()

    # parse coverage file
    coverage = convert_coverage(
        coverage_file, sample_name.replace(
            '/', '-'), len(consensus))

    # load biodata in json format
    vcf_json = convert_vcf(vcf_file)
    gff_map = get_gff_data(gff_directory)
    primers_map = get_primers_data(primers_file, str(consensus))

    # parse the reference name
    reference_name = re.search(
        r"(.+)\.fa.*",
        os.path.basename(reference_file)
    ).group(1)

    embed_code = f"""
        var sample_name = \"{sample_name}\"
        var consensus = \"{consensus}\"
        var coverage = {coverage}
        var vcfData = {vcf_json}
        var gffData = {gff_map}
        var primerData = {primers_map}
        var reference_name = \"{reference_name}\"
    """

    # assemble webpage
    with open(html_file_in) as fd:
        raw_html = fd.read()

    # TODO: make this more robust
    mod_html = raw_html.replace("{EXTERNAL_SNAKEMAKE_CODE_MARKER}", embed_code)

    with open(html_file_out, "w") as fd:
        fd.write(mod_html)


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        sys.argv[5],
        sys.argv[6],
        sys.argv[7],
        sys.argv[8],
        sys.argv[9],
    )
