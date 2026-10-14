"""Convert GenBank features to GFF3, retaining qualifiers and multipart locations."""

import argparse
from pathlib import Path
from urllib.parse import quote

from Bio import SeqIO


def escape(value):
    return quote(str(value), safe=".:^*$@!+_?-|")


def convert(input_path, output_path):
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(input_path) as source, open(output_path, "w") as output:
        output.write("##gff-version 3\n")
        for record_index, record in enumerate(SeqIO.parse(source, "genbank"), 1):
            seqid = escape(record.id)
            output.write(f"##sequence-region {seqid} 1 {len(record.seq)}\n")
            for feature_index, feature in enumerate(record.features, 1):
                if feature.location is None:
                    raise ValueError(f"Feature without a location in {record.id}")
                qualifiers = dict(feature.qualifiers)
                # One ID links all rows of a multipart feature without relying on
                # locus tags being unique across feature types or records.
                qualifiers.setdefault("ID", [f"phold_{record_index}_{feature_index}"])
                attributes = ";".join(
                    f"{escape(key)}={','.join(escape(value) for value in values)}"
                    for key, values in qualifiers.items()
                )
                coding_offset = int(qualifiers.get("codon_start", ["1"])[0]) - 1
                consumed = 0
                # Biopython orders joined parts in translation order, including
                # reverse-strand and origin-spanning CDS features.
                for part in feature.location.parts:
                    phase = str((coding_offset - consumed) % 3) if feature.type == "CDS" else "."
                    strand = {1: "+", -1: "-"}.get(part.strand, ".")
                    fields = [seqid, "phold", escape(feature.type), str(int(part.start) + 1),
                              str(int(part.end)), ".", strand, phase, attributes]
                    output.write("\t".join(fields) + "\n")
                    consumed += len(part)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    convert(args.input, args.output)
