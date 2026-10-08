#!/usr/bin/env python3
"""Plot one Genotate/Pharokka comparison per FASTA contig."""
import argparse
from pathlib import Path
from urllib.parse import quote, unquote
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def read_features(path, lengths):
    features = defaultdict(list)
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected nine GFF columns")
            contig, _, kind, start, end, _, strand, phase, attrs = fields
            if kind != "CDS":
                continue
            start, end = int(start), int(end)
            if contig not in lengths or not 1 <= start <= end <= lengths[contig]:
                raise ValueError(f"{path}:{line_number}: unknown contig or invalid coordinates")
            qualifiers = {}
            for attribute in attrs.split(";"):
                if "=" in attribute:
                    key, value = attribute.split("=", 1)
                    qualifiers[unquote(key)] = [unquote(value)]
            features[contig].append(SeqFeature(
                FeatureLocation(start - 1, end, strand=1 if strand == "+" else -1),
                type="CDS", qualifiers=qualifiers,
            ))
    return features


def color_track(f_cds_feats, cds_track, begin, end):
    data_dict = {
        "nophrog": {"col": "#d3d3d3", "fwd_list": [], "rev_list": []},
        "vfdb_card": {"col": "#FF0000", "fwd_list": [], "rev_list": []},
        "unk": {"col": "#AAAAAA", "fwd_list": [], "rev_list": []},
        "other": {"col": "#4deeea", "fwd_list": [], "rev_list": []},
        "tail": {"col": "#74ee15", "fwd_list": [], "rev_list": []},
        "transcription": {"col": "#ffe700", "fwd_list": [], "rev_list": []},
        "dna": {"col": "#f000ff", "fwd_list": [], "rev_list": []},
        "lysis": {"col": "#001eff", "fwd_list": [], "rev_list": []},
        "moron": {"col": "#8900ff", "fwd_list": [], "rev_list": []},
        "int": {"col": "#E0B0FF", "fwd_list": [], "rev_list": []},
        "head": {"col": "#ff008d", "fwd_list": [], "rev_list": []},
        "con": {"col": "#5A5A5A", "fwd_list": [], "rev_list": []},
    }

    for f in f_cds_feats:
        if ("vfdb_short_name" in f.qualifiers or "AMR_Gene_Family" in f.qualifiers):  # vfdb or CARD
            data_dict["vfdb_card"]["fwd_list"].append(f)
        else:  # no vfdb or card
            if f.qualifiers.get("phrog", ["No_PHROG"])[0] == "No_PHROG":
                data_dict["nophrog"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "unknown function":
                data_dict["unk"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "other":
                data_dict["other"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "tail":
                data_dict["tail"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "transcription regulation":
                data_dict["transcription"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "DNA":
                data_dict["dna"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "lysis":
                data_dict["lysis"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "moron":
                data_dict["moron"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "integration and excision":
                data_dict["int"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "head and packaging":
                data_dict["head"]["fwd_list"].append(f)
            elif f.qualifiers.get("function", ["unknown function"])[0] == "connector":
                data_dict["con"]["fwd_list"].append(f)

            else:
                data_dict["other"]["fwd_list"].append(f)

    for key in data_dict.keys():
        cds_track.genomic_features(
            data_dict[key]["fwd_list"],
            plotstyle="arrow",
            r_lim=(begin, end),
            fc=data_dict[key]["col"],
        )

def new_sign_on_track(gff, cds_track, begin, end, sign, strand):
    f_cds_feats = [feature for feature in gff if feature.location.strand == sign
                   and (strand is None or (int(feature.location.start) % 3 if sign == 1
                        else (int(feature.location.end) - 1) % 3) == strand)]
    color_track(f_cds_feats, cds_track, begin, end)

def overlapping_track(gff, color, sector, begin, end, strand):
    cds_track = sector.add_track((begin, end))
    cds_track.axis(fc=color, ec="none")

    middle=int((begin+end)/2)
    new_sign_on_track(gff, cds_track, begin, middle, 1, strand)
    new_sign_on_track(gff, cds_track, middle, end, -1, strand)
    return cds_track

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--genotate-gff", required=True)
    parser.add_argument("--pharokka-gff", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    records = list(SeqIO.parse(args.fasta, "fasta"))
    lengths = {record.id: len(record.seq) for record in records}
    if len(lengths) != len(records) or any(length == 0 for length in lengths.values()):
        raise ValueError("Duplicate contig IDs or empty sequences in FASTA")
    genotate = read_features(args.genotate_gff, lengths)
    pharokka = read_features(args.pharokka_gff, lengths)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    # Remove stale plots when the set of contigs changes on a rerun.
    for stale in output.glob("*_annotated_by_genotate_pharokka.png"):
        stale.unlink()
    if not records:
        return
    from pycirclize import Circos
    for record in records:
        circos = Circos(sectors={record.id: lengths[record.id]})
        sector = circos.get_sector(record.id)
        sector.text(record.id, r=105, size=10)
        track = overlapping_track(genotate[record.id], "#e3e3e3", sector, 66, 70, 0)
        track.xticks_by_interval(
            interval=min(5000, max(1, lengths[record.id] // 5)), outer=False,
            show_bottom_line=True, label_formatter=lambda value: f"{value / 1000:.1f} Kb",
            label_orientation="vertical", line_kws=dict(ec="grey"),
        )
        overlapping_track(genotate[record.id], "#e3e3e3", sector, 71, 75, 1)
        overlapping_track(genotate[record.id], "#e3e3e3", sector, 76, 80, 2)
        overlapping_track(pharokka[record.id], "#e3e3e3", sector, 84, 90, None)
        filename = quote(record.id, safe="_-") + "_annotated_by_genotate_pharokka.png"
        fig = circos.plotfig()
        try:
            fig.savefig(output / filename, dpi=args.dpi)
        finally:
            plt.close(fig)


if __name__ == "__main__":
    main()
