import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from pycirclize import Circos
from pycirclize.parser import Gff

def extract_features_by_strand(gff, sign, phase):
    filter_gff_records = []
    for rec in gff:
        if rec.strand != sign:
            continue
        if rec.phase != phase:
            continue
        filter_gff_records.append(rec)
    return [rec.to_seq_feature() for rec in filter_gff_records]

def color_track(f_cds_feats, cds_track, begin, end):
    for f in f_cds_feats:
        if ("vfdb_short_name" in f.qualifiers or "AMR_Gene_Family" in f.qualifiers):  # vfdb or CARD
            data_dict["vfdb_card"].append(f)
        else:  # no vfdb or card
            if f.qualifiers.get("phrog")[0] == "No_PHROG":
                data_dict["nophrog"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "unknown function":
                data_dict["unk"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "other":
                data_dict["other"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "tail":
                data_dict["tail"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "transcription regulation":
                data_dict["transcription"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "DNA":
                data_dict["dna"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "lysis":
                data_dict["lysis"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "moron":
                data_dict["moron"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "integration and excision":
                data_dict["int"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "head and packaging":
                data_dict["head"]["fwd_list"].append(f)
            elif f.qualifiers.get("function")[0] == "connector":
                data_dict["con"]["fwd_list"].append(f)

    for key in data_dict.keys():
        cds_track.genomic_features(
            data_dict[key]["fwd_list"],
            plotstyle="arrow",
            r_lim=(begin, end),
            fc=data_dict[key]["col"],
        )

def new_sign_on_track(gff, cds_track, begin, end, sign, strand):
    f_cds_feats = extract_features_by_strand(gff.records, sign, strand)
    color_track(f_cds_feats, cds_track, begin, end)

def overlapping_track(gff, color, sector, begin, end, strand):
    cds_track = sector.add_track((begin, end))
    cds_track.axis(fc=color, ec="none")

    middle=int((begin+end)/2)
    new_sign_on_track(gff, cds_track, begin, middle, 1, strand)
    new_sign_on_track(gff, cds_track, middle, end, -1, strand)
    return cds_track

# metadata phage
sample = sys.argv[1]
size_phage = int(sys.argv[2])

# gff genotate
gff_genotate = Gff(sys.argv[3])
circos = Circos(sectors = {sample: size_phage})
sector = circos.sectors[0]

# gff genotate
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
    
track_init = overlapping_track(gff_genotate, "#e3e3e3", sector, 70, 75, 0)
# Plot xticks & intervals on inner position
track_init.xticks_by_interval(
    interval=5000,
    outer=False,
    show_bottom_line=True,
    label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
    label_orientation="vertical",
    line_kws=dict(ec="grey"),
)
overlapping_track(gff_genotate, "#e3e3e3", sector, 76, 81, 1)
overlapping_track(gff_genotate, "#e3e3e3", sector, 82, 87, 2)

fig = circos.plotfig()

# Add legend
handle_phrogs = [
    Patch(color=data_dict["unk"]["col"], label="Unknown Function"),
    Patch(color=data_dict["other"]["col"], label="Other Function"),
    Patch(
        color=data_dict["transcription"]["col"], label="Transcription Regulation"
    ),
    Patch(
        color=data_dict["dna"]["col"], label="DNA/RNA & nucleotide \n metabolism"
    ),
    Patch(color=data_dict["lysis"]["col"], label="Lysis"),
    Patch(
        color=data_dict["moron"]["col"],
        label="Moron, auxiliary metabolic \n gene & host takeover",
    ),
    Patch(color=data_dict["int"]["col"], label="Integration & excision"),
    Patch(color=data_dict["head"]["col"], label="Head & packaging"),
    Patch(color=data_dict["con"]["col"], label="Connector"),
    Patch(color=data_dict["tail"]["col"], label="Tail"),
    Patch(color=data_dict["vfdb_card"]["col"], label="Virulence Factor/AMR"),
]

phrog_legend_coords = (0.10, 1.185)
phrog_legend = circos.ax.legend(
    handles=handle_phrogs,
    bbox_to_anchor=phrog_legend_coords,
    fontsize=9.5,
    loc="center",
    title="PHROG CDS",
    handlelength=2,
)

circos.ax.add_artist(phrog_legend)

# Save final figure
circos.savefig(sys.argv[4])
plt.close(fig)
