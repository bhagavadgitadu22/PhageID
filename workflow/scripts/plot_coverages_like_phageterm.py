import sys
import pysam
from collections import defaultdict

bam_path = sys.argv[1]

# Open the BAM file
bamfile = pysam.AlignmentFile(bam_path, "rb")

# Initialize dictionaries
coverage = defaultdict(lambda: defaultdict(int))
start_positions = defaultdict(lambda: defaultdict(int))
end_positions = defaultdict(lambda: defaultdict(int))

# Iterate through reads
for read in bamfile.fetch():
    if read.is_unmapped:
        continue

    chrom = bamfile.get_reference_name(read.reference_id)

    # Coverage per position
    for pos in range(read.reference_start, read.reference_end):
        coverage[chrom][pos] += 1

    # Start and end positions (adjust for strand)
    if read.is_reverse:
        end_positions[chrom][read.reference_end - 1] += 1
    else:
        start_positions[chrom][read.reference_start] += 1

bamfile.close()

# Write BEDGRAPH files
def write_bedgraph(filename, data_dict):
    with open(filename, "w") as out:
        for chrom in sorted(data_dict):
            for pos in sorted(data_dict[chrom]):
                value = data_dict[chrom][pos]
                # bedGraph is zero-based, half-open [start, end)
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{value}\n")

# Write the three bedGraph files
write_bedgraph(sys.argv[2], coverage)
write_bedgraph(sys.argv[3], start_positions)
write_bedgraph(sys.argv[4], end_positions)
