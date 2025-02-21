import csv
import re
import sys

gffs = {}
with open(sys.argv[3], newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
    for row in spamreader:
        name = row[0]
        length = row[1]
        phrog = row[2]
        annot = row[3].replace("DNA, RNA and nucleotide metabolism", "DNA").replace("moron, auxiliary metabolic gene and host takeover", "moron")
        category = row[4].replace("DNA, RNA and nucleotide metabolism", "DNA").replace("moron, auxiliary metabolic gene and host takeover", "moron")

        pattern = r"(.*?)_CDS_\[(complement\()?(\d+)\.\.(\d+)(\))?\]$"
        match = re.match(pattern, name)

        if match:
            contig = match.group(1)
            start = match.group(3)
            end = match.group(4)

            strand = "-"
            if match.group(2) is None:
                strand = "+"
            phased = str((int(start)+2)%3)

            gff = [contig, "genotate_0.15", "CDS", start, end, ".", strand, phased, "ID="+contig+"_"+str(start)+"_"+str(end)+";phrog="+str(phrog)+";function="+category+";product="+category]
            str_gff = '\t'.join(gff)
            if contig not in gffs:
                gffs[contig] = []
            gffs[contig].append([str_gff])
        else:
            print(row)

for contig in gffs:
    print(contig)
    contig_renamed = contig.replace("/", "_")
    with open(sys.argv[5], 'w', newline='') as gff_file:
        writer = csv.writer(gff_file)
        writer.writerow(["##gff-version 3"])
        writer.writerow(["##sequence "+str(sys.argv[1])+" "+str(sys.argv[2])])
        writer.writerows(gffs[contig])
        writer.writerow(["##FASTA"])
        with open(sys.argv[4], 'r') as source:
            gff_file.write(source.read())
