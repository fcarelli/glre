import sys

infile = sys.argv[1]
outfile = sys.argv[2]
dinucl = sys.argv[3]

new = open(outfile, "w")

len_seq = set()
for line in open(infile):
    if line[0] != ">":
        len_seq.add(len(line.rstrip()))

for line in open(infile):
    if line[0] == ">":
        chrom=line.split(":")[0]
        start_coord=line.split(":")[1].split("-")[0]
        end_coord=line.split(":")[1].split("-")[1].split("(")[0]
        temp = chrom + "_" + start_coord + "_" + end_coord + "\t"
    else:
        for i in range(max(len_seq)):
            if dinucl in line[i:i+2]:
                temp += "1\t"
            else:
                temp +=	"0\t"
        temp = temp.rstrip() + "\n"
        new.write(temp)

new.close()
