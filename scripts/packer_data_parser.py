import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

new = open(output_file, "w")

header = "bin1_rawTPM\tbin1_adjTPM\tbin1_bootstrapTPM\tbin2_rawTPM\tbin2_adjTPM\tbin2_bootstrapTPM\tbin3_rawTPM\tbin3_adjTPM\tbin3_bootstrapTPM\n"
new.write(header)

for line in open(input_file):
    G = line.rstrip().split("\t")
    if G[1] == "gene.id":
        gene = ""
    else:
        if gene == "":
            gene = G[1]
            temp = G[1] + "\t"
        else:
            if G[1] == gene:
                if "Germline" in G[2]:
                    temp += G[3] + "\t" + G[4] + "\t" + G[5] + "\t"
            else:
                new.write(temp.rstrip() + "\n")
                gene = G[1]
                temp = G[1] + "\t"

new.write(temp.rstrip() + "\n")
new.close()
