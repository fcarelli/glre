import sys

input_file = sys.argv[1]

output_ceylanicum = open(sys.argv[2], "w")
output_contortus = open(sys.argv[3], "w")
output_bacteriophora = open(sys.argv[4], "w")
output_tipulae = open(sys.argv[5], "w")

ceylanicum_ortho = []
contortus_ortho = []
bacteriophora_ortho = []
tipulae_ortho = []

infile = open(input_file)
header = infile.readline()

orthologs = infile.readline()

while orthologs != "":
    orthologs = orthologs[:-1]
    G = orthologs.split("\t")
    if G[3] == "ortholog_one2one":
        if G[1] not in ceylanicum_ortho:
            ceylanicum_ortho.append(G[1])
            output_ceylanicum.write(G[0] + "\t" + G[1] + "\t" + G[2] + "\n")
    if G[6] == "ortholog_one2one":
        if G[4] not in contortus_ortho:
            contortus_ortho.append(G[4])
            output_contortus.write(G[0] + "\t" + G[4] + "\t" + G[5] + "\n")
    if G[9] == "ortholog_one2one":
        if G[7] not in bacteriophora_ortho:
            bacteriophora_ortho.append(G[7])
            output_bacteriophora.write(G[0] + "\t" + G[7] + "\t" + G[8] + "\n")
    if G[12] == "ortholog_one2one":
        if G[10] not in tipulae_ortho:
            tipulae_ortho.append(G[10])
            output_tipulae.write(G[0] + "\t" + G[10] + "\t" + G[11] + "\n")
    orthologs = infile.readline()


output_ceylanicum.close()
output_contortus.close()
output_bacteriophora.close()
output_tipulae.close()
