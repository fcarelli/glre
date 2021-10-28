import sys

input_file = sys.argv[1]

output_inopinata = open(sys.argv[2], "w")
output_briggsae = open(sys.argv[3], "w")
output_nigoni = open(sys.argv[4], "w")
output_remanei = open(sys.argv[5], "w")

inopinata_ortho = []
briggsae_ortho = []
nigoni_ortho = []
remanei_ortho = []

infile = open(input_file)
header = infile.readline()

orthologs = infile.readline()

while orthologs != "":
    orthologs = orthologs[:-1]
    G = orthologs.split("\t")
    if G[3] == "ortholog_one2one":
        if G[1] not in inopinata_ortho:
            inopinata_ortho.append(G[1])
            output_inopinata.write(G[0] + "\t" + G[1] + "\t" + G[2] + "\n")
    if G[6] == "ortholog_one2one":
        if G[4] not in briggsae_ortho:
            briggsae_ortho.append(G[4])
            output_briggsae.write(G[0] + "\t" + G[4] + "\t" + G[5] + "\n")
    if G[9] == "ortholog_one2one":
        if G[7] not in nigoni_ortho:
            nigoni_ortho.append(G[7])
            output_nigoni.write(G[0] + "\t" + G[7] + "\t" + G[8] + "\n")
    if G[12] == "ortholog_one2one":
        if G[10] not in remanei_ortho:
            remanei_ortho.append(G[10])
            output_remanei.write(G[0] + "\t" + G[10] + "\t" + G[11] + "\n")
    orthologs = infile.readline()


output_briggsae.close()
output_nigoni.close()
output_remanei.close()
output_inopinata.close()
