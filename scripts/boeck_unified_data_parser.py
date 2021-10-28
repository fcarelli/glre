import sys

input_gene_names = sys.argv[1]

input_exp_in = open(sys.argv[2])
output_embryo_stage = open(sys.argv[3], "w")

gene_names = {}
for line in open(input_gene_names):
    G=line.rstrip().split(",")
    if G[3] != "":
        if G[1] != "":
            gene_names[G[3]] = G[1]



gene_exp = input_exp_in.readline()
temp = "\t".join(gene_exp.rstrip().split()[1:19]) + "\n"
output_embryo_stage.write(temp)

gene_exp = input_exp_in.readline()
while gene_exp:
    gene_exp_split = gene_exp.rstrip().split()
    if gene_exp_split[0] in gene_names:
        output_embryo_stage.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join(gene_exp_split[1:19]) + "\n")
    gene_exp = input_exp_in.readline()

output_embryo_stage.close()

