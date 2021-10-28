import sys

input_exp = sys.argv[1]
input_gene_names = sys.argv[2]

output_embryo_series_1_out = open(input_exp + ".embryo_series_1.txt", "w")
output_embryo_series_2_out = open(input_exp + ".embryo_series_2.txt", "w")
output_embryo_series_3_out = open(input_exp + ".embryo_series_3.txt", "w")
output_embryo_series_4_out = open(input_exp + ".embryo_series_4.txt", "w")
output_postembryo_stage_out = open(input_exp + ".postembryonic.txt", "w")

output_embryo_series_1_fields = []
output_embryo_series_2_fields = []
output_embryo_series_3_fields = []
output_embryo_series_4_fields = []
output_postembryo_stage_fields = []

postembryo_stages = ["L1_dcpm", "L2_dcpm", "L3_dcpm", "L4_dcpm", "YA_dcpm", "N2_Ad_gonad-1-RZLI_dcpm"]

input_exp_in = open(input_exp + ".txt")
header = input_exp_in.readline()

stages = header.rstrip().split()
for i in range(len(stages)):
    if "N2_EE_50" in stages[i]:
        if "dcpm" in stages[i]:
            output_embryo_series_1_fields.append(i)
    elif "20120223_EMB" in stages[i]:
        if "dcpm" in stages[i]:
            output_embryo_series_2_fields.append(i)
    elif "20120411_EMB" in stages[i]:
        if "dcpm" in stages[i]:
            output_embryo_series_3_fields.append(i)
    elif "20120419_EMB" in stages[i]:
        if "dcpm" in stages[i]:
            output_embryo_series_4_fields.append(i)
    elif stages[i] in postembryo_stages:
        output_postembryo_stage_fields.append(i)

gene_names = {}
for line in open(input_gene_names):
    G=line.rstrip().split(",")
    if G[3] != "":
        if G[1] != "":
            gene_names[G[3]] = G[1]


output_embryo_series_1_out.write("\t".join([stages[i] for i in output_embryo_series_1_fields]) + "\n")
output_embryo_series_2_out.write("\t".join([stages[i] for i in output_embryo_series_2_fields]) + "\n")
output_embryo_series_3_out.write("\t".join([stages[i] for i in output_embryo_series_3_fields]) + "\n")
output_embryo_series_4_out.write("\t".join([stages[i] for i in output_embryo_series_4_fields]) + "\n")
output_postembryo_stage_out.write("\t".join([stages[i] for i in output_postembryo_stage_fields]) + "\n")

gene_exp = input_exp_in.readline()
while gene_exp:
    gene_exp_split = gene_exp.rstrip().split()
    if gene_exp_split[0] in gene_names:
        output_embryo_series_1_out.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join([gene_exp_split[i] for i in output_embryo_series_1_fields]) + "\n")
        output_embryo_series_2_out.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join([gene_exp_split[i] for i in output_embryo_series_2_fields]) + "\n")
        output_embryo_series_3_out.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join([gene_exp_split[i] for i in output_embryo_series_3_fields]) + "\n")
        output_embryo_series_4_out.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join([gene_exp_split[i] for i in output_embryo_series_4_fields]) + "\n")
        output_postembryo_stage_out.write(gene_names[gene_exp_split[0]] + "\t" + "\t".join([gene_exp_split[i] for i in output_postembryo_stage_fields]) + "\n")
    gene_exp = input_exp_in.readline()

output_embryo_series_1_out.close()
output_embryo_series_2_out.close()
output_embryo_series_3_out.close()
output_embryo_series_4_out.close()
output_postembryo_stage_out.close()
