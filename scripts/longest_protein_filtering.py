import sys
import builtins

input_pep=sys.argv[1]
input_long_isoform=sys.argv[2]
out_pep=sys.argv[3]
in_gff=sys.argv[4]
species=sys.argv[5]

longest_isoform=[]
for line in open(input_long_isoform):
    longest_isoform.append(line.split("\t")[0])

if species == "remanei":
    longest_isoform_dict = {}
    for line in open(input_long_isoform):
        longest_isoform_dict[line.split("\t")[0]] = line.rstrip().split("\t")[1]
    longest_isoform_prot = {}
    for line in open(in_gff):
        if line[0] != "#":
            G=line.rstrip().split("\t")
            if G[2] == "CDS":
                gene_detail=G[8].split(";")
                if gene_detail[1].split(":")[1] in longest_isoform:
                    if gene_detail[0].split("-")[1] not in longest_isoform_prot:
                        longest_isoform_prot[gene_detail[0].split("-")[1]] = longest_isoform_dict[gene_detail[1].split(":")[1]]

new=open(out_pep, "w")
flag=0
for line in open(input_pep):
    if line[0] == ">":
        flag=0
        if species == "remanei":
            if line.split(":")[0][1:] in longest_isoform_prot:
                flag = 1
                new.write(">" + longest_isoform_prot[line.rstrip().split(":")[1]] + "\n")
        else:
            if line.split(":")[0][1:] in longest_isoform:
                flag = 1
                new.write(">" + line.split(":")[1])
                longest_isoform.remove(line.split(":")[0][1:])
            elif builtins.any(line.split(":")[0][1:] in transcript for transcript in longest_isoform):
                flag = 1
                result = [transcript for transcript in longest_isoform if line.split(":")[0][1:] in transcript]
                if (len(result) > 1 ):
                    for transcript in result:
                        if transcript.rsplit(".", 1)[0] == line.split(":")[0][1:]:
                            longest_isoform.remove(transcript)
                            new.write(">" + line.split(":")[1])
                else:
                    longest_isoform.remove(result[0])
                    new.write(">" + line.split(":")[1])
    else:
        if flag == 1:
            new.write(line)

new.close()
