import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

hits_dist={}
for line in open(input_file):
    G=line.rstrip().split("\t")
    motif_B=G[6]+"\t"+G[7]+"\t"+G[8]
    if motif_B not in hits_dist:
        hits_dist[motif_B] = []
    hits_dist[motif_B].append(eval(G[12]))

new = open(output_file, "w")
for line in open(input_file):
    G=line.rstrip().split("\t")
    motif_B=G[6]+"\t"+G[7]+"\t"+G[8]
    if eval(G[12]) == min(hits_dist[motif_B]):
        temp = ""
        for i in range(0, len(G)-1):
            temp += G[i] + "\t"
        if eval(G[1]) < eval(G[7]):
            temp += str(eval(G[7])-eval(G[2])-1) + "\n"
        else:
            temp += str(eval(G[1])-eval(G[8])-1) + "\n"
        new.write(temp)

new.close()

