import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
orthologs_file = sys.argv[3]

new = open(output_file, "w")

infile_ortho = open(orthologs_file)
infile_ortho_line = infile_ortho.readline()
ortholog=infile_ortho_line.split("\t")[1]

infile = open(input_file)
infile_line = infile.readline()

domains = {}
while infile_line != "":
    if infile_line[0] != "#":
        if len(infile_line.split()) == 12:
            K = infile_line.split()
            if K[0] == ortholog:
                domains[eval(K[3])] = [K[1], K[7], K[8]]
    infile_line = infile.readline()

print(domains)
if domains != {}:
    for i in range(len(domains)):
        temp = ortholog + "\t" + str(i+1) + "\t" + domains[i+1][0] + "\t" + domains[i+1][1] + "\t" + domains[i+1][2] + "\n"
        new.write(temp)

new.close()
