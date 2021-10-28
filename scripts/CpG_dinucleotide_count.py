import sys

infile=sys.argv[1]
outfile=sys.argv[2]

new = open(outfile, "w")
k = open(infile)
for line in k:
    if ">" in line:
        temp = line.rstrip()+"\t"
    else:
        K = line.upper().count("CG")
        temp += str(round((K/(len(line)-1)), 4)) + "\n"
        new.write(temp)

new.close()
