import sys

input_f = sys.argv[1]
output_f = sys.argv[2]

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

new = open(output_f, "w")
prot=""
for line in open(input_f):
    if line[0] == ">":
        if prot == "":
            new.write(line)
            continue
        else:
            for chunk in chunkstring(prot, 60):
                new.write(chunk+"\n")
            prot = ""
            new.write(line)
            continue
    else:
        prot += line.rstrip()

for chunk in chunkstring(prot, 60):
    new.write(chunk+"\n")

new.close()

            
