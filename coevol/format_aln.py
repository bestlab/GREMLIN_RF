import os, glob

os.system('mkdir -p input_GREMLIN')

F_ALN = glob.glob("aln/*")

for f in F_ALN:
    with open(f) as FIN:
        data = [line for line in FIN if line[0] != ">" and line.strip()]
    print "Loaded {} lines from {}".format(len(data), f)
    f_out = os.path.join("input_GREMLIN", os.path.basename(f))

    # Remove lower case letters in alignment
    seq = [a for a in data if a.isupper()]
    seq = ''.join(seq)

    with open(f_out,'w') as FOUT:
        FOUT.write(seq)
