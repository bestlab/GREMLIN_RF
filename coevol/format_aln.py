import os, glob

os.system('mkdir -p input_GREMLIN')

F_ALN = glob.glob("aln/*")

for f in F_ALN:
    with open(f) as FIN:
        data = [line for line in FIN if line[0] != ">" and line.strip()]
    print "Loaded {} lines from {}".format(len(data), f)
    f_out = os.path.join("input_GREMLIN", os.path.basename(f))

    # Remove lower case letters in alignment
    alignments = [''.join([x for x in seq.strip() if x.isupper() or x=="-"])
                  for seq in data]

    # Sanity check, make sure all alignments are the same size
    sizes = map(len,alignments)
    assert(len(set(sizes)) == 1)
    
    alignments = '\n'.join(alignments)
    
    with open(f_out,'w') as FOUT:
        FOUT.write(alignments)
