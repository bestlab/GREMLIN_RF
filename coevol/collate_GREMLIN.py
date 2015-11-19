import glob, os, shutil

F_ALN = sorted(glob.glob("input_GREMLIN/*.aln"))
os.system('mkdir -p APC')

count = 0

for f_aln in F_ALN:
    pdb = os.path.basename(f_aln).split('.')[0]
    f_gremlin = os.path.join("output_GREMLIN",pdb,pdb+'.gremlin')
    if not os.path.exists(f_gremlin):
        print "Can't find {}.".format(f_gremlin)
        continue

    f_APC = os.path.join("APC", os.path.basename(f_gremlin))
    shutil.copyfile(f_gremlin, f_APC)
    count += 1

print "Found {} GREMLIN files".format(count)
