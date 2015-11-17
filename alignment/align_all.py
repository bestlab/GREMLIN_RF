import glob, os

F_FASTA = glob.glob("fasta/*.fasta")

os.system('mkdir -p aln')

bcmd = "./hhblits -i {f_fasta} -d uniprot20_2013_03/uniprot20_2013_03 -oa3m {f_aln} -n 4 -diff inf -cov 75 -e 0.0000000001 -cpu 20 -maxmem 12"

for f_fasta in F_FASTA:
    pdb = os.path.basename(f_fasta).split('.')[0]
    f_aln = "aln/{}.aln".format(pdb)

    if os.path.exists(f_aln):
        continue

    cmd = bcmd.format(f_aln=f_aln,f_fasta=f_fasta)
    print "starting", f_fasta
    os.system(cmd)
