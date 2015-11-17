import glob, subprocess, os, shutil

def mkdir_p(d):
    try:
        os.mkdir(d)
    except OSError:
        pass

mkdir_p("pdb")
mkdir_p("fasta")
F_PDB = glob.glob("org_pdb/????.pdb")

for f_pdb in F_PDB:
    chain_id = "A"
    
    f_pdb = os.path.abspath(f_pdb)
    pdb = os.path.basename(f_pdb).split('.')[0]

    cmd = "python rosetta_scripts/clean_pdb.py {} {}"
    cmd = cmd.format(f_pdb,chain_id)
    
    os.system(cmd)
    shutil.move(pdb+"_A.pdb","pdb/"+pdb+".pdb")
    shutil.move(pdb+"_A.fasta","fasta/"+pdb+".fasta")

    print "Cleaned", f_pdb
