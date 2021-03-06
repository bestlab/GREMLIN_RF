# Run extract_pdb.py first

import glob, os, itertools, multiprocessing
from glob import glob

os.system('mkdir -p tmscores')

def PDB_ITR():
    PDB_DIRS = sorted(os.listdir('pdb'))
    for pdb_dir in PDB_DIRS:
        os.system('mkdir -p tmscores/{}'.format(pdb_dir))
        FILES = glob(os.path.join("pdb", pdb_dir, "*.pdb"))
        for f_pdb in sorted(FILES):
            yield f_pdb


def process(f_pdb):

    name = f_pdb.split('/')[1]
    pdb  = name.split('_')[0]

    f_out = os.path.basename(f_pdb).replace('.pdb','.tmscore')    
    f_out = os.path.join("tmscores",name,f_out)
    
    if os.path.exists(f_out):
        return None

    # Find "gold standard pdb" and check it exists    
    f_org_pdb = os.path.join('pdb_standard',pdb+'.pdb')
    assert(os.path.exists(f_org_pdb))

    cmd = './TMscore {} {} > {}'.format(f_pdb, f_org_pdb, f_out)
    os.system(cmd)

    return f_out

ITR = itertools.imap(process, PDB_ITR())
P = multiprocessing.Pool()
ITR = P.imap(process, PDB_ITR(), chunksize=20)

for result in ITR:
    if result is not None:
        print result
