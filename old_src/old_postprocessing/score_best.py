import sqlite3, glob, os
import pandas as pd
from sqlalchemy import create_engine # database connection

#FLAG_DEBUG = True
FLAG_DEBUG = False

os.system('mkdir -p pdb tmscores')

f_sqlite_out = "master_rosetta_scores.db"
assert(os.path.exists(f_sqlite_out))

def get_silent_file(pdb,name,group_id):
    
    f = "*{}_{}.out".format(name, group_id)
    f = os.path.join("rosetta_folding_data",pdb,f)
    files = glob.glob(f)
    assert(len(files)==1)
    return files.pop()

def extract_pdb(pdb,name,group_id,structure_id):
    f_silent = get_silent_file(pdb, name, group_id)
    tag = "S_{:08d}".format(structure_id)
    prefix = os.path.join("pdb","{}_{}_{}_".format(pdb,name,group_id))
        
    f_pdb = prefix + tag + '.pdb'
    if os.path.exists(f_pdb):
        return f_pdb

    print "Extract PDB", f_pdb
    
    cmd = ("./extract_pdbs.linuxgccrelease "
           "-in:file:silent {} "
           "-in:file:tags {} "
           "-out:prefix {}")
    cmd = cmd.format(f_silent, tag, prefix)
    os.system(cmd)

    return f_pdb


def score_structure(f_pdb):

    name = f_pdb.split('/')[1]
    pdb  = name.split('_')[0]

    f_TM = os.path.join("tmscores",name.replace('.pdb','.tmscore'))
    
    if os.path.exists(f_TM):
        return f_TM

    print "Scoring", f_TM

    # Find "gold standard pdb" and check it exists    
    f_org_pdb = os.path.join('pdb_standard',pdb+'.pdb')
    assert(os.path.exists(f_org_pdb))

    cmd = './TMscore {} {} > {}'.format(f_pdb, f_org_pdb, f_TM)
    os.system(cmd)

    return f_TM

def process_model(data):
    (pdb,name),data = data

    min_val = data["score"].min()

    X = data["score"] #- data["atom_pair_constraint"]
    min_val = X.min()
    
    p = data[X==min_val].values[0]
    p = dict(zip(cols,p))

    f_pdb = extract_pdb(pdb, name, p["group_id"],p["structure_id"])
    f_tm  = score_structure(f_pdb)
    
    return f_tm
    

print "Loading ROSETTA score database"
engine = create_engine('sqlite:///{}'.format(f_sqlite_out))

if FLAG_DEBUG:
    df = pd.read_sql_query('SELECT * FROM data LIMIT 1000',engine)
else:
    df = pd.read_sql_query('SELECT * FROM data',engine)

# Get column names
cols = pd.read_sql_query('SELECT * FROM data LIMIT 0',engine)
cols = cols.columns.values

INPUT_ITR = df.groupby(["pdb","name"])

import itertools, multiprocessing

if FLAG_DEBUG:
    ITR = itertools.imap(process_model, INPUT_ITR)
else:
    P = multiprocessing.Pool()
    ITR = P.imap(process_model, INPUT_ITR)

for f_TM in ITR:
    print f_TM


