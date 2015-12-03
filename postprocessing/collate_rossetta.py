import sqlite3, glob, os
import pandas as pd
from sqlalchemy import create_engine # database connection

def process_model(f):
    print "Collating rosetta for", f

    name = os.path.basename(f).split('score_')[-1]
    name = '_'.join(name.split('_')[:2])
    
    FIN = open(f)
    
    pdb = os.path.abspath(f).split('/')[-2]
    group_id = int(os.path.basename(f).split('_')[-1].split('.')[0])

    rows = ["pdb","name","group_id","structure_id"]
    rows = rows + FIN.readline().split()[1:-1]
    
    raw_data = []
    for line in FIN:
        line = line.split()[1:]
        structure_id = int(line[-1].split('_')[-1])
        line = [pdb, name, group_id, structure_id] + map(float,line[:-1])
        raw_data.append(line)
        
    df = pd.DataFrame(raw_data, columns=rows)
    FIN.close()
    
    return df

def file_iterator():
    FOLD_DIRS = sorted(glob.glob("rosetta_folding_data/*"))
    for d in FOLD_DIRS:
        F_SCORE = sorted(glob.glob(os.path.join(d,"*.fsc")))
        for k,f in enumerate(F_SCORE):
            yield f

import itertools, multiprocessing

ITR = itertools.imap(process_model, file_iterator())

#P = multiprocessing.Pool()
#ITR = P.imap(process_model, file_iterator())

MASTER_DF = list(ITR)

print "Concatenating"
df = pd.concat(MASTER_DF)

print "Sorting"
df.sort_values(["pdb","name","group_id","structure_id"],inplace=True)
print df.shape

print "Saving"
f_sqlite_out = "master_rosetta_scores.db"
engine = create_engine('sqlite:///{}'.format(f_sqlite_out))

df.to_sql('data', engine, if_exists='replace',chunksize=10**4)

