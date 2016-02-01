import sqlite3, glob, os
import pandas as pd
from sqlalchemy import create_engine # database connection

f_rosetta = "master_rosetta_scores.db"
assert(os.path.exists(f_rosetta))

print "Create ROSETTA score database index"
cmd = '''
CREATE INDEX IF NOT EXISTS item_lookup ON
data (pdb, name, group_id, structure_id)
'''
#conn = sqlite3.connect(f_rosetta)
#conn.execute(cmd)
#conn.commit()
#conn.close()

print "Loading ROSETTA score database"
engine = create_engine('sqlite:///{}'.format(f_rosetta))

def TM_file_iterator():
    for f_tm in glob.glob("tmscores/*"):
        yield f_tm

def find_line(word, lines):
    for line in lines:
        if word in line:
            return line.split()
    return None

def read_TM(f_tm):
    print "Reading", f_tm
    
    with open(f_tm) as FIN:
        raw = FIN.read()
    lines = raw.split('\n')

    f = os.path.basename(f_tm)
    tokens = f.split('.tmscore')[0].split('_')

    data = {
        "pdb" : tokens[0],
        "name" : tokens[1] + '_' + tokens[2],
        "group_id" : tokens[3],
        "structure_id" : int(tokens[5]),
    }

    try:
        data.update( {
            "length_seq" : int(find_line("Structure1:",lines)[-1]),
            "length_org" : int(find_line("Structure2:",lines)[-7]),
            "length_common" : int(find_line("residues in common=",lines)[-1]),
            "RMSD_common"   : float(find_line("RMSD of",lines)[-1]),
            "TMscore" : float(find_line("TM-score  ",lines)[2]),
        })
    except:
        pass
        print "problem with", f_tm


    cmd = '''
    SELECT * FROM data WHERE
    pdb='{pdb}' AND
    structure_id={structure_id} AND
    group_id={group_id} AND
    name='{name}'
    LIMIT 2
    '''.format(**data)

    rdata = pd.read_sql_query(cmd, engine)
    assert( len(rdata) == 1 )

    rdata = dict(rdata.ix[0])
    rdata.pop('index')
    data.update(rdata)

    return data


import itertools
ITR = itertools.imap(read_TM, TM_file_iterator())
data = list(ITR)
df = pd.DataFrame(data)


print "Sorting"
df.sort_values(["pdb","name"],inplace=True)

print "Saving"
f_sqlite_out = "master_tmscores.db"
engine = create_engine('sqlite:///{}'.format(f_sqlite_out))
df.to_sql('data', engine, if_exists='replace')

