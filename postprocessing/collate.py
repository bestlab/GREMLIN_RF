import sqlite3, glob, os
import itertools, multiprocessing

f_db = "scores.db"
conn = sqlite3.connect(f_db)

template = '''
DROP TABLE IF EXISTS tmscore;
CREATE TABLE IF NOT EXISTS tmscore (
    idx INTEGER PRIMARY KEY AUTOINCREMENT,
    length_seq INTEGER,
    length_org INTEGER,
    length_common INTEGER,
    RMSD_common float,
    TMscore float,
    pdb CHAR(4),
    name STRING,
    f_tm STRING,
    UNIQUE (pdb, name, f_tm)
);
'''
conn.executescript(template)

F_TM = sorted(glob.glob("tmscores/*"))


def find_line(word, lines):
    for line in lines:
        if word in line:
            return line.split()
    return None

def read_TM(f_tm):
    with open(f_tm) as FIN:
        raw = FIN.read()
    lines = raw.split('\n')

    name = os.path.basename(f_tm).split('.tmscore')[0]
    pdb = name.split('_')[0]
    structure_idx = '.'.join(name.split('_')[-1].split('.')[1:])
    name = '_'.join(name.split('_')[1:-1])

    data = {
        "name" : name,
        "pdb": pdb,
        "f_tm": os.path.basename(f_tm),
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
        pass #print "problem with", f_tm

    return data


ITR = itertools.imap(read_TM, F_TM)

keys = ("f_tm",'length_seq','length_org','length_common',
        'RMSD_common','TMscore','pdb','name')

cmd_insert = '''
INSERT INTO tmscore ({}) VALUES ({})
'''.format(','.join(keys), ','.join(['?']*len(keys)))


P = multiprocessing.Pool()
ITR = P.imap(read_TM, F_TM)


for data in ITR:

    try:
        conn.execute(cmd_insert,[data[k] for k  in keys])
    except:
        print "Error with", data["f_tm"]

conn.commit()
