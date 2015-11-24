import sqlite3, glob, os
import itertools
#import multiprocessing

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
    rosetta_score float,
    pdb CHAR(4),
    name STRING,
    f_tm STRING,
    cluster_id INTEGER,
    cluster_model_id INTEGER,
    UNIQUE (pdb, name, f_tm)
);
'''
conn.executescript(template)



F_CLUSTERS = glob.glob("clusters/clustered_outfiles/*.out")
ROSETTA_SCORE = {}
for f in F_CLUSTERS:
    print "Loading", f
    with open(f) as FIN:
        for line in FIN:
            if "SCORE:" in line:

                # Skip the header lines
                if '.pdb' not in line:
                    continue
                
                f = os.path.basename(f)
                tokens = line.split()
                score = tokens[1]
                f_tm  = tokens[-1]
                pdb, model, params = f.split('_')
                params = params.replace('.out','')
                cluster_text = '.'.join(f_tm.split('.')[1:3])
                cluster_model_id,cluster_id  = map(int,cluster_text.split('.'))            
                name = model + "_" + params
                key = (pdb, name, cluster_id, cluster_model_id)
                ROSETTA_SCORE[key] = score

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

    cluster_text = '.'.join(os.path.basename(f_tm).split('.')[2:4])
    cluster_model_id,cluster_id  = map(int,cluster_text.split('.'))    
    
    key = (pdb, name, cluster_id, cluster_model_id)
    data["rosetta_score"] = float(ROSETTA_SCORE[key])

    data["cluster_id"] = cluster_id
    data["cluster_model_id"] = cluster_model_id
    
    return data


keys = ("f_tm",'length_seq','length_org','length_common',
        'RMSD_common','TMscore','pdb','name','rosetta_score','cluster_id','cluster_model_id')

cmd_insert = '''
INSERT INTO tmscore ({}) VALUES ({})
'''.format(','.join(keys), ','.join(['?']*len(keys)))


ITR = itertools.imap(read_TM, F_TM)
#P = multiprocessing.Pool()
#ITR = P.imap(read_TM, F_TM)


for data in ITR:
    #try:
    item = [data[k] for k  in keys]
    conn.execute(cmd_insert,item)
    #except:
    #print "Error with", data["f_tm"]

conn.commit()
