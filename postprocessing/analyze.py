import sqlite3, glob, os
import itertools, multiprocessing
import pandas as pd

f_db = "scores.db"
conn = sqlite3.connect(f_db)

cmd_grab = '''
SELECT * FROM tmscore
'''

conn.execute(cmd_grab)
df = pd.read_sql_query(cmd_grab, conn)

print "Looking for mismatched items"
for g,data in df.groupby(["pdb"]):
    common = data["length_common"].values[0]
    lseq   = data["length_seq"].values[0]

    if lseq - common:
        print "MISMATCH FROM PDB FILES!"
        print g, lseq - common
        exit()
        
ITEMS = []

for pdb,data0 in df.groupby("pdb"):
    for name,data1 in data0.groupby("name"):
        
        item = [pdb + '_' + name, data1["TMscore"].max()]
        
        print name, pdb, len(data1)
        
        ITEMS.append([name,pdb,data1["TMscore"].max()])

exit()

df2 = pd.DataFrame(columns=("name","pdb","tmscore"),data=ITEMS)

import seaborn as sns
sns.plt.figure(figsize=(12,6))
sns.barplot(y="tmscore",x="pdb",hue="name",data=df2)
sns.plt.tight_layout()
sns.plt.show()




