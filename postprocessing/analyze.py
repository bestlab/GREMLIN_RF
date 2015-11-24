import sqlite3, glob, os
import itertools, multiprocessing
import numpy as np
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
        
        #item = [pdb + '_' + name, data1["TMscore"].max()]
        TM = data1.groupby("cluster_id").mean()["TMscore"].values
        R = data1.groupby("cluster_id").mean()["rosetta_score"].values

        #TM = data1["TMscore"].values
        #R  = data1["rosetta_score"].values
        idx = np.argsort(R)

        idx = idx[:1]
        print TM[idx].mean()
        #score = TM[idx].mean()
        #score = TM.max()
        #score = TM[idx]

        score = TM[idx].mean()
        datum = [name,pdb,score]
        ITEMS.append(datum)
        

        #for i in idx:
        #    score = TM[i]
        #    datum = [name,pdb,score]
        #    ITEMS.append(datum)

df2 = pd.DataFrame(columns=("name","pdb","tmscore"),data=ITEMS)
names = df2.name.unique()
df3 = pd.DataFrame(columns=names, index=df2.pdb.unique())

for name,pdb,tmscore in ITEMS:
    df3[name][pdb] = tmscore

def ref_line(*args,**kwargs):
    X = np.linspace(0,1,100)
    sns.plt.plot(X,X,'r--',alpha=0.8)
    sns.plt.xlim(0,1)
    sns.plt.ylim(0,1)

def zero_one_range(*args,**kwargs):
    sns.plt.xlim(0,1)
    
import seaborn as sns
g = sns.pairplot(df3)
g.map_upper(ref_line)
g.map_lower(ref_line)
g.map_diag(zero_one_range)

sns.plt.savefig("figures/pairplot.png")

sns.plt.figure(figsize=(16,7))
sns.barplot(y="tmscore",x="pdb",hue="name",data=df2)
sns.plt.tight_layout()
sns.plt.savefig("figures/barchart.png")

sns.plt.show()
