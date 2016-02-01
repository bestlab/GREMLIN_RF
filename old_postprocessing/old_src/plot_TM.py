import sqlite3, glob, os
import itertools, multiprocessing
import numpy as np
import pandas as pd

pdb_filter = "1vfy"

f_db = "scores.db"
conn = sqlite3.connect(f_db)

cmd_grab = '''
SELECT * FROM tmscore WHERE pdb=?
'''

conn.execute(cmd_grab,(pdb_filter,))
df = pd.read_sql_query(cmd_grab, conn)

import seaborn as sns

sns.plt.figure(figsize=(12,8))
sns.plt.title("Protein {} for different models".format(pdb_filter))

for g,data in df.groupby(["name"]):
    print g
    X = data.TMscore[:20000]

    if "exact" in g:
        kde_kws={"color":"k", "ls":'-', "lw":3, "alpha":0.5}
        continue

    if "RF" in g:
        kde_kws={"ls":'-', "lw":2, "alpha":0.95}

    if "GREMLIN" in g:
        kde_kws={"ls":'--', "lw":2, "alpha":0.95}
    
    sns.distplot(X,label=g,hist=False,kde_kws=kde_kws)
    
sns.plt.xlim(0.2,1)
sns.plt.legend()
sns.plt.savefig("figures/tmplot_{}.png".format(pdb_filter))
sns.plt.show()

