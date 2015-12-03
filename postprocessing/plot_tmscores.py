import os
import pandas as pd
import numpy as np
from sqlalchemy import create_engine # database connection

f_db = "master_tmscores.db"
assert(os.path.exists(f_db))

engine = create_engine('sqlite:///{}'.format(f_db))
df = pd.read_sql('data',engine)

#df = df[["name","pdb","score"]]
df2 = pd.DataFrame(
    columns=df.name.unique(),
    index  = df.pdb.unique()
)

for k,row in df.iterrows():
    df2[row["name"]][row["pdb"]] = row.TMscore


def ref_line(*args,**kwargs):
    X = np.linspace(0,1,100)
    sns.plt.plot(X,X,'r--',alpha=0.8)
    sns.plt.xlim(0,1)
    sns.plt.ylim(0,1)

    RMSD_diff(*args,**kwargs)

def RMSD_diff(*args,**kwargs):
    x,y = args
    delta = (y-x).sum() / len(x)
    if delta > 0:
        font_weight = "bold"
        alpha = 1.0
    else:
        font_weight = "ultralight"
        alpha=0.6
        
        
    sns.plt.text(0.05, 0.85, "{:0.3f}".format(delta),
                 fontsize=15, weight=font_weight,alpha=alpha)
    
def zero_one_range(*args,**kwargs):
    sns.plt.xlim(0,1)
    
import seaborn as sns
g = sns.pairplot(df2)
g.map_lower(ref_line)
g.map_upper(ref_line)
g.map_diag(zero_one_range)
sns.plt.savefig("figures/pairplot.png")


sns.plt.figure(figsize=(16,7))
sns.barplot(y="TMscore",x="pdb",hue="name",data=df)
sns.plt.legend(loc=0)
sns.plt.tight_layout()
sns.plt.savefig("figures/barchart.png")

sns.plt.show()
