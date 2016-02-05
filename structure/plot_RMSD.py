import glob, os, subprocess
import itertools, multiprocessing
import numpy as np
import pandas as pd

_PARALLEL = True
#_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))[:]


IGNORE_SET = set(["_0.50", "_1.50", "_3.00","_1.0",])
KEEP_SET = set(["exact","RF_5.0","GREMLIN_3.0"])

D2 = []
for x in D_SYSTEMS:

    FLAGS = [y for y in KEEP_SET if y in x]
    if not FLAGS: continue

    #FLAGS = [y for y in IGNORE_SET if y in x]
    #if FLAGS: continue
    
    D2.append(x)

D_SYSTEMS = D2

def run_system(dir):
    print "Loading", dir

    pdb = dir.split('/')[-1].split("_")[0]
    system = '_'.join(dir.split('/')[-1].split("_")[1:3])
    
    f_traj = os.path.join(dir,"RMSD.txt")
    Q = np.loadtxt(f_traj)

    # Keep the mean of the last 20 frames
    Q = Q[-200:]

    return pdb, system, Q.mean()

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

cols = ['pdb','system','Q']
data = list(ITR)
df = pd.DataFrame(data=data,columns=cols)

#import seaborn as sns
#sns.barplot(x='pdb',y='Q',hue='system',data=df)
#sns.plt.show()
#exit()

import seaborn as sns


df2 = pd.DataFrame(
    columns = df.system.unique(),
    index   = df.pdb.unique()
)

for k,row in df.iterrows():
    df2[row["system"]][row["pdb"]] = row.Q

print df2

def RMSD_diff(*args,**kwargs):
    x,y = args
    delta = (y-x).sum() / len(x)
    if delta > 0:
        font_weight = "bold"
        alpha = 1.0
    else:
        font_weight = "ultralight"
        alpha=0.6

    print delta
        
    sns.plt.text(0.05, 0.85, "{:0.3f}".format(delta),
                 fontsize=15, weight=font_weight,alpha=alpha)

def ref_line(*args,**kwargs):
    X = np.linspace(0,30,100)
    sns.plt.plot(X,X,'r--',alpha=0.8)
    sns.plt.xlim(0,30)
    sns.plt.ylim(0,30)
    #RMSD_diff(*args,**kwargs)

def zero_one_range(*args,**kwargs):

    avg = args[0].mean()
    
    sns.plt.xlim(0,30)
    
    font_weight = "bold"
    alpha = 1.0
    sns.plt.text(0.50, 22.50, "{:0.3f}".format(avg),
                 fontsize=15, weight=font_weight,alpha=alpha)

bins = np.linspace(0, 30, 25)
g = sns.pairplot(df2,diag_kws={"bins":bins})
g.map_lower(ref_line)
g.map_upper(ref_line)
g.map_diag(zero_one_range)

#import matplotlib
#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(3,3)

sns.plt.savefig("figures/pairplot_RMSD.png")

sns.plt.show()
