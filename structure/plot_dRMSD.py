import glob, os, subprocess
import itertools, multiprocessing
import numpy as np
import pandas as pd

_PARALLEL = True
#_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))[:]

IGNORE_SET = set(["_0.50", "_1.50", "_3.00","_1.0",])
KEEP_SET = set(["exact","RF_5.0","GREMLIN_5.0"])

cutoff_N_frames = 50

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
    
    f_traj = os.path.join(dir,"dRMSD.txt")
    RMSD = np.loadtxt(f_traj)

    # Keep the mean of the last 200 frames
    RMSD = RMSD[-cutoff_N_frames:]

    return pdb, system, RMSD.mean()

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

cols = ['pdb','system','RMSD']
data = list(ITR)
df = pd.DataFrame(data=data,columns=cols)

import seaborn as sns

df2 = pd.DataFrame(
    columns = df.system.unique(),
    index   = df.pdb.unique()
)

for k,row in df.iterrows():
    df2[row["system"]][row["pdb"]] = row.RMSD

print df2


def ref_line(*args,**kwargs):
    X = np.linspace(0,30,100)
    sns.plt.plot(X,X,'r--',alpha=0.8)
    #sns.plt.xlim(0,max_view_RMSD)
    sns.plt.ylim(0,14)
    sns.plt.xlim(0,14)

def zero_one_range(*args,**kwargs):

    avg = args[0].mean()
    
    #sns.plt.xlim(0,max_view_RMSD)
    #sns.plt.xlim(0,12)
    
    font_weight = "bold"
    alpha = 1.0
    sns.plt.text(0.50, 22.50, "{:0.3f}".format(avg),
                 fontsize=15, weight=font_weight,alpha=alpha)


max_view_RMSD = int(df.RMSD.max()) + 1

bins = np.linspace(0, max_view_RMSD, 40)
g = sns.pairplot(df2,diag_kws={"bins":bins})
g.map_lower(ref_line)
g.map_upper(ref_line)
g.map_diag(zero_one_range)

#import matplotlib
#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(3,3)

os.system('mkdir -p figures/')
sns.plt.savefig("figures/pairplot_RMSD.png")

sns.plt.show()
