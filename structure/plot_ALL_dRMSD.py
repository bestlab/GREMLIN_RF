import glob, os, subprocess
import itertools, multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
plt = sns.plt

_PARALLEL = True
MP_CORES = 30
#D_SYSTEMS = sorted(glob.glob("systems/1a*"))[:]
D_SYSTEMS = sorted(glob.glob("systems/*"))[:]
cut_frames = 50

def load_residues(pdb):
    f_fasta = os.path.join('fasta',pdb+'.fasta')
    with open(f_fasta) as FIN:
        FIN.readline()
        return FIN.readline().strip()


SEQ = {}
for dir in D_SYSTEMS:
    pdb = dir.split('/')[-1].split("_")[0]
    if pdb not in SEQ:
        print pdb
        SEQ[pdb] = load_residues(pdb)

PDB = [dir.split('/')[-1].split("_")[0] for dir in D_SYSTEMS]
#D_SYSTEMS = [x for x,pdb in zip(D_SYSTEMS,PDB) if len(SEQ[pdb])<100]


def run_system(dir):
    print "Loading", dir

    pdb = dir.split('/')[-1].split("_")[0]
    system = '_'.join(dir.split('/')[-1].split("_")[1:3])

    f_traj = os.path.join(dir,"dRMSD.txt")

    Q = np.loadtxt(f_traj)

    # Keep the mean of the last 20 frames
    Q = Q[-cut_frames:]


    cutoff = float(system.split('_')[-1])
    system = system.split('_')[0]


    return pdb, system, cutoff, Q.mean()

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

cols = ['pdb','system','cutoff', 'Q']
data = list(ITR)
df = pd.DataFrame(data=data,columns=cols)



def plot_label(label,offset=0):
    g = df[df.system==label].groupby("cutoff")
    X, data = zip(*g)
    X = np.array(X)
    Y = [x.Q.mean() for x in data]
    Y_std = [x.Q.std() for x in data]
    #plt.errorbar(X+offset,Y,yerr=Y_std,label=label)
    plt.errorbar(X+offset,Y,label=label)
    plt.scatter(X+offset,Y,color='k',alpha=0.5)
    return X,Y

plot_label("GREMLIN")
plot_label("RF",offset=0.01)

NATIVE_X = np.linspace(0.1, 5, 100)
NATIVE_Y = df[df.system=="exact"].Q.mean()
plt.plot(NATIVE_X,[NATIVE_Y,]*100,'k--',alpha=0.5,label="native")

plt.legend(loc="best",fontsize=18)
plt.xlabel("cutoff",fontsize=18)
plt.xlim(0,5.1)
plt.ylim(0,15)
plt.ylabel("dRMSD average for all proteins",fontsize=18)

os.system("mkdir -p figures")
#plt.savefig("figures/Q_avg.png")
plt.show()
