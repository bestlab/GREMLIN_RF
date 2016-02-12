import glob, os, subprocess
import itertools, multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
plt = sns.plt

_PARALLEL = True
MP_CORES = 30
#D_SYSTEMS  = sorted(glob.glob("systems/*exact*"))[:]
#D_SYSTEMS = sorted(glob.glob("systems/*RF_5.00*"))[:]
D_SYSTEMS = sorted(glob.glob("systems/*GREMLIN_5.00*"))[:]

total_frames = 20

def load_residues(pdb):
    f_fasta = os.path.join('fasta',pdb+'.fasta')
    with open(f_fasta) as FIN:
        FIN.readline()
        return FIN.readline().strip()

    
def run_system(dir):
    print "Loading", dir
    
    pdb = dir.split('/')[-1].split("_")[0]
    system = '_'.join(dir.split('/')[-1].split("_")[1:3])

    seq = load_residues(pdb)
    N = len(seq)

    org_dir = os.getcwd()
    os.chdir(dir)

    f_Q = "Qval.txt"
    f_ENERGY = "energy.xvg"
    f_dRMSD  = "dRMSD.txt"
    f_rg = "rg.txt"
    
    for f in [f_Q,f_ENERGY, f_dRMSD,f_rg]:
        if not os.path.exists(f):
            print "Missing", f, dir
            os.chdir(org_dir)
            return None

    U = np.loadtxt(f_ENERGY,comments=["@","#"])[:,-1][-total_frames:]
    dRMSD = np.loadtxt(f_dRMSD)[-total_frames:]
    Q = np.loadtxt(f_Q)[-total_frames:]
    rg = np.loadtxt(f_rg)[-total_frames:]

    seed = int(dir.split('/')[-1].split("_")[4])

    cutoff = float(system.split('_')[-1])
    system = system.split('_')[0]

    os.chdir(org_dir)

    return pdb, system, seed, cutoff, dRMSD.mean(), N, Q.mean(), rg.mean(), U.mean()


def load_data(D_SYSTEMS):

    INPUT_ITR = D_SYSTEMS
    ITR = itertools.imap(run_system, INPUT_ITR)

    if _PARALLEL:
        import multiprocessing
        MP = multiprocessing.Pool(MP_CORES)
        ITR = MP.imap(run_system, INPUT_ITR)

    cols = ['pdb','system','seed', 'cutoff', 'dRMSD', 'N', 'Q','rg', 'U']
    data = list([x for x in ITR if x is not None])
    df = pd.DataFrame(data=data,columns=cols)

    return df

df = load_data(D_SYSTEMS)

x_key = "N"
y_key = "dRMSD"
X = df.groupby('pdb')[x_key].min()
Y = df.groupby('pdb')[y_key].min()

'''
x_key = "U"
y_key = "dRMSD"
X = df.groupby('pdb')[x_key]

xm = df.groupby('pdb')[x_key].min()
PDB = xm.keys()

X = [xm[pdb] for pdb in PDB]
Y = [df[df.pdb == pdb][df.U == xm[pdb]][y_key].values for pdb in PDB]
Y = np.array(Y)
print Y
'''

#Y = df.groupby('pdb')[y_key].min()
#X = df[x_key]
#Y = df[y_key]


plt.scatter(X,Y)
plt.xlabel(x_key,fontsize=18)
plt.ylabel(y_key,fontsize=18)
plt.ylim(0,12)

plt.show()
exit()
