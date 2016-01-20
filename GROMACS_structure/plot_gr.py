import glob, os, subprocess
import itertools, multiprocessing
import tqdm
import numpy as np

_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    f_traj = os.path.join(dir,"gyrate.xvg")
    data = []
    with open(f_traj) as FIN:
        for line in FIN:
            if line[0] in "#@":
                continue
            line = map(float,line.split())
            data.append(line)
            
    data = np.array(data)
    # Keep only the first two columns
    data = data[:,[0,1]]
    return data

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

data = [result for result in ITR]

print "{} total samples to plot".format(len(data))

import seaborn as sns
plt = sns.plt

for k,row in enumerate(data):
    x,y = row.T
    y /= y[0]
    plt.plot(x,y,alpha=0.5,color='k')

plt.title("Radius of gyration for all proteins/contact maps")   
plt.xlabel('ps')
plt.ylabel('rg/rg(t=0)')
plt.tight_layout()
plt.savefig("figures/rg_sample.png")
plt.show()




