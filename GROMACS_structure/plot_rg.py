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

    if not os.path.exists(f_traj):
        return None
    
    with open(f_traj) as FIN:
        for line in FIN:
            if line[0] in "#@":
                continue
            line = map(float,line.split())
            data.append(line)
            
    data = np.array(data)

    # Keep only the first two columns
    data = data[:,[0,1]]
    return dir,data

INPUT_ITR = D_SYSTEMS[:20]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

data = []
for dir,result in ITR:
    if result is None: continue
    name = os.path.basename(dir)
    tokens = name.split('_')
    pdb, style,length,_,seed = tokens

    # Only take the last 25% of frames
    frame_n = result.shape[0]
    keepframes = int(frame_n*0.25)
    result = result[-keepframes:]

    data.append(result)


print "{} total samples to plot".format(len(data))

import seaborn as sns
plt = sns.plt

for k,row in enumerate(data):
    x,y = row.T
    y /= y[0]
    plt.plot(x,y,alpha=0.05,color='k')

plt.title("Radius of gyration for all proteins/contact maps")   
plt.xlabel('ps')
plt.ylabel('rg/rg(t=0)')
plt.tight_layout()
plt.savefig("figures/rg_sample.png")
plt.show()




