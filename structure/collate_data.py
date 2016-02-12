import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import pandas as pd

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

total_frames = 40

def run_system(dir):
    
    org_dir = os.getcwd()
    os.chdir(dir)
    pdb = os.path.basename(dir).split('_')[0]
    system = '_'.join(dir.split('/')[-1].split("_")[1:3])
    
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

    seed = [seed,]*total_frames
    cutoff = [cutoff,]*total_frames
    pdb = [pdb]*total_frames
    system = [system]*total_frames

    data = zip(*[pdb,system,cutoff,seed,U,Q,rg,dRMSD])
    os.chdir(org_dir)

    return data

            

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)


data = []

for name in tqdm.tqdm(INPUT_ITR):
    result = ITR.next()

    if result is None:
        continue

    for item in result:
        data.append(item)


print "Building dataframe"
    
cols = "pdb","system","cutoff","seed","U","Q","rg","dRMSD"
df = pd.DataFrame(data=data,columns=cols)

print "Saving to disk"
os.system('mkdir -p collated_data')
f_csv = os.path.join('collated_data','system_data.csv')
df.to_csv(f_csv)

print df.shape
