import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import src.utils as utils


_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))
_FORCE = True

NATIVE_CUTOFF_LENGTH = 8.0
kernel_window = 2

def compute_rg(X):
    N = X.shape[0]
    mu = X.mean(axis=0)
    rg2 = np.linalg.norm(X - mu,axis=1).mean()
    return np.sqrt(rg2)

def run_system(dir):

    pdb = os.path.basename(dir).split('_')[0]
    
    org_dir = os.getcwd()
    os.chdir(dir)

    f_coord = "coord.h5"
    f_out   = "rg.txt"
    f_OC = os.path.join("..","..","cmap_coordinates",pdb+'.txt')

    if not os.path.exists(f_OC):
        print "Missing coordinates for cmap", dir
        os.chdir(org_dir)
        return dir
    
    if not os.path.exists(f_coord):
        print "Missing coordinates, extract_coordinates.py first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_out) and not _FORCE:
        print "rg file exists, skipping", dir
        os.chdir(org_dir)
        return dir
    
    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()
    OC = np.loadtxt(f_OC)
    N = OC.shape[0]

    rg_OC = compute_rg(OC)

    RG = []

    for cx in C:
        rg = compute_rg(cx)
        RG.append(rg)

    RG = np.array(RG)

    RG /= rg_OC
    
    #print dir, RMSD[-20:].mean(), org_RMSD[-20:].mean(),RG[-20:].mean()
    print "{} {: 0.4f}".format(dir, RG[-200:].mean())

    #OUTPUT = np.vstack([RG,RG/rg_OC])
    np.savetxt(f_out,RG)
    os.chdir(org_dir)

    return dir


INPUT_ITR = D_SYSTEMS[:]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

#for item in tqdm.tqdm(ITR):
for item in ITR:
    #print "Completed", item
    pass



exit()
#######################################################################
#######################################################################

import glob, os, subprocess
import itertools, multiprocessing
import tqdm

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    f_traj = os.path.join(dir,"gyrate.xvg")
    if os.path.exists(f_traj):
        print "Already computed", dir
        return dir
    
    cmd = ("cd {}; echo 0 | g_gyrate -f traj.xtc -s topol.tpr")
    cmd = cmd.format(dir)
    
    with open(os.devnull, 'w') as shutup:
        try:
            subprocess.check_call(cmd, stdout=shutup,stderr=shutup,
                                shell=True)
        except:
            print "ERROR", dir
            return dir

    return dir

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for k,f in enumerate(tqdm.tqdm(ITR, total=len(INPUT_ITR))):
    pass


#for d in D_SYSTEMS:
#    run_system(d)
#    pass
