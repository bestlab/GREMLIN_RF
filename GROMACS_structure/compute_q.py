import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    
    org_dir = os.getcwd()
    os.chdir(dir)

    f_CMAP  = "contact_map.h5"
    f_Q = "Qval.txt"
    
    if not os.path.exists(f_CMAP):
        print "Missing CMAP, extract_contact_maps first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_Q):
        print "Qval file exists, skipping", dir
        os.chdir(org_dir)
        return dir
    
    h5 = h5py.File(f_CMAP,'r')
    CMAP = h5["CMAP"][:]
    OC   = h5["OC"][:]
    h5.close()

    assert(OC.shape == CMAP[0].shape)    

    Q = []
    for cx in CMAP:
        q = ((cx==OC)*(OC==True)).sum()
        q /= float(OC.sum())
        Q.append(q)

    np.savetxt(f_Q,Q)
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
    print "Completed", item


