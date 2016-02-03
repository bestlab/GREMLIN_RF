import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import src.utils as utils

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))
kernel_window = 2

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
    assert(OC.shape == CMAP[0].shape)
    h5.close()

    N = CMAP.shape[1]
    IDX = utils.generate_matrix_IDX(N,kernel_window)
    VALID_IDX = np.zeros((N,N),dtype=bool)
    for i,j in IDX:
        VALID_IDX[i,j] = VALID_IDX[j,i] = True

    # To ignore VALID_IDX keep this line in
    VALID_IDX = np.ones(VALID_IDX.shape)    

    Q = []
    for cx in CMAP:
        
        q = ((cx==OC)*(OC*VALID_IDX)).sum()

        norm = float((OC*VALID_IDX).sum())
        q /= norm
        
        Q.append(q)

    print dir, np.array(Q).mean()

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
    #print "Completed", item
    pass



