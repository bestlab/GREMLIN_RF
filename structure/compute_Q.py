import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import src.utils as utils
from scipy.spatial.distance import pdist, cdist

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))
_FORCE = True

NATIVE_CUTOFF_LENGTH = 8.0
kernel_window = 2

def run_system(dir):
    
    org_dir = os.getcwd()
    os.chdir(dir)
    pdb = os.path.basename(dir).split('_')[0]
    
    #f_CMAP  = "contact_map.h5"
    f_Q = "Qval.txt"
    f_coord = "coord.h5"
    f_OC = os.path.join("..","..","cmap_coordinates",pdb+'.txt')
    
    #if not os.path.exists(f_CMAP):
    #    print "Missing CMAP, extract_contact_maps first", dir
    #    os.chdir(org_dir)
    #    return dir
    
    if not os.path.exists(f_coord):
        print "Missing coordinates, extract_coordinates.py first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_Q) and not _FORCE:
        print "Qval file exists, skipping", dir
        os.chdir(org_dir)
        return dir

    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()
    OC = np.loadtxt(f_OC)

    #h5 = h5py.File(f_CMAP,'r')
    #CMAP = h5["CMAP"][:]
    #OC   = h5["OC"][:]
    #assert(OC.shape == CMAP[0].shape)
    #h5.close()

    N = OC.shape[0]
    
    # Identify the native contacts
    d_OC = cdist(OC,OC)
    native_map = d_OC < NATIVE_CUTOFF_LENGTH

    # Don't consider a native contacts as i+kernel_window with itself
    for i in range(kernel_window,N-kernel_window):
        for j in range(kernel_window+1):
            native_map[i,i+j] = False
            native_map[i,i-j] = False
            
    # Don't consider the edges for native calculation
    native_map[:kernel_window,:] = False
    native_map[-kernel_window:,:] = False

    native_map[:,:kernel_window] = False
    native_map[:,-kernel_window:] = False
    
    
    native_N = float( native_map.sum() )

    Q = []
    for cx in C:
        
        cx_contacts = cdist(cx,cx) < NATIVE_CUTOFF_LENGTH
        q = (cx_contacts*native_map).sum() / native_N
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



