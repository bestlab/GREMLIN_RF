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

    pdb = os.path.basename(dir).split('_')[0]
    
    org_dir = os.getcwd()
    os.chdir(dir)

    f_coord = "coord.h5"
    f_RMSD  = "dRMSD.txt"
    f_OC = os.path.join("..","..","cmap_coordinates",pdb+'.txt')

    if not os.path.exists(f_OC):
        print "Missing coordinates for cmap", dir
        os.chdir(org_dir)
        return dir
    
    if not os.path.exists(f_coord):
        print "Missing coordinates, extract_coordinates.py first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_RMSD) and not _FORCE:
        print "RMSD file exists, skipping", dir
        os.chdir(org_dir)
        return dir
    
    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()
    OC = np.loadtxt(f_OC)
    N = OC.shape[0]

    # Identify the native contacts
    d_OC = cdist(OC,OC)
    native_map = d_OC < NATIVE_CUTOFF_LENGTH

    # Don't consider a native contacts as i with itself
    for i in range(kernel_window,N-kernel_window):
        for j in range(kernel_window+1):
            native_map[i,i+j] = False
            native_map[i,i-j] = False
            
    # Don't consider the edges for native calculation
    native_map[:kernel_window,:] = False
    native_map[-kernel_window:,:] = False

    native_map[:,:kernel_window] = False
    native_map[:,-kernel_window:] = False

    native_N = native_map.sum()

    # Only consider the native contacts
    d_OC *= native_map

    # Cut coordinates since endpoints aren't contrainted
    #OC = OC[2:-2,:]

    DRMSD = []

    for cx in C:
        
        d_cx = cdist(cx,cx)
        d_cx *= native_map
        
        dist = (d_cx - d_OC)**2
    
        #d_rmsd = (dist.sum() / (N*(N-1)))**(0.5)
        d_rmsd = (dist.sum() / native_N)**(0.5)

        DRMSD.append(d_rmsd)

    #rot, tran = sup.get_rotran()
    #cx = np.dot(cx, rot) + tran
    #org_RMSD = np.array(org_RMSD)
    #RG = np.array(RG)

    DRMSD = np.array(DRMSD)
    
    #print dir, RMSD[-20:].mean(), org_RMSD[-20:].mean(),RG[-20:].mean()
    print "{} {: 0.4f}".format(dir, DRMSD[-200:].mean())
    
    np.savetxt(f_RMSD,DRMSD)
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



