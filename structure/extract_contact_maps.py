import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
from scipy.spatial.distance import cdist, pdist, squareform

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

distance_cutoff = 8.0

def run_system(dir):

    pdb = dir.split('/')[-1].split("_")[0]
    f_cmap = os.path.join('cmap',pdb+'.cmap')
    OC = np.loadtxt(f_cmap,dtype=bool)
    
    org_dir = os.getcwd()
    os.chdir(dir)

    f_coord = "coord.h5"
    f_CMAP  = "contact_map.h5"
    
    if not os.path.exists(f_coord):
        print "Missing coords, extract_coord first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_CMAP):
        print "Coord file exists, skipping", dir
        os.chdir(org_dir)
        return dir

    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()

    sample_N = C.shape[0]
    N = C.shape[1]

    CMAP = np.zeros((sample_N,N,N),dtype=bool)

    for i,cx in enumerate(C):
        dr  = squareform(pdist(cx))
        idx = dr < distance_cutoff
        CMAP[i] = idx
    

    '''
    for cx in C[::10]:
        dr  = squareform(pdist(cx))
        idx = dr < distance_cutoff
        print cx
        import seaborn as sns
        sns.heatmap(dr)
        sns.plt.show()
        #exit()
        print 
    '''

    h5 = h5py.File(f_CMAP,'w')
    h5.create_dataset("CMAP",data=CMAP, compression="gzip")
    h5.create_dataset("OC",data=OC, compression="gzip")
    h5.close()
    
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


