import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import MDAnalysis

_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    org_dir = os.getcwd()
    os.chdir(dir)

    f_movie = "movie.pdb"
    f_coord = "coord.h5"
    
    if not os.path.exists(f_movie):
        print "Missing movie, extract_pdb first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_coord):
        print "Coord file exists, skipping", dir
        os.chdir(org_dir)
        return dir


    P = MDAnalysis.Universe(f_movie, multiframe=True)
    T = P.trajectory

    C = []
    for k,x in enumerate(T):
        y = np.array(x.positions)
        C.append(y)
        #if k>10: break
    C = np.array(C)

    test_val = ((C[0][0] - C[1][0])**2).sum()
    print "Test val", test_val
    assert( test_val )

    #####################################################################
        
    h5 = h5py.File(f_coord,'w')
    h5["coord"] = C
    h5.close()
    
    os.chdir(org_dir)
    return dir


INPUT_ITR = D_SYSTEMS[:]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for item in tqdm.tqdm(ITR):
    print "Completed", item
    pass

