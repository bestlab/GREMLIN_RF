import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import src.utils as utils
import Bio.PDB.QCPSuperimposer as QC

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))
kernel_window = 2


def run_system(dir):

    pdb = os.path.basename(dir).split('_')[0]
    
    org_dir = os.getcwd()
    os.chdir(dir)

    f_coord = "coord.h5"
    f_RMSD  = "RMSD.txt"
    f_OC = os.path.join("..","..","cmap_coordinates",pdb+'.txt')

    if not os.path.exists(f_OC):
        print "Missing coordinates for cmap", dir
        os.chdir(org_dir)
        return dir
    
    if not os.path.exists(f_coord):
        print "Missing coordinates, extract_coordinates.py first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_RMSD):
        print "RMSD file exists, skipping", dir
        os.chdir(org_dir)
        return dir
    
    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()
    OC = np.loadtxt(f_OC)

    # Move the coordinates to something sensible
    C  -= C.mean(axis=0)
    OC -= OC.mean(axis=0)

    assert(C[0].shape == OC.shape)

    RMSD = []

    Q = QC.QCPSuperimposer()    
    for cx in C:
        Q.set(OC,cx)
        Q.run()
        val = Q.get_rms()
        RMSD.append(val)

    RMSD = np.array(RMSD)
    print dir, RMSD[-20:].mean()

    np.savetxt(f_RMSD,RMSD)
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



