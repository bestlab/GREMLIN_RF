import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import src.utils as utils
#import Bio.PDB.QCPSuperimposer as QC
from Bio.SVDSuperimposer import SVDSuperimposer 

_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

_FORCE = True

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

    if os.path.exists(f_RMSD) and not _FORCE:
        print "RMSD file exists, skipping", dir
        os.chdir(org_dir)
        return dir
    
    h5 = h5py.File(f_coord,'r')
    C = h5["coord"][:]
    h5.close()
    OC = np.loadtxt(f_OC)

    # Move the coordinates to something sensible
    #C  -= C.mean(axis=0)
    #OC -= OC.mean(axis=0)

    median_OC = np.median([np.linalg.norm(a-b)
                           for a,b in zip(OC,OC[1:])])
    median_C  = np.median([np.linalg.norm(a-b)
                           for a,b in zip(C[-1],C[-1][1:])])

    assert(C[0].shape == OC.shape)
    RMSD = []
    org_RMSD = []

    sup = SVDSuperimposer()

    RG = []
    OC -= OC.mean(axis=0)
    OC_RG = ((np.linalg.norm(OC,axis=1)**2).sum()/len(OC)) ** 0.5

    for cx in C:
        cx -= cx.mean(axis=0)

        rg_cx = ((np.linalg.norm(cx,axis=1)**2).sum()/len(cx)) ** 0.5
        RG.append(rg_cx)
        
        sup.set(OC,cx)
        sup.run()
        RMSD.append(sup.get_rms())
        org_RMSD.append(sup.get_init_rms())


    rot, tran = sup.get_rotran()
    cx = np.dot(cx, rot) + tran

    RMSD = np.array(RMSD)
    org_RMSD = np.array(org_RMSD)
    RG = np.array(RG)
    
    #print dir, RMSD[-20:].mean(), org_RMSD[-20:].mean(),RG[-20:].mean()
    print "{} {: 0.4f} {: 0.4f}".format(dir, RMSD[-20:].mean(),
                                      RG[-20:].mean() / OC_RG)
    

    '''
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(OC[:,0],OC[:,1],OC[:,2],'b')
    #ax.plot(OC[:,0],OC[:,1],OC[:,2],'k',alpha=0.5)

    ax.scatter(cx[:,0],cx[:,1],cx[:,2],color='r')
    #ax.plot(cx[:,0],cx[:,1],cx[:,2],'k',alpha=0.5)
    plt.show()
    exit()

    print OC
    #exit()
    '''
    
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



