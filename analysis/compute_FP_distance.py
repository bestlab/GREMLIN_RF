from __future__ import division
from src.utils import fixedL_cut
from src.utils import generate_matrix_IDX
import numpy as np
import pandas as pd
import itertools,os,glob
import src.utils as utils

from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

kernel_window = 2
cut_iterations = 500

_PARALLEL = True
MP_CORES = 30

def build_cmap_from_index(N, IDX, contacts):
    C = np.zeros((N,N))

    for idx_idx in np.where(contacts)[0]:
        i,j = IDX[idx_idx]
        C[i,j] = C[j,i] = 1
    return C
    
def compute_cmap(f_rank):
    base = os.path.basename(f_rank)
    model_name = os.path.dirname(f_rank)
    pdb  = base.split('.')[0]

    func = {
        "G2" : utils.load_improved_GREMLIN_dataset,
        "APC": utils.load_GREMLIN_dataset,
    }
        
    G = func[model_name](pdb)
    NATIVE_MATRIX = utils.load_contact_map(pdb)

    N   = G.shape[0]
    IDX = list(generate_matrix_IDX(N,kernel_window))

    # Load the native contacts
    native = np.array([NATIVE_MATRIX[idx] for idx in IDX])
    g      = np.array([G[idx] for idx in IDX])

    data = {}
    CUT_IDX = np.linspace(1,5*N,cut_iterations)
    
    for cut_idx in CUT_IDX:
        contacts = fixedL_cut(g,native,cut_idx,index_only=True)
        C = build_cmap_from_index(N, IDX,contacts)
        data[cut_idx/N] = C

    return data,NATIVE_MATRIX

#import pylab as plt
#import matplotlib.pylab as plt
#import seaborn as sns

os.system('mkdir -p FP_distance')
F_CMAP = glob.glob("cmap/*.cmap")
PDB = [os.path.basename(f).split('.')[0] for f in F_CMAP]

def compute_FP(pdb):

    f_FP = "FP_distance/{}_FP_distance.txt".format(pdb)

    if os.path.exists(f_FP):
        print "f_FP exists, skipping", f_FP
        return pdb

    print "Starting", pdb
            
    args = {"pdb" : pdb}
    f_APC = "APC/{pdb}.gremlin".format(**args)
    f_G2  = "G2/{pdb}.gremlin".format(**args)

    if not os.path.exists(f_APC) or not os.path.exists(f_G2):
        print "Missing gremlin", f_APC
        return pdb
    

    C1,NATIVE = compute_cmap(f_APC)
    C2,NATIVE = compute_cmap(f_G2)

    

    NATIVE_COORDS = np.array(np.where(NATIVE)).T
    T_NATIVE = cKDTree(NATIVE_COORDS)

    data = []

    for key in sorted(C2.keys()):
        cx1_coords = np.array(np.where(C1[key])).T
        cx2_coords = np.array(np.where(C2[key])).T
        dist1,_ = T_NATIVE.query(cx1_coords)
        dist2,_ = T_NATIVE.query(cx2_coords)

        data.append( [key, dist1.mean(), dist2.mean()] )

    np.savetxt(f_FP, data)

    return pdb


ITR = itertools.imap(compute_FP, PDB)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(compute_FP, PDB)

for pdb in ITR:
    print "Completed", pdb


