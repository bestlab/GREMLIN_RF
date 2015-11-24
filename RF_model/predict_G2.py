from __future__ import division
import itertools, os, glob
import numpy as np
from src.utils import APC_L2, generate_matrix_IDX
from src.utils import load_dataset_lookup
from src.utils import generate_feature_vectors
import src.utils as utils
from sklearn.externals import joblib

clf_dir = 'clf'
MP_CORES = 20
kernel_window = 2

def compute_predictions(pdb,clf):

    # Load the GREMLIN data
    g   = utils.load_GREMLIN_dataset(pdb)
    fasta,seq = utils.load_seq(pdb)

    # Sanity checks    
    assert(len(seq) == g.shape[0])

    # APC corrections and whiten
    g = APC_L2(g)
    N = g.shape[0]

    IDX = generate_matrix_IDX(N,kernel_window)

    # Load the native contacts
    #NATIVE_MATRIX = utils.load_contact_map(pdb)
    #native = [NATIVE_MATRIX[idx] for idx in IDX]

    #################################################################

    X = generate_feature_vectors(g,seq,IDX,kernel_window)
    Yp = clf.predict_proba(X)[:,1]

    g2 = np.zeros(g.shape)
    for (i,j),y in zip(IDX,Yp):
        g2[i,j] = g2[j,i] = y

    '''
    # Reorder based off of ranking
    order = np.argsort(Yp)[::-1]
    IDX0 = np.array(map(list,IDX))
    IDX0 = IDX0[order]
    
    W = np.array([G[i,j] for i,j in IDX])
    order = np.argsort(W)[::-1]
    IDX1 = np.array(map(list,IDX))
    IDX1 = IDX1[order]
    '''

    return g2
        
#############################

CLF_DIR = {}

def compute_and_save_prediction((pdb,f_clf)):
    global CLF_DIR
    
    f_save = os.path.join('G2',pdb+'.g2.gremlin')

    if os.path.exists(f_save):
        print "Already computed", pdb
        #return pdb
    
    if f_clf not in CLF_DIR:
        # Unload other models
        CLF_DIR = {}
        
        print "Loading model", f_clf
        CLF_DIR[f_clf] = joblib.load(f_clf)
        
    clf = CLF_DIR[f_clf]
    g2 = compute_predictions(pdb,clf)
    
    np.savetxt(f_save, g2, fmt="%0.6f")
    return pdb

os.system('mkdir -p G2')
#os.system('mkdir -p predictions/RF')
#os.system('mkdir -p predictions/org')

PDB_CLF = load_dataset_lookup(clf_dir, [])
INPUT_ITR = PDB_CLF.items()

#import multiprocessing
#P = multiprocessing.Pool(MP_CORES)

func = compute_and_save_prediction
ITR = itertools.imap(func, INPUT_ITR)
#ITR = P.imap(func, INPUT_ITR)

for pdb in ITR:
    print "Finished", pdb
