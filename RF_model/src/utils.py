from __future__ import division
import numpy as np
import itertools
import h5py, os, glob, json

#V = np.loadtxt("kernel_test/V.txt")

def APC_fro(X):
    score  = X
    row    = np.mean(score,axis=0)
    weight = np.outer(row,row.T) / score.mean()
    final = score-weight
    final = (final+final.T)/2
    return final

def APC_L2(G, project_N_N=False):
    if project_N_N:
        G2 = np.linalg.norm(G,ord="fro",axis=(2,3))
    else:
        G2  = G
        
    G2 = APC_fro(G2)

    # This is data whitening!
    G2 -= G2.mean()
    G2 /= G2.std()

    N = G.shape[0]
    G2 /= np.sqrt(N)

    return G2

def generate_matrix_IDX(N,dx):
    '''
    Returns the matrix indices of the upper diagonal, not
    near an edge or diagonal by a distance of dx.
    Also returns indices that are "close-range".
    '''
    IDX = set()
    IDX_close = set()
    diag_offset = 0
    diag_inner_offset = 0

    for i,j in itertools.product(range(dx,N-dx),repeat=2):
    
        # Don't take ultra close range contacts near
        # or below the diagonal
        
        if i<=j+dx + diag_offset: continue
        if abs(i-j) < diag_inner_offset: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        IDX.add((i,j))

    return IDX

def generate_feature_vectors(G,seq,IDX,kernel_window):

    dx = kernel_window # This is symmetric so image size is 2*dx + 1
    image_shape  = (2*dx+1,2*dx+1)
    X = []

    # Project down to the basic gaussian kernel
    #v0 = V[0]#.reshape(image_shape)
    #v1 = V[1]#.reshape(image_shape)

    for i,j in IDX:
        # This loads the block of size (2*dx+1,2*dx+1,20,20)
        # Unravel the 20x20 residue channel
        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()

        #two_hot_encoding = [0,]*20
        #two_hot_encoding[seq[i]] += 1
        #two_hot_encoding[seq[j]] += 1
        #image = image.tolist() + two_hot_encoding
        
        X.append(image)
        
    return X

def compute_conf(Y,native, IDX_CLOSE=None):
    '''
    Computes the confusion matrix for a given prediction and the native state.
    '''
    
    TP = ((native==1)&(Y==1)).sum()
    FP = ((native==0)&(Y==1)).sum()
    FN = ((native==1)&(Y==0)).sum()
    TN = ((native==0)&(Y==0)).sum()

    conf = {}
    conf["sensitivity"] = TP/(TP+FN)
    conf["precision"]   = TP/(TP+FP)
    conf["specificity"] = TN/(TN+FP)
    conf["negative_predictive_value"] = TN/(TN+FN)

    return conf

def fixedL_cut(W,native,cut_idx):
    '''
    Intermediate step when computing the confusion matrix. Given a fixed
    cut index (say the Nth top ranked score), this will find all matching
    residues and compute the sensitivity/precision.
    '''

    cutoff = W[np.argsort(W)[::-1]][:cut_idx]
    if cutoff.size:
        cutoff = cutoff[-1]
    else:
        cutoff = W.max()

    W2 = np.zeros(W.shape)
    W2[W>=cutoff] = 1

    return compute_conf(W2,native)
    

def load_dataset_lookup(clf_dir, bad_protein_list=[],debug_cutoff=10**10):
    '''
    Loads the list of proteins and their corresponding clf files.
    Useful since the models are stored into separate files for the folds.
    '''

    h5 = h5py.File("data/gremlin_data.h5",'r')
    PDB = h5["PDB"].keys()[:debug_cutoff]
    h5.close()

    # Remove a few bad proteins
    PDB = np.array([pdb for pdb in PDB if pdb
                    not in bad_protein_list])

    # For the PDB list, find out which classifier they belong to
    #f_clf_base = os.path.join(clf_dir,"*window_{}_*_ratio_{}.json")
    #f_clf_glob = f_clf_base.format(kernel_window,ratio_TP_to_FP)
    
    f_clf_glob = os.path.join(clf_dir,"*.json")
    F_CLF_JSON = glob.glob(f_clf_glob)

    PDB_CLF = {}

    for f_json in F_CLF_JSON:
        with open(f_json) as FIN:
            js = json.loads(FIN.read())
            f_clf = js["f_clf"]
        for pdb in js["test_pdb"]:
            PDB_CLF[pdb] = f_clf

    return PDB_CLF
