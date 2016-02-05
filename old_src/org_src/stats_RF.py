from __future__ import division
import h5py, itertools, os, json, gc, glob, argparse
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.externals import joblib
from src.utils import APC_L2

#clf_dir = 'clf2'
clf_dir = 'clf'

diag_window = 5
debug_cutoff = 10**10
#debug_cutoff = 10

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-r", "--ratio_TP_to_FP", type=int,default=40)
args = parser.parse_args()

kernel_window = args.kernel_window
ratio_TP_to_FP = args.ratio_TP_to_FP

args = {
    "kernel_window":args.kernel_window,
    "ratio_TP_to_FP":args.ratio_TP_to_FP,
    "clf_name":"stats_ExtraTreesClassifier",
}



def compute_conf(Y,native):
    TP = ((native==1)&(Y==1)).sum()
    #TN = ((native==0)&(Y==0)).sum()
    FP = ((native==0)&(Y==1)).sum()
    FN = ((native==1)&(Y==0)).sum()

    conf = {}
    conf["sensitivity"] = TP/(TP+FN)
    conf["precision"]   = TP/(TP+FP)
    return conf

def fixedL_cut(W,native,cut_idx):

    cutoff = W[np.argsort(W)[::-1]][:cut_idx]
    if cutoff.size:
        cutoff = cutoff[-1]
    else:
        cutoff = W.max()

    W2 = np.zeros(W.shape)
    W2[W>=cutoff] = 1

    return compute_conf(W2,native)
    

CLF_MODEL = {}

def compute_RF_stats(pdb):

    global CLF_MODEL
    
    G = h5["PDB"][pdb]["gremlin_matrix"]

    # Load the contact map
    N = G.shape[0]

    NATIVE = np.zeros(shape=(N,N))
    dx = kernel_window

    for i,j in h5["PDB"][pdb]["native_contacts"][:]:
        i,j = sorted([i,j])[::-1]
    
        # Don't take ultra close range contacts near or below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        NATIVE[i,j] = NATIVE[j,i] =1


    f_clf = PDB_F_CLF[pdb]
    
    if f_clf not in CLF_MODEL:
        print "Loading model", f_clf

        # Remove an old model for memory sake
        if CLF_MODEL:
            key = CLF_MODEL.keys().pop()
            del CLF_MODEL[key]
            gc.collect()

        CLF_MODEL[f_clf] = joblib.load(f_clf)
    
    clf = CLF_MODEL[f_clf]
        
    #print "Loading GREMLIN data"
    G = h5["PDB"][pdb]["gremlin_matrix"][:]
    G = APC_L2(G)

    #################################################################

    #print "Segmenting GREMLIN images"
    IDX = []; X = []; W = []
    for i,j in itertools.product(range(dx,N-dx),repeat=2):
        # Don't take ultra close range contacts near or below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()
        IDX.append((i,j))
        X.append(image)
        W.append(G[i,j])

    print "Running prediction"
    native = np.array([NATIVE[idx] for idx in IDX])

    #Y = clf.predict(X)
    Yp = clf.predict_proba(X)[:,1]
        
    # sorting by scores
    for cut_idx in range(1,len(IDX)):
        conf = fixedL_cut(Yp,native,cut_idx)
        yield cut_idx, conf["sensitivity"], conf["precision"]
    


#############################

# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

# Remove a few bad proteins
PDB = np.array([pdb for pdb in PDB if pdb not in ["1jbe","1hdo","1tqh","1dqg","1fk5","1gbs","2phy","1xkr","1kw4","1i1j","1kqr"]])

# For the PDB list, find out which classifier they belong to
F_CLF_JSON = glob.glob("{}/*window_{}_*_ratio_{}.json".format(clf_dir,
                                                              kernel_window,
                                                              ratio_TP_to_FP))
PDB_F_CLF = {}

for f_json in F_CLF_JSON:
    with open(f_json) as FIN:
        js = json.loads(FIN.read())

    f_clf = js["f_clf"]
    for pdb in js["test_pdb"]:
        PDB_F_CLF[pdb] = f_clf


os.system('mkdir -p stats')
os.system('mkdir -p stats/RF')

for pdb in PDB:
    f_save = os.path.join('stats','RF',pdb+'.txt')
    if os.path.exists(f_save):
        print "Already computed", pdb
        continue
    
    with open(f_save,'w') as FOUT:
        for item in compute_RF_stats(pdb):
            print "{} {:d} {:0.3f} {:0.3f}".format(pdb, *item)
            FOUT.write("{:0.8f} {:0.8f}\n".format(*item[1:]))

