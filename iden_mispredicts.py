from __future__ import division
import h5py, itertools, os, json, gc, glob, argparse
import numpy as np
import argparse
from sklearn.externals import joblib
from src.utils import APC_L2

os.system('mkdir -p missed_contacts')

diag_window = 5
debug_cutoff = 10**10

cutoff_N_count = 10

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-r", "--ratio_TP_to_FP", type=int,default=40)
args = parser.parse_args()

kernel_window = args.kernel_window
ratio_TP_to_FP = args.ratio_TP_to_FP
        
CLF_MODEL = {}
pdb_cutoff = 10**10

def load_image_data(pdb):

    f_missed = os.path.join("missed_contacts",pdb+'.txt')
    if os.path.exists(f_missed):
        #print f_missed, "already computed"
        return False
    print "Working on", f_missed
    
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

    #print "Segmenting GREMLIN images"
    IDX = []; X = []
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

    # Running prediction
    Y = clf.predict(X)
    score = clf.predict_proba(X)
    del X

    Y0 = score[Y==0,0]
    Y1 = score[Y==1,0]
    native = np.array([NATIVE[idx] for idx in IDX])

    TP = ((native==1)&(Y==1))
    TN = ((native==0)&(Y==0))
    FP = ((native==0)&(Y==1))
    FN = ((native==1)&(Y==0))

    #select_class = (FN|FP)
    select_class = (TP|FP|FN|TN)
    select_pos_score = score[select_class][:,1]

    NATIVE_N = NATIVE.sum()
    cutoff_N = cutoff_N_count*int(NATIVE_N)

    print "{} TP+FP, {} TN+FN scores".format(TP.sum()+FP.sum(),FN.sum()+TN.sum())
    print "{} native contacts".format(NATIVE_N)
    print "{} cutoff_N".format(cutoff_N)

    # Find the challenging scores, those closest to 0.5
    diff_score = np.abs(select_pos_score-0.5)
    
    challenge_args = np.argsort(diff_score)
    delta_score_cutoff = diff_score[challenge_args][:cutoff_N].max()
    print "Score cutoff +/- ", delta_score_cutoff

    NEW_SET = []
    NEW_labels = []

    for (i,j),is_predicted,s in zip(IDX,Y,score[:,1]):

        is_native = bool(NATIVE[i,j])

        is_new = False

        if np.abs(s-0.5) < delta_score_cutoff:
            is_new = True
            #if is_native and is_predicted:
            #    is_new = True
            #if not is_native and is_predicted:
            #    is_new = True
            #if is_native and not is_predicted:
            #    is_new = True
                
        if is_new:
            NEW_SET.append((i,j))
            NEW_labels.append(is_native)

    pos = sum(NEW_labels)
    neg = len(NEW_labels) - pos
    print "Labels {} P {} N".format(pos, neg)
    
    #exit()

    '''    
    TP = ((native==1)&(Y==1)).sum()
    TN = ((native==0)&(Y==0)).sum()
    FP = ((native==0)&(Y==1)).sum()
    FN = ((native==1)&(Y==0)).sum()

    conf = {}
    conf["sensitivity"] = TP/(TP+FN)
    conf["specificity"] = TN/(TN+FP)
    conf["precision"]   = TP/(TP+FP)
    conf["accuracy"]    = (TP+TN)/(TP+FP+FN+TN)

    js = json.dumps(conf,indent=2)
    print js
    '''
    
    np.savetxt(f_missed, NEW_SET,fmt="%i")

    return None


#############################

# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

#np.random.shuffle(PDB)

# Remove a few bad proteins
# 1fk5, 1dqg, 2phy, 1xkr 1kw4 1i1j 1kqr have no TP!?
PDB = np.array([pdb for pdb in PDB if pdb not in ["1jbe","1hdo","1tqh","1dqg","1fk5","1gbs","2phy","1xkr","1kw4","1i1j","1kqr"]])

# For the PDB list, find out which classifier they belong to
F_CLF_JSON = glob.glob("clf/*window_{}_*_ratio_{}.json".format(kernel_window,ratio_TP_to_FP))
PDB_F_CLF = {}

for f_json in F_CLF_JSON:
    with open(f_json) as FIN:
        js = json.loads(FIN.read())

    f_clf = js["f_clf"]
    for pdb in js["test_pdb"]:
        PDB_F_CLF[pdb] = f_clf

PX, PY = [],[]
FX, FY = [],[]

for pdb in PDB[:pdb_cutoff]:
    try:
        load_image_data(pdb)
    except ValueError:
        print "ERRRRRRRRRROR", pdb

