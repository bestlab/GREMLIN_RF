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

#FLAG_SWITCH_Y = True
FLAG_SWITCH_Y = False

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



CLF_MODEL = {}

def load_image_data(pdb):

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
    #Y = clf.predict_proba(X)[:,1]
    Y = clf.predict(X)
    Yp = clf.predict_proba(X)[:,1]
    cutoff_Yp = Yp[np.argsort(Yp)[::-1]][:int((3/2)*N)][-1]
    Y = np.array([1 if yp>cutoff_Yp else 0 for yp in Yp])

    native = np.array([NATIVE[idx] for idx in IDX])

    def compute_conf(Y):
        TP = ((native==1)&(Y==1))
        TN = ((native==0)&(Y==0))
        FP = ((native==0)&(Y==1))
        FN = ((native==1)&(Y==0))

        TP = TP.sum()
        FP = FP.sum()
        FN = FN.sum()
        TN = TN.sum()
        
        conf = {}
        conf["sensitivity"] = TP/(TP+FN)
        #conf["specificity"] = TN/(TN+FP)
        conf["precision"]   = TP/(TP+FP)
        #conf["accuracy"]    = (TP+TN)/(TP+FP+FN+TN)

        return conf

    conf_RF = compute_conf(Y)

    total_native = NATIVE.sum()
    fixed_L = int((3.0/2)*N)

    #idx = np.argsort(W)[::-1]
    #NATIVEX = np.array([NATIVE[i,j] for i,j in IDX])
    
    # sorting by scores
    W = np.array(W)
    cutoff = W[np.argsort(W)[::-1]][:int((3/2)*N)][-1]
    GM = np.zeros((N,N))
    GM[G>=cutoff] = 1
    Y_L = np.array([GM[i,j] for i,j in IDX])

    conf_L = compute_conf(Y_L)

    # All L values
    #FRACTION_NATIVE = np.cumsum(NATIVEX)/total_native
    #ACCURACY = np.cumsum(NATIVEX)/(np.arange(len(NATIVEX))+1)

    # Fixed L prediction
    #F_NATIVE = NATIVEX[:fixed_L].sum()/total_native
    #F_ACC    = NATIVEX[:fixed_L].sum()/float(fixed_L)
    
    # RF prediction
    #P_NATIVE = TP/total_native
    #P_ACC = conf["accuracy"]

    #print "Total native/total predicted", NATIVE.sum()
    #print "RF Accuracy/P_native", P_ACC, P_NATIVE
    #return FRACTION_NATIVE, ACCURACY, P_NATIVE, P_ACC, F_NATIVE, F_ACC

    js = json.dumps(conf_RF,indent=2)
    print "RF", pdb, js

    js = json.dumps(conf_L,indent=2)
    print "L-fixed", pdb, js

    return conf_RF, conf_L
    


#############################

# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

#np.random.shuffle(PDB)

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

PX, PY = [],[]
FX, FY = [],[]

for pdb in PDB:
    print pdb

    conf_RF, conf_L = load_image_data(pdb)
    
    PX.append(conf_L["sensitivity"])
    PY.append(conf_RF["sensitivity"])

    FX.append(conf_L["precision"])
    FY.append(conf_RF["precision"])

    
    #plt.plot(X,Y,alpha=.225,color='k')


line = np.linspace(0,1,100)

figsize = (8,6)

plt.figure(figsize=figsize)
plt.scatter(PX,PY,color='r')#, label=label.format(**args))
plt.legend()
plt.title("Sensitivity")
plt.ylabel("Random Forest Model")
plt.xlabel("3/2 fixed L model")
plt.xlim(0,1.0)
plt.ylim(0,1)
plt.plot(line,line,'k--',alpha=.75,lw=1)


plt.figure(figsize=figsize)
plt.scatter(FX,FY,color='r')#, label=label.format(**args))
plt.legend()
plt.title("Precision")
plt.ylabel("Random Forest Model")
plt.xlabel("3/2 fixed L model")
plt.xlim(0,1.0)
plt.ylim(0,1)
plt.plot(line,line,'k--',alpha=.75,lw=1)

plt.show()

