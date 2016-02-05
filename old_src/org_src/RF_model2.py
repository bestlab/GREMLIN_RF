from __future__ import division
import h5py, itertools, os, json, gc
import numpy as np
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import ExtraTreesClassifier as Classifier
import argparse
from src.subsampling import balanced_subsample
from src.utils import APC_L2

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-n", "--n_estimators", type=int,default=200)
parser.add_argument("-r", "--ratio_TP_to_FP", type=int,default=2)
args = parser.parse_args()

k_fold_n = 4
kernel_window = args.kernel_window
n_estimators = args.n_estimators
clf_name =  Classifier.__name__

ratio_TP_to_FP = args.ratio_TP_to_FP

debug_cutoff = 10**10
diag_window = 5

debug_cutoff = 10*100
#debug_cutoff = 10
subsample_tasks = 1

#######################################################################

def load_image_data(pdb):

    f_missed = os.path.join("missed_contacts",pdb+'.txt')
    missed_contacts = np.loadtxt(f_missed)
    NC = missed_contacts[:,[0,1]].astype(int)
    #score = missed_contacts[:,2]

    G = h5["PDB"][pdb]["gremlin_matrix"]
    N = G.shape[0]
    
    dx = kernel_window # This is symmetric so image size is 2*dx + 1

    # Load the contact map
    NATIVE = np.zeros(shape=(N,N))

    for i,j in h5["PDB"][pdb]["native_contacts"][:]:
        i,j = sorted([i,j])[::-1]
    
        # Don't take ultra close range contacts near or
        # below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        NATIVE[i,j] = NATIVE[j,i] =1
        
    Y = np.array([NATIVE[idx[0],idx[1]] for idx in NC],dtype=int)

    #print pdb, "label balance", Y.sum()/len(Y)
    
    X = []
    G = APC_L2(G[:])

    dx = kernel_window # This is symmetric so image size is 2*dx + 1
    image_shape  = (2*dx+1,2*dx+1)

    for i,j in NC:
        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()
        X.append(image)

    del G
    gc.collect()

    # Swap axes so the dim is (N samples, K channels, width, height)
    X = np.array(X)

    return X,Y

def load_fold_dataset(PDB_LIST):
    print "Loading {} protein GREMLIN data".format(len(PDB_LIST))
    X,Y = [], []
    for pdb in PDB_LIST:
        x,y = load_image_data(pdb)
        X.append(x)
        Y.append(y)

    return np.vstack(X), np.hstack(Y)


# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

# 1fk5, 1dqg, 2phy, 1xkr 1kw4 1i1j 1kqr have no TP!?
PDB = np.array([pdb for pdb in PDB if pdb not in ["1jbe","1hdo","1tqh","1dqg","1fk5","1gbs","2phy","1xkr","1kw4","1i1j","1kqr"]])


# Partition the proteins into cross-validation sets
folds_idx = list(KFold(len(PDB), n_folds = k_fold_n))

# Turn the folds into PDB lists
PDB_folds = [{
    "test"  :PDB[train],
    "train" :PDB[test]
    } for test,train in folds_idx]


for k_fold in range(k_fold_n):

    # Build a dataset
    print "Starting k-fold", k_fold
    fold = PDB_folds[k_fold]

    # Load training dataset
    X_train,Y_train = load_fold_dataset(fold["train"])

    # Subsampling to achive even rates
    XS,YS = [],[]
    print "Building balanced subsampling"
    ##for n in range(ratio_TP_to_FP/2):
    for n in range(subsample_tasks):
        x_subsampled, y_subsampled = balanced_subsample(X_train,Y_train)
        XS.append(x_subsampled)
        YS.append(y_subsampled)
    X_train = np.vstack(XS)
    Y_train = np.hstack(YS)

    stats = {
        "k_fold_n": k_fold_n,
        "k_fold":k_fold,
        "kernel_window":kernel_window,
        "n_estimators":n_estimators,
        "clf_name":clf_name,
        "test_pdb":fold["test"].tolist(),
        "train_pdb":fold["train"].tolist(),
        "diag_window":diag_window,
        "ratio_TP_to_FP":ratio_TP_to_FP,
    }
    print fold["test"]

    print "Training ExtraTrees classifier"
    clf = Classifier(n_estimators=n_estimators,n_jobs=28,
                     #class_weight='subsample')
                     class_weight="auto") # ExtraTrees
    clf.fit(X_train,Y_train)
    stats["train_acc"] = clf.score(X_train, Y_train)

    print "Training complete"
    print 'Training Accuracy: %.3f'%stats["train_acc"]
    
    del X_train, Y_train
    gc.collect()


    X_test,Y_test = load_fold_dataset(fold["test"])

    stats["test_acc"] = clf.score(X_test, Y_test)
    print 'Testing Accuracy: %.3f'%stats["test_acc"]

    X_test_TP = X_test[Y_test==1]
    Y_test_TP = Y_test[Y_test==1]
    stats["test_acc_TP"] = clf.score(X_test_TP, Y_test_TP)
    print 'Testing Accuracy TP: %.3f'%stats["test_acc_TP"]

    X_test_FP = X_test[Y_test==0]
    Y_test_FP = Y_test[Y_test==0]
    stats["test_acc_FP"] = clf.score(X_test_FP, Y_test_FP)
    print 'Testing Accuracy FP: %.3f'%stats["test_acc_FP"]
        
    pred_probas = clf.predict_proba(X_test)[:,1]
    Y_predict = clf.predict(X_test)
    
    total_contacts = Y_test.sum()
    predicted_contacts = Y_predict[Y_test==1].sum()
    print 'Total contacts predicted %i/%i'%(predicted_contacts,total_contacts)
    
    fpr,tpr,_ = roc_curve(Y_test, pred_probas)
    stats["ROC_AUC"] = auc(fpr,tpr)
    print "ROC area under the curve", stats["ROC_AUC"]

    f_name = "{clf_name}_{n_estimators}_window_{kernel_window}_kfold_{k_fold}_ratio_{ratio_TP_to_FP}"
    f_name = f_name.format(**stats)
    
    clf_dir = os.path.join("clf2",f_name)
    clf_dir = clf_dir.format(**stats)
    os.system("rm -rf {dir}; mkdir {dir}".format(dir=clf_dir))

    # Save the model
    f_clf = os.path.join(clf_dir, "model.clf")
    joblib.dump(clf, f_clf)

    stats["f_clf"] = f_clf

    # Save the params
    f_json = os.path.join("clf2","stats_"+f_name+".json")
    with open(f_json,'w') as FOUT:
        js = json.dumps(stats,indent=2)
        FOUT.write(js)

    del X_test, Y_test, clf
    gc.collect()
