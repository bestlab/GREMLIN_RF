import h5py, itertools, os, json, gc
import numpy as np
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import ExtraTreesClassifier as Classifier
import argparse
from src.subsampling import balanced_subsample

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-n", "--n_estimators", type=int,default=200)
#parser.add_argument("-r", "--ratio_TP_to_FP", type=float,default=5.)0
parser.add_argument("-r", "--ratio_TP_to_FP", type=int,default=2)
args = parser.parse_args()

k_fold_n = 4
kernel_window = args.kernel_window
n_estimators = args.n_estimators
clf_name =  Classifier.__name__

ratio_TP_to_FP = args.ratio_TP_to_FP

debug_cutoff = 10**10
diag_window = 5


#subsample_tasks = args.ratio_TP_to_FP
subsample_tasks = 1

#######################################################################

def APC_correction_factor(score):
    row    = np.mean(score,axis=0)
    weight = np.outer(row,row.T) / score.mean()
    return weight

def APC_fro(score):
    weight = APC_correction_factor(score)
    final = score-weight

    # Symmetrize the matrix
    final = (final+final.T)/2
    return final

def load_image_data(pdb,load_all=False):
    G = h5["PDB"][pdb]["gremlin_matrix"]
    N = G.shape[0]
    
    dx = kernel_window # This is symmetric so image size is 2*dx + 1

    # Load the contact map
    NATIVE = np.zeros(shape=(N,N))
    idx_true_pos = set()
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
        idx_true_pos.add((i,j))

    idx_true_neg = set()
    # Find the set of all possible neg_idx
    for i,j in itertools.product(range(N),repeat=2):

        # Don't take ultra close range contacts near or
        # below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue

        # Don't select native contacts
        if NATIVE[i,j]: continue

        idx_true_neg.add((i,j))

    idx_true_neg = list(idx_true_neg)
    np.random.shuffle(idx_true_neg)

    if load_all:
        FP_choosen = len(idx_true_neg)
    else:
        FP_choosen = int(ratio_TP_to_FP*len(idx_true_pos))
        
    ratio = float(len(idx_true_neg))/ len(idx_true_pos)
    
    print "{} {:5d} {:5d} {:0.4f}".format(pdb, len(idx_true_pos), FP_choosen, ratio)
    
    idx_true_neg = idx_true_neg[:FP_choosen]

    X,Y = [],[]
    G = G[:]
    G = np.linalg.norm(G,ord="fro",axis=(2,3))
    G = APC_fro(G)

    dx = kernel_window # This is symmetric so image size is 2*dx + 1
    image_shape  = (2*dx+1,2*dx+1)
    
    for i,j in idx_true_pos:
        # This loads the block of size (2*dx+1,2*dx+1,20,20)
        # Unravel the 20x20 residue channel
        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()

        X.append(image)
        Y.append(1)

    for i,j in idx_true_neg:
        # This loads the block of size (2*dx+1,2*dx+1,20,20)
        # Unravel the 20x20 residue channel
        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()
        
        X.append(image)
        Y.append(0)

    del G
    gc.collect()

    # Swap axes so the dim is (N samples, K channels, width, height)
    sample_n = len(Y)
    X = np.array(X)
    
    return X,Y

def load_fold_dataset(PDB_LIST, load_all=False):
    print "Loading {} protein GREMLIN data (loadall={})".format(len(PDB_LIST),load_all)
    X,Y = [], []
    for pdb in PDB_LIST:    
        x,y = load_image_data(pdb,load_all)
        X.append(x)
        Y.append(y)

    return np.vstack(X), np.hstack(Y)


# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

# Remove a few bad proteins
PDB = np.array([pdb for pdb in PDB if pdb not in ["1jbe","1hdo","1tqh"]])

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
    #XS,YS = [],[]
    #print "Building balanced subsampling"
    #for n in range(subsample_tasks):
    #    x_subsampled, y_subsampled = balanced_subsample(X_train,Y_train)
    #    XS.append(x_subsampled)
    #    YS.append(y_subsampled)
    #X_train = np.vstack(XS)
    #Y_train = np.hstack(YS)

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
    clf = Classifier(n_estimators=n_estimators,n_jobs=28,)
                     #class_weight='subsample')
                     #class_weight="auto") # ExtraTrees
    clf.fit(X_train,Y_train)
    stats["train_acc"] = clf.score(X_train, Y_train)

    print "Training complete"
    print 'Training Accuracy: %.3f'%stats["train_acc"]
    
    del X_train, Y_train
    gc.collect()

    # For testing, now load the entire dataset!
    X_test,Y_test = load_fold_dataset(fold["test"],load_all=True)

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
    
    #pred_prob = clf.predict(X_test)
    #import seaborn as sns
    #sns.distplot(pred_prob)
    #sns.plt.show()
    #print pred_probas
    #exit()
    
    fpr,tpr,_ = roc_curve(Y_test, pred_probas)
    stats["ROC_AUC"] = auc(fpr,tpr)
    print "ROC area under the curve", stats["ROC_AUC"]

    f_name = "{clf_name}_{n_estimators}_window_{kernel_window}_kfold_{k_fold}_ratio_{ratio_TP_to_FP}"
    f_name = f_name.format(**stats)
    clf_dir = os.path.join("clf",f_name)
    clf_dir = clf_dir.format(**stats)
    os.system("rm -rf {dir}; mkdir {dir}".format(dir=clf_dir))

    # Save the model
    f_clf = os.path.join(clf_dir, "model.clf")
    joblib.dump(clf, f_clf)

    stats["f_clf"] = f_clf

    # Save the params
    f_json = os.path.join("clf","stats_"+f_name+".json")
    with open(f_json,'w') as FOUT:
        js = json.dumps(stats,indent=2)
        FOUT.write(js)

    del X_test, Y_test, clf
    gc.collect()
