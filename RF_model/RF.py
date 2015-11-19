import itertools, os, json, gc, multiprocessing, glob
import numpy as np
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import ExtraTreesClassifier as Classifier
import argparse
from src.utils import APC_L2, generate_matrix_IDX
from src.utils import generate_feature_vectors

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-n", "--n_estimators", type=int,required=True)
parser.add_argument("-r", "--ratio_TP_to_TN", type=int,required=True)
parser.add_argument("-k", "--k_folds", type=int,default=4)
args = parser.parse_args()

k_fold_n = args.k_folds
kernel_window = args.kernel_window
n_estimators = args.n_estimators
clf_name =  Classifier.__name__
ratio_TP_to_TN = args.ratio_TP_to_TN

nodesize = 2
debug_cutoff = 10**10
#debug_cutoff = 10

MP_CORES = 25

#######################################################################

# GREMLIN'S sequence encoding
GREMLIN_seq_encoding = 'ARNDCQEGHILKMFPSTWYV'
seq_encoding = dict(zip(GREMLIN_seq_encoding,range(20)))

def load_GREMLIN_dataset(pdb):
    f_GREMLIN = os.path.join("APC",pdb+'.gremlin')
    return np.loadtxt(f_GREMLIN)

def load_seq(pdb):
    f_fasta = os.path.join("fasta",pdb+'.fasta')
    with open(f_fasta) as FIN:
        FIN.readline()
        return FIN.readline().strip()

def load_contact_map(pdb):
    f_cmap = os.path.join("cmap",pdb+'.cmap')
    return np.loadtxt(f_cmap).astype(int)

def load_all_image_data(pdb):
    return load_image_data(pdb,load_all=True)

def load_image_data(pdb,load_all=False):

    g   = load_GREMLIN_dataset(pdb)
    fasta = load_seq(pdb)
    seq = [seq_encoding[aa] for aa in fasta]
    NATIVE = load_contact_map(pdb)

    # Sanity checks    
    assert(len(seq) == g.shape[0])
    assert(len(seq) == NATIVE.shape[0])

    # APC corrections and whiten
    g = APC_L2(g)

    N = g.shape[0]
    IDX = generate_matrix_IDX(N,kernel_window)

    idx_true_pos = set()
    idx_true_neg = set()

    # Load the set of TP and TN
    for i,j in IDX:
        if NATIVE[i,j]:
            idx_true_pos.add((i,j))
        else:
            idx_true_neg.add((i,j))

    # Shuffle the contacts
    idx_true_neg = list(idx_true_neg)
    np.random.shuffle(idx_true_neg)

    # If we are only loading a subset of TN, truncate here
    if load_all:
        FP_choosen = len(idx_true_neg)
    else:
        FP_choosen = int(ratio_TP_to_TN*len(idx_true_pos))
        
    ratio = float(len(idx_true_neg))/ len(idx_true_pos)

    status_str = "{} {:5d} {:5d} {:0.4f}"
    print status_str.format(pdb, len(idx_true_pos), FP_choosen, ratio)
    
    idx_true_neg = idx_true_neg[:FP_choosen]

    X0 = generate_feature_vectors(g,seq,idx_true_pos,kernel_window)
    Y0 = [1,]*len(X0)
    
    X1 = generate_feature_vectors(g,seq,idx_true_neg,kernel_window)
    Y1 = [0,]*len(X1)

    # Concatenate the two samples and make them a numpy array    
    X = np.array(X0+X1)
    Y = np.array(Y0+Y1)

    return X,Y

def load_fold_dataset(PDB_LIST, load_all=False):
    status_msg = "Loading {} protein GREMLIN data (loadall={})"
    print status_msg.format(len(PDB_LIST),load_all)
    
    X,Y = [], []

    func = load_image_data
    if load_all: func = load_all_image_data

    #ITR = itertools.imap(func, PDB_LIST)
    ITR = P.imap(func, PDB_LIST)
    
    for x,y in ITR:
        X.append(x)
        Y.append(y)

    return np.vstack(X), np.hstack(Y)

def train_model(stats, X_train, Y_train, X_test=None, Y_test=None):
        
    print "Training ExtraTrees classifier"
    clf = Classifier(n_estimators=n_estimators,n_jobs=30,
                     min_samples_leaf=nodesize,
                     #class_weight='balanced_subsample',
                     )
    clf.fit(X_train,Y_train)
    stats["train_acc"] = clf.score(X_train, Y_train)

    print "Training complete"
    print 'Training Accuracy: %.3f'%stats["train_acc"]
    
    # Breakout early if no test set is given
    if X_test is None:
        return clf, stats

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

    return clf, stats

def generate_clf_f_name(stats):
    stats["f_name"] = ("{clf_name}_{n_estimators}_window_{kernel_window}_" +
                       "kfold_{k_fold}_ratio_{ratio_TP_to_TN}").format(**stats)
    stats["clf_dir"] = os.path.join("clf",stats["f_name"]).format(**stats)
    stats["f_clf"] = os.path.join(stats["clf_dir"], "model.clf")
    return stats

def save_model(stats,clf):
    
    # Save the model
    os.system("rm -rf {clf_dir}; mkdir {clf_dir}".format(**stats))
    joblib.dump(clf, stats["f_clf"])

    # Save the params
    f_json = os.path.join("clf","stats_"+stats["f_name"]+".json")
    with open(f_json,'w') as FOUT:
        js = json.dumps(stats,indent=2)
        FOUT.write(js)



# Load the GREMLIN data files
F_GREMLIN = glob.glob("APC/*.gremlin")
PDB = np.array([os.path.basename(f).split('.')[0] for f in F_GREMLIN])

# If debugging, truncate the number of samples
PDB = PDB[:debug_cutoff]

# Multiprocessing Pool
P = multiprocessing.Pool(MP_CORES)

# Remove a few bad proteins (removed for now)
# PDB = np.array([pdb for pdb in PDB if pdb not in bad_protein_list])

# Partition the proteins into cross-validation sets
folds_idx = list(KFold(len(PDB), n_folds = k_fold_n))

# Turn the folds into PDB lists
PDB_folds = [{
    "test"  :PDB[train],
    "train" :PDB[test]
    } for test,train in folds_idx]


os.system('mkdir -p clf')

for k_fold in range(k_fold_n):

    # Build a dataset
    print "Starting k-fold", k_fold
    fold = PDB_folds[k_fold]

    stats = {
        "k_fold_n": k_fold_n,
        "k_fold"  : k_fold,
        "kernel_window":kernel_window,
        "n_estimators":n_estimators,
        "clf_name":clf_name,
        "test_pdb":fold["test"].tolist(),
        "train_pdb":fold["train"].tolist(),
        "ratio_TP_to_TN":ratio_TP_to_TN,
        "nodesize" : nodesize,
    }

    stats = generate_clf_f_name(stats)

    if os.path.exists(stats["f_clf"]):
        print stats["f_clf"], "already computed! Skipping"
        continue

    # Load training dataset
    X_train,Y_train = load_fold_dataset(fold["train"])

    # For testing, now load the entire dataset!
    X_test,Y_test = load_fold_dataset(fold["test"],load_all=True)

    clf, stats = train_model(stats, X_train, Y_train, X_test, Y_test)
    save_model(stats,clf)

    del clf
    gc.collect()

############################################################################\
# Load a best-fit model using all proteins
############################################################################

print "Starting best-fit model"

stats = {
    "k_fold_n": 1,
    "k_fold"  : 0,
    "kernel_window":kernel_window,
    "n_estimators":n_estimators,
    "clf_name":"GlobalFit",
    "test_pdb":[],
    "train_pdb":PDB.tolist(),
    "ratio_TP_to_TN":ratio_TP_to_TN,
    "nodesize" : nodesize,
}

stats = generate_clf_f_name(stats)

if os.path.exists(stats["f_clf"]):
    print stats["f_clf"], "already computed! Finished!"
    exit()

# Load training dataset
X_train,Y_train = load_fold_dataset(PDB)

clf, stats = train_model(stats, X_train, Y_train)

save_model(stats,clf)
