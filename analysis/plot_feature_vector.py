from sklearn.externals import joblib
import h5py, glob, os,json
import sqlite3, itertools, argparse
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier as Classifier
from src.utils import APC_L2, generate_matrix_IDX
from src.utils import generate_feature_vectors
import scipy.linalg
import seaborn as sns
plt = sns.plt

clf_dir = "clf"

args = {}
args["kernel_window"]  = 2
args["n_estimators"]   = 200
args["ratio_TP_to_FP"] = 20

k_fold_n = 4
kernel_window = args["kernel_window"]
n_estimators = args["n_estimators"]
ratio_TP_to_FP = args["ratio_TP_to_FP"]


clf_name = args["clf_name"] = Classifier.__name__


debug_cutoff = 10**10

bad_protein_list = ["1jbe","1hdo","1tqh","1dqg","1fk5","1gbs",
                    "2phy","1xkr","1kw4","1i1j","1kqr"]

bad_protein_list = []

#######################################################################

def load_feature_importance(k_fold=0):
    stats = args.copy()
    stats["kernel_window"] = kernel_window
    stats["k_fold"] = k_fold

    f_name = ("{clf_name}_{n_estimators}_window_{kernel_window}_" +
              "kfold_{k_fold}_ratio_{ratio_TP_to_FP}")
    f_name = f_name.format(**stats)
    
    f_model = os.path.join(clf_dir,f_name)
    f_model = f_model.format(**stats)
    f_clf = os.path.join(f_model, "model.clf")

    print "Loading model", f_clf
    assert(os.path.exists(f_clf))
    
    clf   = joblib.load(f_clf)
    
    IMP = np.array([tree.feature_importances_ for tree in
                    clf.estimators_])
    return IMP


f_kernel = "kernel_test/U.txt"

if not os.path.exists(f_kernel):

    ITR = itertools.imap(load_feature_importance, xrange(k_fold_n))
    import multiprocessing 
    MP = multiprocessing.Pool()
    ITR = MP.imap(load_feature_importance, xrange(k_fold_n))
    IMP = np.vstack(list(ITR))

    print IMP.shape

    U,s,V = scipy.linalg.svd(IMP)
    print U.shape, V.shape

    # Save the kernel
    os.system('mkdir -p kernel_test')
    np.savetxt("kernel_test/U.txt",U)
    np.savetxt("kernel_test/s.txt",s)
    np.savetxt("kernel_test/V.txt",V)


U = np.loadtxt("kernel_test/U.txt")
s = np.loadtxt("kernel_test/s.txt")
V = np.loadtxt("kernel_test/V.txt")

cut = 6

#sq_cut = int(np.sqrt(cut))
plot_Xn, plot_Yn = 2,3

f, axes = plt.subplots(plot_Xn, plot_Yn,
                       figsize=(10,6),
                       sharex=True, sharey=True)

axes = axes.ravel()

for k,(row,eigenvalue) in enumerate(zip(V[:cut,],s)):
    print k
    ax = axes[k]
    if row.sum() < 1: row *= -1
    row = row.reshape((5,5))

    print row, row.max(), row.min()

    sns.heatmap(row,ax=ax,cbar=False, vmax=0.5, vmin=-0.5)
    ax.set_title("$V_{}$, $\lambda_{}={:0.3f}$".format(k,k,eigenvalue),
                 fontsize=18)

#plt.suptitle("RF importance vector decomp")
plt.tight_layout()

os.system('mkdir -p figures/')
plt.savefig('figures/avg_features.pdf')
plt.savefig('figures/avg_features.png')

#plt.figure()
#plt.plot(s)

plt.show()
