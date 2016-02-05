from __future__ import division
from sklearn.externals import joblib
import h5py, glob, os,json
import sqlite3, itertools
import numpy as np
from src.utils import APC_L2

pdb = '1a3a'
clf_dir = "clf2"
kernel_window = 2
diag_window = 5
ratio_TP_to_FP = 10

k_fold = 0
args = {
    "kernel_window":kernel_window,
    "ratio_TP_to_FP":ratio_TP_to_FP,
    "clf_name":"stats_ExtraTreesClassifier",
    "k_fold":k_fold,
    "n_estimators":200,
}
        
f_name = "{clf_name}_{n_estimators}_window_{kernel_window}_kfold_{k_fold}_ratio_{ratio_TP_to_FP}.json"
f_json = os.path.join(clf_dir, f_name.format(**args))

with open(f_json) as FIN:
    js = json.load(FIN)

assert( pdb in js["test_pdb"] )


f_clf = js["f_clf"]
clf   = joblib.load(f_clf)
h5 = h5py.File("data/gremlin_data.h5",'r')
h5_group = h5['PDB'][pdb]

G = h5_group["gremlin_matrix"]
N = G.shape[0]

NATIVE = np.zeros((N,N))
for k,(i,j) in enumerate(h5_group["native_contacts"]):
    NATIVE[i,j] = NATIVE[j,i] = 1

G = APC_L2(G[:])

dx = kernel_window # This is symmetric so image size is 2*dx + 1
image_shape  = (2*dx+1,2*dx+1,20*20)
total_images = (N-2*dx)**2

#print "Segmenting GREMLIN images"
IDX = []; X = []; GRANK = []
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
    GRANK.append(G[i,j])

def clean_diag(X):
    for i,j in IDX:
        X[(i,j)] += 1
        X[(j,i)] += 1
    X[X==1] = 0
    X[X>1 ] = 1

def confusion_stats(X,Y):
    X = np.triu(X,k=diag_window).ravel()
    Y = np.triu(Y,k=diag_window).ravel()

    TP = ((X==1)&(Y==1)).sum()
    TN = ((X==0)&(Y==0)).sum()
    FP = ((X==0)&(Y==1)).sum()
    FN = ((X==1)&(Y==0)).sum()


    conf = {}
    conf["sensitivity"] = TP/(TP+FN)
    conf["precision"]   = TP/(TP+FP)
    return conf

clean_diag(NATIVE)

# Find top 3/2 contacts
GRANK  = np.array(GRANK)
cutoff = GRANK[np.argsort(GRANK)[::-1]][:int((3/2)*N)][-1]
GM = np.zeros((N,N))
GM[G>=cutoff] = 1
clean_diag(GM)

# Running prediction
Y  = clf.predict(X)
Yp = clf.predict_proba(X)[:,1]

cutoff_Yp = Yp[np.argsort(Yp)[::-1]][:int((3/2)*N)][-1]
print cutoff_Yp

CMAP = np.zeros((N,N))
for k,(i,j) in enumerate(IDX):
    
    if Yp[k] > cutoff_Yp:
        CMAP[i,j] = CMAP[j,i] = 1

    #CMAP[i,j] = CMAP[j,i] = Y[k]


############## Stats
print "Fixed 3/2 L"
print confusion_stats(NATIVE, GM)
print "RF model"
print confusion_stats(NATIVE, CMAP)
##############
#exit()



import pylab as plt
import matplotlib.pylab as plt
import seaborn as sns

sns.set_style("white")
f, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)
axes = axes.ravel()

origin_loc = "upper"

CMAP_list = np.sort(CMAP.ravel())[::-1]
axes[0].matshow(CMAP,vmin=0,vmax=1,origin=origin_loc)
axes[0].set_title("Random Forest contact map")


axes[1].matshow(NATIVE,origin=origin_loc)
axes[1].set_title("Native contact map")

axes[2].matshow(GM,origin=origin_loc)
axes[2].set_title("GREMLIN 3/2 contact map")

plt.suptitle("{}".format(pdb))

plt.show()


