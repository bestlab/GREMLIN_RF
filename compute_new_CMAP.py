from sklearn.externals import joblib
import h5py, glob, os,json
import sqlite3, itertools
import numpy as np
import argparse

import pylab as plt
import matplotlib.pylab as plt
import seaborn as sns


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--f_stats",type=str,required=True)
parser.add_argument("-m", "--f_matrix", type=str,required=True)
parser.add_argument("-n", "--f_native", type=str,required=True)
args = parser.parse_args()


def APC_fro(X):
    score  = X
    row    = np.mean(score,axis=0)
    weight = np.outer(row,row.T) / score.mean()
    final = score-weight
    final = (final+final.T)/2
    return final


with open(args.f_stats) as FIN:
    js = json.load(FIN)


GO = np.loadtxt(args.f_matrix)
G = APC_fro(GO)

N = G.shape[0]
print G.mean(), GO.mean()
GO = G

NATIVE_X, NATIVE_Y = [],[]
with open(args.f_native) as FIN:
    for line in FIN:
        if "#" not in line:
            line = line.split()
            i,j = map(int,line[:2])
            NATIVE_X.append(i)
            NATIVE_Y.append(j)

            NATIVE_X.append(j)
            NATIVE_Y.append(i)

clf = joblib.load(js["f_clf"])


dx = js["kernel_window"] # This is symmetric so image size is 2*dx + 1
image_shape  = (2*dx+1,2*dx+1,20*20)

dataX,IDX = [],[]
for i,j in itertools.product(range(dx,N-dx),repeat=2):
        
    # Don't take ultra close range contacts near or below the diagonal
    if i<=j+dx: continue

    IDX.append((i,j))
    image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()
    dataX.append(image)

print "Reshaping data"
X = np.array(dataX)

print "Running prediction"
Y2 = clf.predict_proba(X)[:,1]
Y1 = clf.predict(X)
print "Done with prediction"

CMAP1 = np.zeros((N,N))
CMAP2 = np.zeros((N,N))
for k,(i,j) in enumerate(IDX):
    CMAP1[i,j] = CMAP1[j,i] = Y1[k]
    CMAP2[i,j] = CMAP2[j,i] = Y2[k]

sns.set_style("white")
f, axes = plt.subplots(2, 2, figsize=(10, 10),
                       sharex=True, sharey=True)
axes = axes.ravel()
origin_loc = "upper"

#CMAP_list = np.sort(CMAP.ravel())[::-1]

axes[3].set_title("Random Forest enchanced GREMLIN")
axes[2].set_title("Random Forest contacts")

# Hide diag_window
for i,j in itertools.product(range(dx,N-dx),repeat=2):
    # Don't take ultra close range contacts near or below the diagonal
    if abs(i-j)<=js['diag_window']:
        G[i,j] = 0


axes[0].set_title("GREMLIN contact map L cutoff")

go_list = np.sort(GO.ravel())[::-1]
cutoff = go_list[int(2*(3/2.0)*N)]
GO_32L = GO.copy()
GO_32L[GO<cutoff] = 0
GO_32L[GO>=cutoff] = 1

axes[1].set_title("GREMLIN contact map (3/2) L cutoff")

scatter_args = {"s":25, "alpha":0.80, "marker":'s', "lw":0}

TP,TN,FP,FN = [],[],[],[]
NATIVE = zip(*[NATIVE_X,NATIVE_Y])
for i,j in itertools.product(range(N),repeat=2):
    key = i,j
    if GO_32L[i,j]:
        if key in NATIVE: TP.append(key)
        if key not in NATIVE: FP.append(key)
    else:
        if key in NATIVE: FN.append(key)
        if key not in NATIVE: TN.append(key)

axes[1].scatter(zip(*TP)[0], zip(*TP)[1], color='g', **scatter_args)
axes[1].scatter(zip(*FP)[0], zip(*FP)[1], color='r', **scatter_args)
axes[1].scatter(zip(*FN)[0], zip(*FN)[1], color='m', **scatter_args)

TP,TN,FP,FN = [],[],[],[]
NATIVE = zip(*[NATIVE_X,NATIVE_Y])
for i,j in itertools.product(range(N),repeat=2):
    key = i,j
    if CMAP1[i,j]:
        if key in NATIVE: TP.append(key)
        if key not in NATIVE: FP.append(key)
    else:
        if key in NATIVE: FN.append(key)
        if key not in NATIVE: TN.append(key)

axes[3].scatter(zip(*TP)[0], zip(*TP)[1], color='g', **scatter_args)
axes[3].scatter(zip(*FP)[0], zip(*FP)[1], color='r', **scatter_args)
axes[3].scatter(zip(*FN)[0], zip(*FN)[1], color='m', **scatter_args)

#plt.tight_layout()

#ax = axes[3]
#sns.distplot(TP_FP,label="True Native Contacts",ax=ax)
#sns.distplot(TN_FN,label="Non-Native Contacts",ax=ax)
#ax.set_legend(loc="best")
#ax.set_xlabel("Random Forest classifier score")

title = "{}, w={}, TP/FP ratio={}"
plt.suptitle(title.format(args.f_matrix,
                          js["kernel_window"],
                          js["ratio_TP_to_FP"]))

np.savetxt("matrix.txt",CMAP1,fmt="%d")

#axes[3].matshow(CMAP1,vmin=0,vmax=1,origin=origin_loc)
axes[0].matshow(G,vmin=0,vmax=1,origin=origin_loc)
#axes[1].matshow(GO_32L,vmin=0,vmax=1,origin=origin_loc)
axes[2].matshow(CMAP2,origin=origin_loc)

plt.ylim(N,0)
plt.xlim(0,N)

plt.show()


