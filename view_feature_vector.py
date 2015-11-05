from sklearn.externals import joblib
import h5py, glob, os,json
import sqlite3, itertools
import numpy as np

kernel_window = 2
ratio_TP_to_FP = 40
#ratio_TP_to_FP = 10.0

k_fold = 0
args = {
    "kernel_window":kernel_window,
    "ratio_TP_to_FP":ratio_TP_to_FP,
    "clf_name":"stats_ExtraTreesClassifier",
    "k_fold":k_fold,
    "n_estimators":200,
}
        
f_name = "{clf_name}_{n_estimators}_window_{kernel_window}_kfold_{k_fold}_ratio_{ratio_TP_to_FP}.json"
f_json = os.path.join("clf", f_name.format(**args))
with open(f_json) as FIN:
    js = json.load(FIN)


f_clf = js["f_clf"]
clf = joblib.load(f_clf)

IMP = np.array([tree.feature_importances_ for tree in clf.estimators_])
print IMP
print IMP.shape
import scipy.linalg
U,s,V = scipy.linalg.svd(IMP)
print U.shape, V.shape

import seaborn as sns
plt = sns.plt
cut = 16
ax = sns.heatmap(V[:cut,:])

labels = ["{:0.3f}".format(eigenvalue/s[0]) for eigenvalue in s][::-1]
ax.set_yticklabels(labels)

sq_cut = int(np.sqrt(cut))
f, axes = plt.subplots(sq_cut, sq_cut, figsize=(7, 7), sharex=True, sharey=True)
axes = axes.ravel()

for k,row in enumerate(V[:cut,]):
    print k
    ax = axes[k]
    if row.sum() < 1: row *= -1
    row = row.reshape((5,5))

    sns.heatmap(row,ax=ax,cbar=False)
    ax.set_title("{}".format(k))

plt.suptitle("RF importance vector decomp")
plt.tight_layout()

plt.show()
exit()

print clf
F = clf.feature_importances_
import seaborn as sns
plt = sns.plt

fig, axes = sns.plt.subplots(1, 2, figsize=(16, 7))

dx = kernel_window*2 + 1
F = F.reshape((dx,dx))

ax0, ax1 = axes.ravel()

ax0.set_title("features for RF, window {}".format(kernel_window))
sns.heatmap(F,ax=ax0,annot=True)

ax1.set_title("log(features) for RF, window {}".format(kernel_window))
sns.heatmap(np.log(F),ax=ax1)
#plt.suptitle("Long-range contacts (>20)")
plt.suptitle("short & long range contacts (>5)")
sns.plt.tight_layout()
sns.plt.show()
