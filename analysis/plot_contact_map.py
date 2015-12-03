from __future__ import division
from src.utils import fixedL_cut
from src.utils import generate_matrix_IDX
import numpy as np
import pandas as pd
import itertools,os,glob
import src.utils as utils

args = {"pdb" : "1a3a"}

os.system('mkdir -p png_cmap')

kernel_window = 2

def build_cmap_from_index(N, IDX, contacts):
    C = np.zeros((N,N))

    for idx_idx in np.where(contacts)[0]:
        i,j = IDX[idx_idx]
        C[i,j] = C[j,i] = 1
    return C
    
def compute_cmap(f_rank):
    base = os.path.basename(f_rank)
    model_name = os.path.dirname(f_rank)
    pdb  = base.split('.')[0]

    func = {
        "G2" : utils.load_improved_GREMLIN_dataset,
        "APC": utils.load_GREMLIN_dataset,
    }
        
    G = func[model_name](pdb)
    NATIVE_MATRIX = utils.load_contact_map(pdb)

    N   = G.shape[0]
    IDX = list(generate_matrix_IDX(N,kernel_window))

    # Load the native contacts
    native = np.array([NATIVE_MATRIX[idx] for idx in IDX])
    g      = np.array([G[idx] for idx in IDX])

    data = {}
    CUT_IDX = range(1,5*N,20)
    
    for cut_idx in CUT_IDX:
        contacts = fixedL_cut(g,native,cut_idx,index_only=True)
        C = build_cmap_from_index(N, IDX,contacts)
        data[cut_idx/N] = C

    return data,NATIVE_MATRIX

import pylab as plt
import matplotlib.pylab as plt
import seaborn as sns


C1,NATIVE = compute_cmap("APC/{pdb}.gremlin".format(**args))
C2,NATIVE = compute_cmap("G2/{pdb}.gremlin".format(**args))

sns.set_style("white")
f, axes = plt.subplots(1, 2, figsize=(14, 8), sharex=True, sharey=True)
axes = axes.ravel()
plt.tight_layout()

for key in sorted(C1.keys()):
    origin_loc = "upper"
    
    N = C1[key].shape[0]

    axes[1].matshow(C2[key],vmin=0,vmax=1,origin=origin_loc)
    axes[0].matshow(C1[key],origin=origin_loc)

    XS,YS = np.where(NATIVE)
    axes[0].scatter(XS,YS,alpha=0.35,color='r',s=8,edgecolor=None,marker='s')
    axes[1].scatter(XS,YS,alpha=0.35,color='r',s=8,edgecolor=None,marker='s')

    args["Lcut"] = key
    axes[0].text(5,-2,"{pdb} GREMLIN, $L={Lcut:0.3f} N$".format(**args),
                 fontsize=14)
    axes[1].text(5,-2,"{pdb} Random Forest, $L={Lcut:0.3f} N$".format(**args),
                 fontsize=14)


    plt.xlim(0,N)
    plt.ylim(N,0)

    f_png = 'png_cmap/{pdb}_{Lcut:0.3f}.png'.format(**args)
    print f_png
    plt.savefig(f_png)

    axes[0].cla()
    axes[1].cla()


#plt.show()


