import numpy as np
import pandas as pd
import collections, os, itertools
from glob import glob
import multiprocessing
import seaborn as sns
plt = sns.plt

FLAG_plot_fill = True
reference_dir = "stats/APC"

col1,col2 = "precision","sensitivity"
#col1,col2 = "sensitivity","specificity"

xname = "precision $TP/(TP+FP)$"
yname = "sensitivity $TP/(TP+FN)$"

glob_str = "????.txt"


testing = {
#    "RF_10" : "RF baseline",
#    "RF_balanced_subtree" : "Random Forests model",
    "stats/G2":"RF model",
}


data = {}
bin_n = 30
bins = np.linspace(0,1,bin_n)

pdb_cutoff = 10**10

def figure_options():
    
    plt.legend(loc=0,fontsize=18)
    plt.xlabel(yname,fontsize=18)
    plt.ylabel(xname,fontsize=18)

    plt.xlim(0,1)
    plt.ylim(0,1)


def load_single_stat_file(f_txt):
    global col1, col2
    
    #print "Loading", f_txt
    cols = ("sensitivity","precision",
            "specificity","negative_predictive_value")

    try:
        data = np.loadtxt(f_txt)
        df = pd.read_csv(f_txt, sep=" ",index_col=0,names=cols)
    except ValueError:
        print "Problem with", f_txt
        return None

    SX = df[col1]
    SY = df[col2]

    mu_Y = SY.groupby(pd.cut(SX,bins)).agg("mean")
    mu_Y = mu_Y.interpolate()

    return mu_Y


def load_glob(F_TXT):

    YA = collections.defaultdict(list)

    #ITR = itertools.imap(load_single_stat_file, F_TXT)
    P = multiprocessing.Pool()
    ITR = P.imap(load_single_stat_file, F_TXT)
    
    for mu_Y in ITR:

        if mu_Y is None:
            continue

        for x,y in zip(bins,mu_Y):

            if not np.isnan(y):
                YA[x].append(y)

    for x in bins: YA[x]
    keys = sorted(YA.keys())

    print "Converting to numpy"
    for x in keys:
        YA[x] = np.array(YA[x])

    print "Computing means, std"
    MU_X = np.zeros(bins.size)
    MU_Y = np.zeros(bins.size)
    STD_X = np.zeros(bins.size)
    STD_Y = np.zeros(bins.size)


    for i,x in enumerate(keys):
        MU_X[i] = x
        STD_X[i] = x
        if YA[x].any(): 
            MU_Y[i] = np.mean(YA[x])
            STD_Y[i] = np.std(YA[x])
    
        else:
            MU_Y[i]  = 0
            STD_Y[i] = 0


    return MU_X, MU_Y, STD_X, STD_Y


f_dir = os.path.join(reference_dir, glob_str)
F_FIXED = sorted(glob(f_dir))[:pdb_cutoff]


MU_X, MU_Y, STD_X, STD_Y = load_glob(F_FIXED)

colors = sns.color_palette("hls", len(testing)+1)
alpha = 0.3
color = colors[0]

plt.plot(MU_X,MU_Y,label="GREMLIN",
         linestyle='--',
         color=color)

plt.fill_between(STD_X, MU_Y-STD_Y, MU_Y+STD_Y,
    alpha=alpha, facecolor=color,
    edgecolor=None,
    linewidth=0, antialiased=True)

figure_options()
plt.savefig("figures/GREMLIN_only_Acc_Pre.png")

for n,testing_dir in enumerate(testing):
    label = testing[testing_dir]    

    f_dir = os.path.join(testing_dir, glob_str)
    F_RF = sorted(glob(f_dir))[:pdb_cutoff]

    color = colors[n+1]

    MU_X, MU_Y, STD_X, STD_Y = load_glob(F_RF)
    plt.plot(MU_X,MU_Y,label=label,color=color)

    if not FLAG_plot_fill:
        continue

    plt.fill_between(STD_X, MU_Y-STD_Y, MU_Y+STD_Y,
        alpha=alpha, facecolor=color,
        linewidth=0, antialiased=True)


#######################################################

figure_options()
plt.savefig("figures/GREMLIN_RF_Acc_Pre.png")
plt.show()
#exit()
