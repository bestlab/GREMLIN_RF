import matplotlib.pylab as plt
import seaborn as sns

import collections
import numpy as np
import pandas as pd
import itertools,os,glob

F_FP = glob.glob("FP_distance/*.txt")[:]
data = np.vstack([np.loadtxt(f) for f in F_FP])

idx = np.argsort(data[:,0])
data = data[idx]

cols = ["L-cut","GREMLIN","RF"]
df = pd.DataFrame(data=data,columns=cols)

cut = pd.cut(df["L-cut"],np.linspace(0,5,50))

current_palette = itertools.cycle(sns.color_palette())

def plot_label(label,offset=0):
    Y = df.groupby(cut)[label].mean()
    YERR = df.groupby(cut)[label].std()
    X = df.groupby(cut)["L-cut"].mean()

    color = current_palette.next()
    plt.plot(X+offset,Y,label=label,lw=2,color=color)
    plt.fill_between(X+offset,Y+YERR, Y-YERR,alpha=0.15,color=color)
    return X,Y

plot_label("GREMLIN")
plot_label("RF")

plt.xlabel("cutoff",fontsize=18)
plt.xlim(0,5.0)
plt.ylabel("avg. CMAP distance to TP (True positive)",fontsize=18)

plt.legend(loc=2,fontsize=18)
plt.ylim(ymin=0)

plt.savefig("figures/FP_distance.png")
plt.savefig("figures/FP_distance.pdf")
plt.show()
