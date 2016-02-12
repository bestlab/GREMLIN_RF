import glob, os, subprocess
import itertools, multiprocessing
import numpy as np
import pandas as pd
import seaborn as sns
plt = sns.plt

y_name = "Q"
#y_name = "rg"

print "Loading CSV"
f_csv = os.path.join('collated_data','system_data.csv')
df_master = pd.read_csv(f_csv)

print "Loading averages"
data = []
for key,df in df_master.groupby(['pdb','system','cutoff']):
    item = {"pdb":key[0],
            "system":key[1],
            "cutoff":key[2],}

    # Identify the seed with minimum energy
    U_avg = df.groupby('seed').U.mean()
    idx_U_min_seed = U_avg.idxmin()

    # Select the seed with minimum
    df = df[df.seed == idx_U_min_seed]

    item["Y"] = df[y_name].mean()
    item["Y_err"] = df[y_name].std()

    data.append(item)


df = pd.DataFrame(data)

def plot_label(label,offset=0):
    g = df[df.system==label].groupby("cutoff")
    X,Y,Y_std = [],[],[]
    
    for x, data in g:
        X.append(x)
        Y.append(data.Y.mean())
        Y_std.append(data.Y_err.mean())

    X = np.array(X)
        
    #X, data = zip(*g)
    #X = np.array(X)
    #Y = [x[y_name].mean() for x in data]
    #Y_std = [x.Q.std() for x in data]
    #plt.errorbar(X+offset,Y,yerr=Y_std,label=label)
    
    plt.errorbar(X+offset,Y,label=label)
    plt.scatter(X+offset,Y,color='k',alpha=0.5)
    return X,Y

plot_label("GREMLIN")
plot_label("RF",offset=0.01)

NATIVE_X = np.linspace(0.1, 5, 100)
NATIVE_Y = df_master[df_master.system=="exact"][y_name].mean()
plt.plot(NATIVE_X,[NATIVE_Y,]*100,'k--',alpha=0.5,label="native")

plt.legend(loc="best",fontsize=18)
plt.xlabel("cutoff",fontsize=18)
plt.xlim(0,5.1)
plt.ylabel("{} for all proteins".format(y_name), fontsize=18)

os.system("mkdir -p figures")
plt.savefig("figures/Q_avg.png")
plt.show()
