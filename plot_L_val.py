import h5py, itertools, os, json, gc, glob, argparse
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.externals import joblib
from src.utils import APC_L2

#kernel_window  = 2
#ratio_TP_to_FP = 1.0
#ratio_TP_to_FP = 5.0
#ratio_TP_to_FP = 10.0

diag_window = 5
debug_cutoff = 10**10

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--kernel_window",type=int,required=True)
parser.add_argument("-r", "--ratio_TP_to_FP", type=int,default=40)
args = parser.parse_args()

kernel_window = args.kernel_window
ratio_TP_to_FP = args.ratio_TP_to_FP
k_fold = 0

args = {
    "kernel_window":args.kernel_window,
    "ratio_TP_to_FP":args.ratio_TP_to_FP,
    "clf_name":"stats_ExtraTreesClassifier",
    "k_fold":k_fold,
    "n_estimators":200,
}
        

CLF_MODEL = {}
pdb_cutoff = 10**10

def load_image_data(pdb):
    global CLF_MODEL
    
    G = h5["PDB"][pdb]["gremlin_matrix"]

    key = 'weighted_scores_sum'
    W = h5["PDB"][pdb]['scores'][key][:]
    N = G.shape[0]

    data_dir = "/media/hoppeta/9cb7fc61-4955-439c-a637-52c9b6d1eafd/research/fei_research/backup_files/gremlin_ranking_scores"

    W = np.zeros((N,N),dtype=float)
    with open(os.path.join(data_dir, pdb+".ij.ranking.dat")) as FIN:
        for line in FIN:
            i,j,g = map(float,line.split())
            W[i-1][j-1] = W[j-1][i-1] = g

    # Load the contact map
    NATIVE = np.zeros(shape=(N,N))
    dx = kernel_window

    for i,j in h5["PDB"][pdb]["native_contacts"][:]:
        i,j = sorted([i,j])[::-1]
    
        # Don't take ultra close range contacts near or below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        NATIVE[i,j] = NATIVE[j,i] =1

    f_clf = PDB_F_CLF[pdb]
    
    print f_clf
    if f_clf not in CLF_MODEL:
        print "Loading model", f_clf

        # Remove an old model for memory sake
        if CLF_MODEL:
            key = CLF_MODEL.keys().pop()
            del CLF_MODEL[key]
            gc.collect()

        CLF_MODEL[f_clf] = joblib.load(f_clf)
    
    clf = CLF_MODEL[f_clf]
        
    print "Loading GREMLIN data"
    G = h5["PDB"][pdb]["gremlin_matrix"][:]
    G = APC_L2(G)
    #image_shape  = (2*dx+1,2*dx+1,20*20)

    print "Segmenting GREMLIN images"
    IDX = []; dataX = []
    for i,j in itertools.product(range(dx,N-dx),repeat=2):
        # Don't take ultra close range contacts near or below the diagonal
        if i<=j+diag_window: continue

        # Edges
        if i-dx < 0: continue
        if j-dx < 0: continue

        if i+dx >= N: continue
        if j+dx >= N: continue
        
        #image = G[i-dx:i+dx+1,j-dx:j+dx+1]
        #image = image.reshape(image_shape)

        image = G[i-dx:i+dx+1,j-dx:j+dx+1].ravel()
        IDX.append((i,j))
        dataX.append(image)

    print "Reshaping data", len(IDX)
    dataX = np.array(dataX)
    #dataX = dataX.swapaxes(1,3).swapaxes(2,3)
    #X = dataX.reshape(dataX.shape[0],(20*20)*(2*kernel_window+1)**2)
    X = dataX

    print "Running prediction"
    #Y = clf.predict_proba(X)[:,1]
    Y = clf.predict(X)

    GMATRIX = np.zeros(shape=(N,N))
    for (i,j),y in zip(IDX,Y):
        GMATRIX[i,j] = GMATRIX[j,i] = y

    #print GMATRIX[NATIVE==1].sum()/2
    #print W.shape, NATIVE.shape, GMATRIX.shape

    #sns.heatmap(GMATRIX)
    #sns.plt.figure()
    #sns.heatmap(NATIVE)
    #sns.plt.figure()
    #sns.heatmap(W)
    #sns.plt.show()
    
    # Upper diagonal
    idx = np.triu_indices(N, k=2)
    W = W[idx]
    NATIVE = NATIVE[idx]
    GMATRIX = GMATRIX[idx]

    total_native = NATIVE.sum()

    fixed_L = int((3.0/2)*N)

    # sorting by scores
    idx = np.argsort(W)[::-1]
    W = W[idx]
    NATIVEX = NATIVE.copy()[idx]

    # All L values
    FRACTION_NATIVE = np.cumsum(NATIVEX)/total_native
    ACCURACY = np.cumsum(NATIVEX)/(np.arange(len(NATIVEX))+1)

    # Fixed L prediction
    F_NATIVE = NATIVEX[:fixed_L].sum()/total_native
    F_ACC    = NATIVEX[:fixed_L].sum()/float(fixed_L)
    
    # RF prediction
    PRED = NATIVE[np.where(GMATRIX)[0]]
    P_NATIVE = PRED.sum()/total_native
    P_ACC = PRED.sum()/float(PRED.size)

    print "Total native/total predicted", NATIVE.sum(), GMATRIX.sum()
    print "RF Accuracy/P_native", P_ACC, P_NATIVE
    return FRACTION_NATIVE, ACCURACY, P_NATIVE, P_ACC, F_NATIVE, F_ACC

#############################

# Load the GREMLIN dataset 
h5 = h5py.File("data/gremlin_data.h5",'r')
PDB = h5["PDB"].keys()[:debug_cutoff]

# Remove a few bad proteins
PDB = np.array([pdb for pdb in PDB if pdb not in ["1jbe","1hdo","1tqh"]])

# For the PDB list, find out which classifier they belong to
F_CLF_JSON = glob.glob("clf/*window_{}_*_ratio_{}.json".format(kernel_window,ratio_TP_to_FP))
PDB_F_CLF = {}

for f_json in F_CLF_JSON:
    with open(f_json) as FIN:
        js = json.loads(FIN.read())

    f_clf = js["f_clf"]
    for pdb in js["test_pdb"]:
        PDB_F_CLF[pdb] = f_clf

PX, PY = [],[]
FX, FY = [],[]

for pdb in PDB[:pdb_cutoff]:
    print pdb

    #X,Y,px,py,fx,fy = load_image_data(pdb)
    print pdb,px,py,fx,fy
    
    PX.append(px)
    PY.append(py)
    FX.append(fx)
    FY.append(fy)
    
    plt.plot(X,Y,alpha=.225,color='k')


label = "RF, kernel_window {kernel_window}, ratio_TP_to_FP {ratio_TP_to_FP}"
plt.scatter(PX,PY,color='r', label=label.format(**args))
plt.legend()

plt.xlabel("fraction of native contacts")
plt.ylabel("accuracy")
plt.xlim(0,1.0)
plt.ylim(0,1)

line = np.linspace(0,1,100)
plt.figure()
plt.title("faction of native contacts; kernel_window {kernel_window}, ratio_TP_to_FP {ratio_TP_to_FP}".format(**args))
plt.scatter(FX,PX)
plt.plot(line,line,'r--',alpha=.75,lw=1)
plt.ylabel("prediction from RF")
plt.xlabel("prediction from (3/2)L cutoff")
plt.xlim(0,1)
plt.ylim(0,1)

plt.figure()
plt.title("accuracy; kernel_window {kernel_window}, ratio_TP_to_FP {ratio_TP_to_FP}".format(**args))
plt.scatter(FY,PY)
plt.plot(line,line,'r--',alpha=.75,lw=1)
plt.ylabel("prediction from RF")
plt.xlabel("prediction from (3/2)L cutoff")
plt.xlim(0,1)
plt.ylim(0,1)

plt.show()

