from src.utils import fixedL_cut
from src.utils import generate_matrix_IDX
import numpy as np
import pandas as pd
import itertools,os,glob
import src.utils as utils

kernel_window = 2

#############################

def compute_rank(f_rank):
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
    IDX = generate_matrix_IDX(N,kernel_window)

    # Load the native contacts
    native = np.array([NATIVE_MATRIX[idx] for idx in IDX])
    g      = np.array([G[idx] for idx in IDX])

    # sorting by scores
    cols = ("sensitivity","precision",
            "specificity","negative_predictive_value",
            "false_positive_rate", "true_positive_rate")

    print "Sorting score", pdb, model_name

    CUT_IDX = range(1,len(g))
    
    data = []
    for cut_idx in CUT_IDX:  
        item = fixedL_cut(g,native,cut_idx)
        data.append([item[k] for k in cols])
    df = pd.DataFrame(columns=cols,index=CUT_IDX,data=data)

    return pdb, model_name, df

def model_iterator(model_names):
    os.system('mkdir -p stats')

    for name in model_names:
        os.system('mkdir -p {}'.format(os.path.join("stats",name)))
        FILES = glob.glob(os.path.join(name,"*.gremlin"))
        for f in sorted(FILES):
            yield f

if __name__ == "__main__":

    MODELS = ["G2","APC"]
    func = compute_rank

    import multiprocessing
    P = multiprocessing.Pool()

    ITR = itertools.imap(func, model_iterator(MODELS))
    ITR = P.imap(func, model_iterator(MODELS))

    for item in ITR:

        pdb, model_name, df = item

        line = ("{} {sensitivity:0.8f} {precision:0.8f} " +
                "{specificity:0.8f} {negative_predictive_value:0.8f} "
                "{false_positive_rate:0.8f} {true_positive_rate:0.8f}\n")

        f_save = os.path.join('stats',model_name,pdb+'.txt')
        with open(f_save,'w') as FOUT:
            for row in df.T:
                values = df.ix[row]
                s = line.format(row,**values)
                FOUT.write(s)

        print "Completed", pdb, model_name
