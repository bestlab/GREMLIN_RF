import glob, os
import src.utils as utils
import numpy as np

# Fixed L cutoffs to determine the predictions
L_SET = [0.5, 1.0, 1.5, 3.0, 3.5]

known_models = sorted(glob.glob("G2/*.gremlin"))
PDB = [os.path.basename(f).split('.')[0] for f in known_models]
kernel_window = 2

print "Found {} known models to predict.".format(len(PDB))

os.system('mkdir -p predictions')

def score_ordering(IDX, A):
    a  = np.array([A[idx] for idx in IDX])
    order = np.argsort(a)[::-1]
    return np.array(map(list,IDX))[order]

for pdb in PDB:
    G1   = utils.load_GREMLIN_dataset(pdb)
    G2  = utils.load_improved_GREMLIN_dataset(pdb)
    NATIVE = utils.load_contact_map(pdb)

    N = NATIVE.shape[0]
    IDX = utils.generate_matrix_IDX(N,kernel_window)

    sidx1 = score_ordering(IDX, G1)
    sidx2 = score_ordering(IDX, G2)

    args = {"model" : "GREMLIN", "pdb": pdb}
    f_save = "predictions/{pdb}_{model}_{L:0.2f}.txt"

    for L in L_SET:
        args["L"] = L
        cut_idx = int(N*args["L"])
        contacts = sidx1[:cut_idx]

        upper_diag_contacts = np.array([contacts[:,1],contacts[:,0]]).T
        contacts = np.vstack([contacts, upper_diag_contacts])
        
        np.savetxt(f_save.format(**args), contacts,fmt="%d")


    args["model"] = "RF"
    for L in L_SET:
        args["L"] = L
        cut_idx = int(N*args["L"])
        contacts = sidx2[:cut_idx]

        upper_diag_contacts = np.array([contacts[:,1],contacts[:,0]]).T
        contacts = np.vstack([contacts, upper_diag_contacts])
        
        np.savetxt(f_save.format(**args), contacts,fmt="%d")

    # Save an exact model for a lower bound
    sidx3 = score_ordering(IDX, NATIVE)
    total_contacts = sum([NATIVE[i,j] for i,j in sidx3])
    contacts = sidx3[:total_contacts]

    upper_diag_contacts = np.array([contacts[:,1],contacts[:,0]]).T
    contacts = np.vstack([contacts, upper_diag_contacts])

    args["model"] = "exact"
    args["L"] = 1
    
    np.savetxt(f_save.format(**args), contacts,fmt="%d")

    print "Saved predictions for", pdb



