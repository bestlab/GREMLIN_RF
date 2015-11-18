import glob, subprocess, os, shutil

from scipy.spatial.distance import cdist
import Bio.PDB as bio
import numpy as np

_CB_residue_cutoff = 8.0


def mkdir_p(d):
    try:
        os.mkdir(d)
    except OSError:
        pass

mkdir_p("cmap")


F_PDB = sorted(glob.glob("pdb/????.pdb"))



def get_residue_pos(res, f_pdb):
    name = res.get_resname()
    
    # Glycine atoms use CA
    if name == "GLY":
        return res["CA"].coord

    # All other atoms use CB
    try:
        return res["CB"].coord
    except KeyError:
        print "Protein",f_pdb, "missing CB atom for", res, "using CA instead"
        return res["CA"].coord

    raise ValueError

def contact_map(f_pdb, r_cutoff=_CB_residue_cutoff):

    parser = bio.PDBParser()
    model  = parser.get_structure("model", f_pdb)
    protein = model.get_chains().next()

    R = []
    for residue in protein.get_residues():
        xyz = get_residue_pos(residue, f_pdb)
        R.append(xyz)
        
    R = np.array(R)
    DR = cdist(R,R)
    CMAP = (DR <= r_cutoff).astype(int)

    print "Completed CMAP for", f_pdb
    
    return CMAP



for f_pdb in F_PDB:
    pdb = os.path.basename(f_pdb).split('.')[0]
    f_cmap = os.path.join("cmap",pdb+'.cmap')
    
    if os.path.exists(f_cmap): continue
    
    CMAP = contact_map(f_pdb)
    
    np.savetxt(f_cmap,CMAP,fmt="%i")

    
