import glob, subprocess, os, shutil
from scipy.spatial.distance import cdist
import Bio.PDB as bio
import numpy as np

_CB_residue_cutoff = 8.0

# If True, always return CA distance
_USE_CA_ONLY = True

def mkdir_p(d):
    try:
        os.mkdir(d)
    except OSError:
        pass

mkdir_p("cmap")
mkdir_p("cmap_coordinates")

F_PDB = sorted(glob.glob("pdb/????.pdb"))

def get_residue_pos(res, f_pdb):
    name = res.get_resname()

    if _USE_CA_ONLY:
        return res["CA"].coord
    
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
    
    return R, CMAP


for f_pdb in F_PDB:
    pdb = os.path.basename(f_pdb).split('.')[0]
    f_cmap = os.path.join("cmap",pdb+'.cmap')
    f_coord = os.path.join("cmap_coordinates",pdb+'.txt')
    
    if os.path.exists(f_cmap) and os.path.exists(f_coord):
        continue

    R, CMAP = contact_map(f_pdb)
    np.savetxt(f_coord,R)
    np.savetxt(f_cmap,CMAP,fmt="%i")

    
