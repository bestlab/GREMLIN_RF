import glob, os, itertools
import numpy as np
import src.utils as utils

index_offset = 1

def choose_atom(res_name):
    if res_name == "GLY":
        return "CA"
    if res_name == "G":
        return "CA"
    return "CB"

FADE_line = ('AtomPair {atom1:s} {res1:d} {atom2:s} {res2:d} '
             'FADE -10 19 10 {score:0.2f} 0')

def build_constraint_text(pdb, f_prediction):
    IDX = np.loadtxt(f_prediction).astype(int)
    fasta,seq = utils.load_seq(pdb)
    constraint_text = []
    
    for i,j in IDX:

        atom_i = seq[i]
        atom_j = seq[j]
        
        fade = {"res1":i+index_offset,
                "res2":j+index_offset,
                "score":15.0}
        
        fade["atom1"] = choose_atom(atom_i)
        fade["atom2"] = choose_atom(atom_j)
        line = FADE_line.format(**fade)
        constraint_text.append(line)

    constraint_text = '\n'.join(constraint_text)
    return constraint_text

def process_pdb(pdb):
    print "Processing", pdb
    F_PREDICTIONS = glob.glob("predictions/"+pdb+"*")
    for f_pred in F_PREDICTIONS:
        text = build_constraint_text(pdb, f_pred)
        f_FADE = os.path.basename(f_pred).replace('.txt','.constraints')
        f_FADE = os.path.join("FADE",f_FADE)
        with open(f_FADE,'w') as FOUT:
            FOUT.write(text)

os.system('mkdir -p FADE')

F_PRED = glob.glob("predictions/*.txt")
PDB = set([os.path.basename(f).split('_')[0]
           for f in F_PRED])

ITR = itertools.imap(process_pdb, PDB)
for f_save in ITR: pass

