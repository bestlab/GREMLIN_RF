import os
import argparse
import itertools
import time

from glob import glob

import numpy as np
import pandas as pd
import seaborn as sns

def render_pdb(f_pdb, f_png, mol):
    
    mol.reinitialize()
    mol.viewport(600,600)
    
    mol.load(f_pdb)

    mol.set("seq_view",1)
    mol.set("seq_view_format",3)
    mol.bg_color("white")
    
    mol.show("cartoon")
    mol.hide("lines")

    mol.remove("hydrogens")
    mol.remove("solvent")

    mol.set("cartoon_fancy_helices",1)
    mol.set("ray_trace_mode",1)
    mol.set("antialias",4)
    mol.set("ray_opaque_background", 0)
    

    chains = mol.get_chains()
    colors = sns.color_palette("Set1", len(chains))
    
    for chain_name,color in zip(chains,colors):
        color_name = "color_{}".format(chain_name)
        mol.set_color(color_name,color)
        try:
            mol.color(color_name, "chain {}".format(chain_name))
        except:
            print "Can't set color for", chain_name

    time.sleep(1)
    mol.count_frames()
    mol.sync()

    print "Rendering ", f_png
    mol.ray()

    mol.png(f_png)

    mol.delete('all')
    print "Finished", f_png

    return f_png


__single_instance = False
def launch_pymol():
    if not __single_instance:
        import pymol
        pymol.finish_launching()
        pymol.cmd.sync()
        mol = pymol.cmd
        _single_instance = True
    return mol


def PDB_ITR():
    for pdb_dir in os.listdir('pdb'):
        for f_pdb in glob(os.path.join("pdb", pdb_dir, "*.pdb")):
            yield f_pdb

#data_desc = [x for x in cargs.data_directory.split('/') if x][-1]
os.system('mkdir -p png')

#png_dir = 'png/{}'.format(data_desc)
#os.system('mkdir -p {}'.format(png_dir))
#print "Saving PNG's in directory {}".format(png_dir)

import __main__
__main__.pymol_argv = [ 'pymol', '-q']
#__main__.pymol_argv = [ 'pymol']

mol = None

for f_pdb in PDB_ITR():
    name = f_pdb.split('/')[1]
    pdb  = name.split('_')[0]

    f_out = name + '_' +os.path.basename(f_pdb).replace('.pdb.pdb','.png')
    f_png = os.path.join("png",f_out)

    if not os.path.exists(f_png):
        print "Starting render", f_pdb
        mol = launch_pymol()
        render_pdb(f_pdb, f_png, mol)
    else:
        print "Skipping ", f_pdb

    #if not mol is None:
    #    mol.quit(2)
