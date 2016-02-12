import glob, os, subprocess, h5py
import itertools, multiprocessing
import tqdm
import numpy as np
import MDAnalysis

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

# Infer the numnber of projected timesteps saved
f_config = "md_config.mdp"
with open(f_config) as FIN:
    for line in FIN:
        if "nsteps" in line:
            nsteps = int(line.split()[2])
        if "nstxtcout" in line:
            output_times = int(line.split()[2])
        
expected_output_steps = (nsteps / output_times) + 1

def run_system(dir):

    org_dir = os.getcwd()
    os.chdir(dir)

    f_movie = "movie.pdb"
    f_coord = "coord.h5"

    if not os.path.exists(f_movie):
        print "Missing movie, extract_pdb first", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_coord):
        print "Coord file exists, skipping", dir
        os.chdir(org_dir)
        return dir

    P = MDAnalysis.Universe(f_movie, multiframe=True)
    T = P.trajectory

    C = []
    for k,x in enumerate(T):
        y = np.array(x.positions)
        C.append(y)
    C = np.array(C)

    if C.shape[0] != expected_output_steps:
        print "WARNING timestep mismatch! Exiting", dir
        print C.shape[0], expected_output_steps
        os.chdir(org_dir)
        return dir


    test_val = ((C[0][0] - C[1][0])**2).sum()
    #print "Test val", test_val
    if not test_val:
        print "TEST_VAL problem (should not be zero!)", test_val
        os.chdir(org_dir)
        return dir

    ##################################################################
        
    h5 = h5py.File(f_coord,'w')
    h5["coord"] = C
    h5.close()

    # Remove the movie file
    os.system("rm movie.pdb")
    
    os.chdir(org_dir)
    return dir


INPUT_ITR = D_SYSTEMS[:]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for item in tqdm.tqdm(ITR):
    #print "Completed", item
    pass

