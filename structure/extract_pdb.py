import glob, os, subprocess
import itertools, multiprocessing
import tqdm
import numpy as np

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    org_dir = os.getcwd()
    os.chdir(dir)

    f_movie = "movie.pdb"
    f_traj  = "traj.xtc"
    f_coord = "coord.h5"
    
    if os.path.exists(f_movie):
        print "Already computed movie", dir
        os.chdir(org_dir)
        return dir

    if os.path.exists(f_coord):
        print "Already computed coord.h5", dir
        os.chdir(org_dir)
        return dir

    if not os.path.exists(f_traj):
        print "Missing traj", f_traj, dir
        os.chdir(org_dir)
        return dir

    print "Starting", dir

    cmd = ("echo 0 0 | "
           "trjconv -f {} -o {} "
           "-center -pbc mol -ur compact -conect").format(f_traj,f_movie)

    try:
        with open(os.devnull, 'w') as shutup:
            subprocess.check_call(cmd, stdout=shutup,stderr=shutup,
                                shell=True)
    except:
        print "Error with", dir
    
    os.chdir(org_dir)
    return dir


INPUT_ITR = D_SYSTEMS[:]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for item in ITR: pass

