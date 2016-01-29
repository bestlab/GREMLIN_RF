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
    
    #if os.path.exists(f_movie):
    #    print "Already computed", dir
    #    os.chdir(org_dir)
    #    return dir

    print "Starting", dir

    cmd = ("echo 0 0 | "
           "trjconv -f traj.xtc -o movie.pdb "
           "-center -pbc mol -ur compact")

    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd, stdout=shutup,stderr=shutup,
                              shell=True)
    
    os.chdir(org_dir)
    return dir


INPUT_ITR = D_SYSTEMS[:]
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for item in ITR: pass

