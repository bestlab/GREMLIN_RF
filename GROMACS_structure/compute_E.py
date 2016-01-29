import glob, os, subprocess
import itertools, multiprocessing
import tqdm

_PARALLEL = True
_PARALLEL = False
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):

    f_xtc = os.path.join(dir,"traj.xtc")
    if not os.path.exists(f_xtc):
        print "INCOMPLETE", dir
        return dir
    
    f_traj = os.path.join(dir,"energy.xvg")
    if os.path.exists(f_traj):
        print "Already computed", dir
        return dir
    
    cmd = ("cd {}; echo 7 | g_energy -f ener.edr -s topol.tpr")
    cmd = cmd.format(dir)
    
    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd, stdout=shutup,stderr=shutup,
                              shell=True)

    return dir

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for k,f in enumerate(tqdm.tqdm(ITR, total=len(INPUT_ITR))):
    pass


#for d in D_SYSTEMS:
#    run_system(d)
#    pass
