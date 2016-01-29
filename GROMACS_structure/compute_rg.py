import glob, os, subprocess
import itertools, multiprocessing
import tqdm

_PARALLEL = True
MP_CORES = 30
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    f_traj = os.path.join(dir,"gyrate.xvg")
    if os.path.exists(f_traj):
        print "Already computed", dir
        return dir
    
    cmd = ("cd {}; echo 0 | g_gyrate -f traj.xtc -s topol.tpr")
    cmd = cmd.format(dir)
    
    with open(os.devnull, 'w') as shutup:
        try:
            subprocess.check_call(cmd, stdout=shutup,stderr=shutup,
                                shell=True)
        except:
            print "ERROR", dir
            return dir

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
