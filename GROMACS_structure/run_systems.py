import glob, os
import itertools, multiprocessing, subprocess
import tqdm

# echo 0 | trjconv -f traj.xtc -o movie.pdb && pymol movie.pdb

_PARALLEL = 1
MP_CORES = 28
D_SYSTEMS = sorted(glob.glob("systems/*"))

def run_system(dir):
    f_traj = os.path.join(dir,"traj.xtc")
    if os.path.exists(f_traj):
        print "SKIPPING", dir
        return dir
    
    cmd = ("cd {}; mdrun -table ../../energy_table/TABLE.xvg -s topol.tpr")
    cmd = cmd.format(dir)
    try:
        with open(os.devnull, 'w') as shutup:
            subprocess.check_call(cmd,
                                  stdout=shutup,
                                  stderr=shutup, shell=True)
    except Exception as Ex:
        print "FAILED ON ", dir, Ex
        os.remove(f_traj)

    return dir

INPUT_ITR = D_SYSTEMS
ITR = itertools.imap(run_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(run_system, INPUT_ITR)

for f in tqdm.tqdm(ITR, total=len(INPUT_ITR)):
    pass

#for d in D_SYSTEMS:
#    run_system(d)
#    pass
