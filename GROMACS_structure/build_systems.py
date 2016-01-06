import os, glob, itertools, random, subprocess
import tqdm

contact_strength = 4.0
total_samples = 1

_PARALLEL = 1
MP_CORES = 30

# System models are stored in here
os.system('mkdir -p systems')

# Energy table is generated here for all simulations
os.system('mkdir -p energy_table')
f_energy = 'energy_table/TABLE.xvg'

if not os.path.exists(f_energy):
    print "Building custom energy table"
    os.system('python write_tables.py 0.5 0.7 > '+f_energy)

def build_system(item):
    f,seed_n = item
    
    base = os.path.basename(f).rstrip('.txt')
    base = base + "_seed_{}".format(seed_n)
    
    os.system('mkdir -p systems/'+base)

    f_contacts = os.path.join('systems',base,'contacts.dat')
    f_sequence = os.path.join('systems',base,'sequence.dat')
    f_mdp = os.path.join('systems',base,'config.mdp')
    
    if os.path.exists(f_contacts):
        return f

    #print "Building config file", f_mdp
    with open("md_config.mdp") as FIN, open(f_mdp,'w') as FOUT:
        for line in FIN:
            if "ld-seed" in line:
                line = "ld-seed    = {}\n".format(seed_n)
            FOUT.write(line)

    #print "Building contact map", f_contacts
    with open(f) as FIN, open(f_contacts,'w') as FOUT:
        for line in FIN:
            i,j = map(int,line.split())
            s = "{:<8d} {:<8d} {:0.2f}".format(i,j, contact_strength)
            FOUT.write(s+'\n')

    #print "Building sequence file", f_sequence
    pdb = base.split('_')[0]
    f_fasta = os.path.join('fasta',pdb+'.fasta')
    with open(f_fasta) as FIN:
        FIN.readline()
        seq = FIN.readline().strip()

    with open(f_sequence,'w') as FOUT:
        FOUT.write(seq)

    #print "Building gmx file"
    cmd = ("python ./go_builder/go_builder_hack.py -k {base_dir} "
           "--seq=`cat {f_seq}` -l {f_contacts} --gmx")
    cmd = cmd.format(base_dir=os.path.join('systems',base),
                     f_seq=f_sequence,
                     f_contacts=f_contacts)

    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd, stdout=shutup,stderr=shutup, shell=True)

    #print "Generating TPR file"
    org_dir = os.getcwd()
    working_dir = os.path.join("systems",base)
    os.chdir(working_dir)
    cmd = "grompp -f config.mdp -c GO.gro -p GO_gromacs_go.top"    

    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd, stdout=shutup,stderr=shutup, shell=True)

    os.chdir(org_dir)

    return f


F_CONTACT_MAPS = sorted(glob.glob("predictions/*.txt"))
SEEDS = range(total_samples)
    
INPUT_ITR = list(itertools.product(F_CONTACT_MAPS, SEEDS))

ITR = itertools.imap(build_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(build_system, INPUT_ITR)


for f in tqdm.tqdm(ITR, total=len(INPUT_ITR)):
    pass



