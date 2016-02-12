import os, glob, itertools, random, subprocess
import numpy as np
import tqdm
from Bio.PDB.Polypeptide import one_to_three
from Bio.PDB import *

lower_U_rg = 0.4
upper_U_rg = 0.6

contact_strength = 4.0
contact_cutoff_scalepoint = 2.0

total_samples = 10

_USE_EXACT_ONLY = False
_USE_LIMITED_SUBSET = 0

_PARALLEL = True
MP_CORES = 30

# System models are stored in here
os.system('mkdir -p systems')

# Energy table is generated here for all simulations
os.system('mkdir -p energy_table')
f_energy = 'energy_table/TABLE.xvg'

#if not os.path.exists(f_energy):

# Always build a new energy table just in case
if True:
    print "Building custom energy table"
    cmd = 'python write_tables.py {} {} > {}'
    cmd = cmd.format(lower_U_rg,upper_U_rg,f_energy)
    os.system(cmd)

def fix_PDB(seq, f_pdb):
    #print "Fixing CA atom labels", f_pdb    
    seq3 = map(one_to_three, seq)
    SEQ_ITR = iter(seq3)

    output = []
    with open(f_pdb) as FIN:
        for line in FIN:

            if "ALA" in line:
                line = line.replace("ALA", SEQ_ITR.next())
            output.append(line)

    with open(f_pdb,'w') as FOUT:
        FOUT.write(''.join(output))

def check_PDB(f_pdb):
    parser1 = PDBParser()
    model1  = parser1.get_structure('A',f_pdb)

    X = []
    for res in model1.get_residues():
        try:
            atm = res["CA"].coord
        except Exception as Ex:
            print "Exception with CA atom", Ex, f_pdb
            return False

        X.append(atm)
    X = np.array(X)
    dist = [np.linalg.norm(x0-x1) for x0,x1 in zip(X,X[1:])]
    if max(dist) > 4.5:
        print "Gapped protein, skipping", f_pdb
        return False

    return True    

def match_PDB_coords(f_pdb1,f_pdb2):
    parser1 = PDBParser()
    model1  = parser1.get_structure('A',f_pdb1)

    parser2 = PDBParser()
    model2 = parser2.get_structure('A',f_pdb2)

    for res1,res2 in zip(model1.get_residues(),
                         model2.get_residues()):

        assert( res1.resname == res2.resname )
        res1["CA"].coord = res2["CA"].coord

    # Center the protein
    CA = np.array([x.coord for x in model1.get_atoms()])
    CA -= CA.mean(axis=0)
    for ca,atom1 in zip(CA,model1.get_atoms()):
        atom1.coord = ca

    io = PDBIO()
    io.set_structure(model1)
    io.save(f_pdb1)

def build_system(item):
    f,seed_n = item
    
    base = os.path.basename(f).rstrip('.txt')
    base = base + "_seed_{}".format(seed_n)

    pdb = base.split('_')[0]
    f_pdb = os.path.join('pdb',pdb+'.pdb')
    
    f_contacts = os.path.join('systems',base,'contacts.dat')
    f_sequence = os.path.join('systems',base,'sequence.dat')
    f_mdp = os.path.join('systems',base,'config.mdp')

    local_contact_strength  = contact_strength
    Lcut = float(base.split('_')[2])
    
    local_contact_strength *= contact_cutoff_scalepoint/Lcut


    if os.path.exists(f_mdp):
        print "System already created! Skipping", f
        return f

    if not check_PDB(f_pdb):
        return f

    os.system('mkdir -p systems/'+base)
    
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
            s = "{:<8d} {:<8d} {:0.2f}".format(i,j,
                                               local_contact_strength)
            FOUT.write(s+'\n')

    #print "Building sequence file", f_sequence
    
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
        subprocess.check_call(cmd, stdout=shutup,
                              stderr=shutup, shell=True)

    #print "Fixing input PDB"
    fix_PDB(seq, os.path.join('systems',base,'GO.pdb'))

    # Uncomment if you want to start with native structure
    #match_PDB_coords(os.path.join('systems',base,'GO.pdb'),
    #                 os.path.join('pdb',pdb+'.pdb'))

    #print "Generating coordinate file"
    org_dir = os.getcwd()
    working_dir = os.path.join("systems",base)
    os.chdir(working_dir)
    
    cmd = "editconf -f {f_pdb} -o {f_gro} -box 80 80 80"
    cmd = cmd.format(f_pdb="GO.pdb", f_gro="GO.gro")
    
    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd,shell=True,
                              stdout=shutup,stderr=shutup)
        
    #print "Generating TPR file"
    cmd = "grompp -f config.mdp -c GO.gro -p GO_gromacs_go.top -maxwarn 2" 

    with open(os.devnull, 'w') as shutup:
        subprocess.check_call(cmd, stdout=shutup,stderr=shutup, shell=True)

    os.chdir(org_dir)

    return f


F_CONTACT_MAPS = sorted(glob.glob("predictions/*.txt"))
SEEDS = range(total_samples)

if _USE_EXACT_ONLY:
    F_CONTACT_MAPS = [x for x in F_CONTACT_MAPS if "_exact" in x]

if _USE_LIMITED_SUBSET:
    F_CONTACT_MAPS = F_CONTACT_MAPS[:_USE_LIMITED_SUBSET]

    
INPUT_ITR = list(itertools.product(F_CONTACT_MAPS, SEEDS))
ITR = itertools.imap(build_system, INPUT_ITR)

if _PARALLEL:
    import multiprocessing
    MP = multiprocessing.Pool(MP_CORES)
    ITR = MP.imap(build_system, INPUT_ITR)


for f in tqdm.tqdm(ITR, total=len(INPUT_ITR)):
    pass

#for f in ITR:
#    pass



