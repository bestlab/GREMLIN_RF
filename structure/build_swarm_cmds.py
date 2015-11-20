import glob, os, json

# We only want to test small proteins
max_protein_length = 82

total_decoy_n = 20000
decoy_n_per   = 250

F_PRED = glob.glob("FADE/*.constraints")

def get_model_name(f):
    name = os.path.basename(f).split('_')[1:]
    a,b = name[0], '.'.join(name[1].split('.')[:-1])
    model = '_'.join([a,b])
    return model

# Determine the types of models
models = set([get_model_name(f) for f in F_PRED])

PDB = list(set([os.path.basename(f).split('_')[0]
            for f in F_PRED]))

# Measure the size of each protein
def load_fasta(pdb):
    f_fasta = os.path.join("fasta",pdb+'.fasta')
    with open(f_fasta) as FIN:
        FIN.readline()
        fasta = FIN.readline().strip()
    return fasta

p_length = dict(zip(PDB,[len(load_fasta(pdb)) for pdb in PDB]))

# Remove large proteins
PDB = [pdb for pdb in PDB if p_length[pdb] <= max_protein_length]

print "Building commands for {} proteins".format(len(PDB))

os.system('mkdir -p rosetta_folding_data')

'''
print "Checking fragment generation ..."
F_FASTA = [x for x in glob.glob("fasta/*.fasta")
           if os.path.basename(x).split('.')[0] in target_proteins]

PDB_FRAGMENTS = []

for f in F_FASTA:
    pdb = os.path.basename(f).split('.')[0]
    local_dir = os.getcwd()
    new_dir   = os.path.join(local_dir, "fragments",pdb)

    check_file = os.path.join(new_dir,"t001_.make_fragments.success")
    if not os.path.exists(check_file):
        print " + PDB {} fragment information not complete".format(pdb)
    else:
        PDB_FRAGMENTS.append((pdb, new_dir))
'''

print "Allocating locations"
workdir = os.path.join(os.getcwd(),"rosetta_folding_data")


SWARM_CMDS = []
print "Building swarm file"

constrain_dir = os.path.join(os.getcwd(),"FADE")

for simname in models:

    for pdb in PDB:
        frag_dir = os.path.join("fragments",pdb)

        f_3mer = os.path.join(frag_dir, "t001_.200.3mers")
        f_9mer = os.path.join(frag_dir,"t001_.200.9mers")
        f_fasta= os.path.join(frag_dir,"t001_.fasta")
        f_psipred = os.path.join(frag_dir,"t001_.psipred_ss2")

        f_constraints = pdb + "_" + simname  + '.constraints'
        f_constraints = os.path.join(constrain_dir,f_constraints)

        try:
            assert(os.path.exists(f_constraints))
        except:
            logging.warning("Missing constraints for {}".format(pdb))
            continue

        working_dir = os.path.join(os.getcwd(), "rosetta_folding_data", pdb)

        basecmd = (
            "AbinitioRelax "
            "-in:file:fasta {} "
            "-in:file:frag3 {} "
            "-in:file:frag9 {} "
            "-use_filters true "
            "-psipred_ss2 {} "
            "-seed_offset {} "
            "-out:nstruct {} "
            "-out:file:silent {} "
            "-out:sf {} "
            "-constraints:cst_file {}"
            )

        if not os.path.isdir(working_dir):
            os.system("mkdir {}".format(working_dir))
            print "Created directory", working_dir

        for decoy_n in xrange(0, total_decoy_n, decoy_n_per):
            seed_n = decoy_n // decoy_n_per

            f_outfile = os.path.join(working_dir,
                                     "default_{}_{}.out".format(simname,seed_n))
            f_outfile_sc = os.path.join(working_dir,
                                        "score_{}_{}.fsc".format(simname,seed_n))


            transform_vars = [f_fasta, f_3mer, f_9mer, f_9mer, f_outfile,
                              f_outfile_sc, f_constraints]

            # Transform each of the files into a abspath
            f_fasta = os.path.abspath(f_fasta)
            f_3mer  = os.path.abspath(f_3mer)
            f_9mer  = os.path.abspath(f_9mer)
            f_outfile = os.path.abspath(f_outfile)
            f_outfile_sc = os.path.abspath(f_outfile_sc)
            f_constraints = os.path.abspath(f_constraints)
            

            cmd = basecmd.format(f_fasta,f_3mer,f_9mer,f_psipred, seed_n, 
                                 decoy_n_per, f_outfile, f_outfile_sc,
                                 f_constraints)

            if not os.path.exists(f_outfile):
                SWARM_CMDS.append((working_dir, cmd,))

print "Found {} AbinitioRelax commands to run".format(len(SWARM_CMDS))
swarm_template = "cd {} && module load rosetta && time {}\n"

f_swarm = "swarm_folding.sh"
print "Building {} now".format(f_swarm)

with open(f_swarm, 'w') as FOUT:
    for f_dir,f in SWARM_CMDS:
        cmd = swarm_template.format(f_dir,f)
        FOUT.write(cmd)

