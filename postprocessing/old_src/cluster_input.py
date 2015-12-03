import glob, os, itertools
import pandas as pd

#local_dir = os.getcwd()
# Use less cores or system will top out
CPU_CORES = 10

os.system("mkdir -p clusters")
os.system("mkdir -p clusters/clustered_outfiles")

def project_iterator():

    PDB_DIR = os.listdir("rosetta_folding_data")

    for pdb in PDB_DIR:
        FILES = glob.glob(os.path.join("rosetta_folding_data",pdb,"*.out"))
        data = []
        for f in FILES:
            tokens = f.split('_')
            model_type = tokens[-3]
            L_cut = tokens[-2]
            data.append([f,model_type,L_cut])
        df = pd.DataFrame(data=data,columns=("name","model_type","L"))

        groups = df.groupby(["model_type","L"]).groups

        for key in groups:
            idx = groups[key]
            filelist = df.ix[idx].name.values
            model_type,L = key

            yield pdb, model_type,L,filelist


def cluster_project((pdb,model_type,L,filelist)):
    
    f_name = "{}_{}_{}".format(pdb,model_type,L)
    f_out = "clusters/combined/" + f_name + '.out'
    
    cmd = ("./cluster.linuxgccrelease "
           "-in:file:silent {f_combined} "
           "-out:file:silent {f_cluster} "
           "-sort_groups_by_energy true "
           "-limit_cluster_size 20")

    f_cluster = "clusters/clustered_outfiles/" + f_name + '.out'
    cmd = cmd.format(f_cluster=f_cluster,f_combined=f_out)

    os.system(cmd)

    return cmd

INPUT_ITR = project_iterator()
ITR = itertools.imap(cluster_project, INPUT_ITR)

import multiprocessing
P = multiprocessing.Pool(CPU_CORES)
ITR = P.imap(cluster_project, INPUT_ITR)

for item in ITR:
    print item

