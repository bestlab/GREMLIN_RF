import glob, os

### COMMAND TO RUN!
### GREMLIN/run_gremlin.sh /opt/mcr/v81/bin  input/1a3a.aln



docker_image_name = "matlab-mcl-r2013a"

F_ALN = sorted(glob.glob("input_GREMLIN/*.aln"))

F_ALN = [x for x in F_ALN if "1guu" if x]
print F_ALN

output_dir = 'output_GREMLIN'
os.system('mkdir -p {}'.format(output_dir))

docker_start = ("docker run -it "
                "--volume={input_dir}:/input "
                "--volume={output_dir}:/output "
                "--volume={gremlin_dir}:/GREMLIN "
                "--name matlab_gremlin_{pdb} "
                ) + docker_image_name + ' '

def run_aln(f_aln):
    pdb = os.path.basename(f_aln).split('.')[0]
    f_out = "{}/{}/{}.GREMLIN.txt".format(output_dir,pdb,pdb)
    os.system('mkdir -p {}/{}'.format(output_dir,pdb))

    args = {
        "pdb" : pdb,
        "input_dir" : os.path.abspath("input_GREMLIN"),
        "output_dir" : os.path.abspath("output_GREMLIN"),
        "gremlin_dir" : os.path.abspath("GREMLIN"),
    }

    if os.path.exists(f_out):
        return f_aln

    print "Starting", f_aln
    cmd = (docker_start +
           "GREMLIN/run_gremlin.sh /opt/mcr/v81 input/{pdb}.aln")
    cmd = cmd.format(**args)
    
    os.system(cmd)
    return f_aln


import itertools
ITR = itertools.imap(run_aln, F_ALN)

#import multiprocessing
#P = multiprocessing.Pool()
#ITR = P.imap(run_aln, F_ALN)

for f in ITR:
    print "Completed", f


