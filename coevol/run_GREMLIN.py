import glob, os

docker_image_name = "matlab-mcl-r2013a"
F_ALN = sorted(glob.glob("input_GREMLIN/*.aln"))

output_dir = 'output_GREMLIN'
os.system('mkdir -p {}'.format(output_dir))

docker_start = ("docker run --rm "
                "--volume={input_dir}:/input "
                "--volume={output_dir}:/output "
                "--volume={gremlin_dir}:/GREMLIN "
                "--name matlab_gremlin_{pdb} "
                ) + docker_image_name + ' '

def run_aln(f_aln):
    pdb = os.path.basename(f_aln).split('.')[0]
    os.system('mkdir -p {}/{}'.format(output_dir,pdb))

    args = {
        "pdb" : pdb,
        "input_dir"   : os.path.abspath("input_GREMLIN"),
        "output_dir"  : os.path.join(os.path.abspath("output_GREMLIN"),pdb),
        "gremlin_dir" : os.path.abspath("GREMLIN"),
        "f_out"       : "output_GREMLIN/{pdb}/{pdb}.gremlin".format(pdb=pdb),
    }

    if os.path.exists(args["f_out"]):
        return f_aln

    cmd = (docker_start +
           "GREMLIN/run_gremlin.sh /opt/mcr/v81 input/{pdb}.aln {f_out}")
    cmd = cmd.format(**args)

    print "Starting", f_aln
    os.system(cmd)
    
    return f_aln


#import itertools
#ITR = itertools.imap(run_aln, F_ALN)

import multiprocessing
P = multiprocessing.Pool()
ITR = P.imap(run_aln, F_ALN)

for f in ITR:
    pass
    #print "Completed", f


