import glob, os, itertools, multiprocessing

#F_SILENT = glob.glob("clusters/clustered_outfiles/*.out")
F_SILENT = glob.glob("clusters/combined/*.out")

def extract(f):
    
    name = os.path.basename(f).split('.out')[0]
    dir_name = os.path.join('pdb',name) + '/'
    os.system('mkdir -p '+dir_name)

    base_pdb = os.path.join(dir_name, 'c.0.0.pdb.pdb')
    if os.path.exists(base_pdb):
        return None

    cmd = ("./extract_pdbs.linuxgccrelease "
           "-in::file::silent {input} "
           "-out::prefix {outfile}")

    
    os.system( cmd.format(input=f,outfile=dir_name) )
    return f


os.system('mkdir -p pdb')

ITR = itertools.imap(extract,F_SILENT)
P = multiprocessing.Pool()
ITR = P.imap(extract,F_SILENT)

for f in ITR:
    print "Completed", f
