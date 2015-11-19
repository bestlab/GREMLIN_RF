Use GREMLIN to compute the coevol information matrices.

First use `format_aln.py` to convert the `aln` directory into a minimial alignment, this creates single line alignments in `input_GREMLIN`.

GREMLIN, a MATLAB program, provides a stand-alone precompiled version, but this still requires the MATLAB MCL libraries to be installed locally. We provide a [docker](docker_matlab) script that builds the libraries into a container.

With the container built, `run_GREMLIN.py` will compute the coevolutationary alignments in parallel. 

When this is complete, run `collate_GREMLIN.py` to save and finish the APC-corrected results (an NxN matrix).
