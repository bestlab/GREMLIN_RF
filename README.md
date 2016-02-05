# GREMLIN Random Forests
Enchancing coevolutionary results with random forests.

------------------------------------------------------------

This protocol assumes you have HHBlits and GREMLIN installed. Run `make build` in each of the directories in order

### `preprocessing/`

Follow the instructions in [/preprocessing/pdb](/preprocessing/org_pdb) to prepare a set of fasta and pdb files.

### `alignment/`

### `coevol/`

GREMLIN is packaged with this program. The entire GREMLIN directory is released under a Attribution-NonCommercial-ShareAlike 3.0 Unported licence.

### `RF_model/`

Train the extremely random forests with scikit-learn.

### `structure/`

Test structure files with a simple C-alpha folding model.

### `analysis/`

Run post-processing analysis.

----------------------------------------------------------------------

Presention on research given 2/4/2016
[ENHANCING THE COEVOLUTIONARY SIGNAL](http://thoppe.github.io/Presentation_Enhancing_coevol_signal/Brooks_presentation.html#/)


