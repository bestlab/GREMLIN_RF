
build_all:
	make clean
	python collate_rossetta.py
	make compute

compute:
	python score_best.py
	python collate_tmscores.py
	python plot_tmscores.py


clean:
	rm -rvf pdb tmscores
	rm -vf  master_tmscores.db master_rosetta_scores.db
