all:
	echo "pass"

full_clean:
	rm -rvf systems/ energy_table/

build:
	make full_clean
	python build_systems.py
	python run_systems.py
	python compute_rg.py
	python compute_E.py

	python plot_gr.py &
	python plot_E.py &

	cd systems/1a3a_exact_1.00_seed_0/ && echo 0 | trjconv -f traj.xtc -o movie.pdb && pymol -d "show spheres; set sphere_scale, 4.0" movie.pdb


extract:
	python extract_pdb.py
	python extract_coordinates.py
