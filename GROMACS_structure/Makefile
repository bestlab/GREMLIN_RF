all:
	echo "pass"

full_clean:
	rm -rvf systems/ energy_table/

build:
	python build_systems.py
	python run_systems.py
	python compute_rg.py
	python compute_E.py

	python plot_gr.py

cd systems/1a3a_exact_1.00_seed_0/
