all:
	echo "pass"

full_clean:
	rm -rvf systems/ energy_table/

build:
	python build_systems.py
	python run_systems.py
	python compute_rg.py
	python plot_gr.py
