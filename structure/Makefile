all:
	echo "pass"

full_clean:
	rm -rvf systems/ energy_table/

build:
#	make full_clean
#	python build_scaled_systems.py
	python run_systems.py
	make extract

extract:
	python extract_pdb.py
	python extract_coordinates.py
#	python extract_contact_maps.py

stats:
	python compute_Q.py
	python compute_rg.py
	python compute_E.py
	python compute_RMSD.py
	python compute_dRMSD.py

full:
	make build
	make extract
	make stats

dRMSD:
	python compute_dRMSD.py  && python plot_ALL_dRMSD.py

Q:
	python compute_Q.py  && python plot_ALL_Q.py
