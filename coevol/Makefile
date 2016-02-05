.PHONY: clean

all:
	@echo "Options make [build/clean]"

build:
	python format_aln.py
	python run_GREMLIN.py
	python collate_GREMLIN.py

clean:
	rm -rvf APC GREMLIN
