.PHONY: clean

all:
	@echo "Options make [build/clean]"

build:
	python RF.py
	python predict_G2.py
	python build_contacts.py

clean:
	rm -rvf clf G2 predictions
