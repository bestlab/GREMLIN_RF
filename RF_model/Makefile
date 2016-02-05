.PHONY: clean

RF_ARGS = --kernel_window 2 --n_estimators 200 --ratio_TP_to_TN 20 --k_folds 5

all:
	@echo "Options make [build/clean]"

build:
	python RF.py $(RF_ARGS)
	python predict_G2.py
	python build_contacts.py

clean:
	rm -rvf clf G2 predictions
