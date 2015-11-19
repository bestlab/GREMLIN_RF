Compute the random forest model over the APC GREMLIN matrices.

Models are trained with 

    python RF.py --kernel_window 2 --n_estimators 200 --ratio_TP_to_TN 20 --k_folds 5

After training run `predict_G2.py` to generate the new updated matrices. This populates the directory `G2`.

From here run `build_contacts.py` to build the ordered set of contact maps for given cutoff values.