# Linking drug target and pathway activation for effective therapy using multi-task learning

## Drug response prediction in 4 different settings: 

Setting 1: macau_GDSC.py

Setting 2: macau_GDSC_pred_new_drug.py

Setting 3: macau_GDSC_transduction.py

Setting 4: macau_GDSC_pred_new_drug_new_cell.py

## Generate the interaction matrix: 

Step 1: use macau_GDSC_pred_new_drug_new_cell.py (for one tissue) or macau_GDSC_interaction_tissue_QC.py (for all tissues) to check the quality of the features. make sure the performance is greater than 0.3

Step 2: Generate the interaction matrix using macau_GDSC_transduction.py

Step 3: use macau_GDSC_interaction_tissue_PERMUTATION.py to do random permutation of the cell lines side features. This gives p-value for each dara point of the interaction matrix


