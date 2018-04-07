# Linking drug target and pathway activation for effective therapy using multi-task learning

## Introduction to Macau

Macau technical paper: http://ieeexplore.ieee.org/document/8168143/

Macau tutorial can be found here: http://macau.readthedocs.io/en/latest/source/examples.html#

To understand more about MCMC Sampling: http://twiecki.github.io/blog/2015/11/10/mcmc-sampling/

![Macau factorization model: (a) The drug response (IC50) is computed by 2 latent matrices. Each of them is being sampled by a Gibbs sampler. In presence of additional information (side information), the latent matrix is predicted by a multiplication of a link matrix and the side information matrix. Arrows in this figure indicate the matrix multiplication. (b) By multiplying the 2 link matrices, we obtain the interaction matrix, which is the interaction between the features of the drugs with the features of the cell lines.
](https://github.com/Katan5555/Macau_project_1/blob/master/image/Keynote_Figure_1.001.png)

## Drug response prediction in 4 different settings: 

Setting 1: macau_GDSC.py

Setting 2: macau_GDSC_pred_new_drug.py

Setting 3: macau_GDSC_transduction.py

Setting 4: macau_GDSC_pred_new_drug_new_cell.py

You may need a certain environment to submit your job (depending on the cluster): cluster_script_GDSC.sh


## Generate the interaction matrix: 

Step 1: use macau_GDSC_pred_new_drug_new_cell.py (for one tissue) or macau_GDSC_interaction_tissue_QC.py (for all tissues) to check the quality of the features. make sure the performance is greater than 0.3

Step 2: Generate the interaction matrix using macau_GDSC_transduction.py

Step 3: use macau_GDSC_interaction_tissue_PERMUTATION.py to do random permutation of the cell lines side features. This gives p-value for each dara point of the interaction matrix

## Plot the results:

Analyse the drug response prediction result: Drug_response_analysis.Rmd

Analysis feature interaction result: GDSC_interaction_analysis_Tissue.Rmd

Scatter plot of the result using: Scatter_plots.Rmd
