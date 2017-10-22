#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 16:50:34 2016
@author: miyang
"""

import macau
import scipy.io
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import mean_squared_error
path = "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/"
# path = "/home/my871390/MI_YANG/RWTH_Aachen/macau_work_dir/macau_test_sanger/"


## loading data
name_of_response = 'IC50'
name_of_cell_feature = 'GEX'
name_of_drug_feature = 'target'
RES = pd.read_csv(path + 'DATA/' + name_of_response, index_col = False) ; RES = RES.drop(RES.columns[0], axis=1) 
cell_feature = pd.read_csv(path + 'DATA/'  + name_of_cell_feature, index_col = False) ; cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False) ; drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
drug_feature = scipy.sparse.coo_matrix(drug_feature)

nsamples   = 600
num_latent = 10

RES = scipy.sparse.coo_matrix(RES)
cell_feature = scipy.sparse.coo_matrix(cell_feature)

repetition = 10

for ite in range (0, repetition):    
    ## running factorization (Macau)
    result = macau.macau(Y = RES,
                     Ytest      = 0.7,
                     side       = [None,None],
                     num_latent = num_latent,
                     precision  = "adaptive",
                     burnin     = 400,
                     nsamples   = nsamples,
                     univariate = True)
    
    result
                           
    pcorr = stats.pearsonr(result.prediction.y, result.prediction.y_pred)[0]
    RMSE = mean_squared_error(result.prediction.y, result.prediction.y_pred)

    if ite == 0:
        pcorr_iteration = pcorr
        RMSE_iteration = RMSE
    else:
        pcorr_iteration = np.vstack((pcorr_iteration,pcorr)) 
        RMSE_iteration = np.vstack((RMSE_iteration,RMSE)) 

pcorr = round( np.mean(pcorr_iteration) , 3 ) 
sd_pcorr = round( np.std(pcorr_iteration) , 5 ) 

RMSE = round( np.mean(RMSE_iteration) , 3 ) 
sd_RMSE = round( np.std(RMSE_iteration) , 5 ) 

