#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 19:45:20 2016
@author: miyang
"""

import macau
import scipy.io
import glob, os, datetime
import pandas as pd
import numpy as np
import os.path
import random
#  path = "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/"
path = "/home/my871390/MI_YANG/RWTH_Aachen/macau_work_dir/macau_test_sanger/"


## loading data
name_of_response = 'IC50'    ##   IC50_ECFP4  IC50_target_Leiden  IC50
name_of_cell_feature = 'GEX'  ## GEX progeny11
name_of_drug_feature = 'target'   ##   pathway , target , target_pathway , Target_Expansion  ECFP4  target_Leiden
## loading data
RES = pd.read_csv(path + 'DATA/' + name_of_response, index_col = False) ; RES = RES.drop(RES.columns[0], axis=1) 
cell_feature = pd.read_csv(path + 'DATA/'  + name_of_cell_feature, index_col = False) ; cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
cell_feature = scipy.sparse.coo_matrix(cell_feature)
drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False) ; drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)

repetition = 1  ## should always be 1 or we have to add a new list
nsamples   = 600
nfolds = 10
num_latent = 10

# result_folder = 'newDrug_by_cell_'+name_of_drug_feature+'_'+name_of_cell_feature+'_'+name_of_response+'_'+'rep'+str(repetition)+'_fold'+str(nfolds)+'_sample'+str(nsamples)+'_latent'+str(num_latent)
result_folder = 'newDrug_by_cell_'+name_of_drug_feature+'_'+name_of_response+'_'+'rep'+str(repetition)+'_fold'+str(nfolds)+'_sample'+str(nsamples)+'_latent'+str(num_latent)

os.chdir(path + 'DATA_RESULT_STORAGE') 
if os.path.exists(path + 'DATA_RESULT_STORAGE/'+result_folder ):
    os.chdir(result_folder)
else:
    os.makedirs(result_folder) ; os.chdir(result_folder)


## cross validation, split the data
from sklearn import cross_validation

for ite in range (0, repetition):
    
    kf = cross_validation.KFold(len(RES.iloc[:,1]), n_folds = nfolds, shuffle = True,random_state = random.randint(1, 10000) )
    
    pcorr_by_cell = []  ## to collect all the rsult from the cross validations
    for train_index, test_index in kf:
      
        drug_feature_train, drug_features_test = drug_feature.iloc[ train_index  , : ], drug_feature.iloc[ test_index , : ]
        drug_feature_train = scipy.sparse.coo_matrix(drug_feature_train)

        RES_train = RES.iloc[ train_index , : ]
        RES_train = scipy.sparse.coo_matrix(RES_train)
        
        ## setting directory 
        os.chdir(path + 'DATA_RESULT_STORAGE/'+result_folder)
        mydir = os.getcwd()+"/"+str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))+"-"+ str(random.randint(1, 10000))
        os.makedirs(mydir) ;  os.chdir(mydir)
        
        result = macau.macau(Y = RES_train,
                         Ytest      = None,
                         side       = [drug_feature_train, None ],
                         num_latent = num_latent,
                         precision  = "adaptive",
                         burnin     = 400,
                         nsamples   = nsamples,
                         univariate = True,
                         save_prefix = "GDSC")
    
        meanvalue = np.loadtxt("GDSC-meanvalue.csv").tolist()
    
        Ystore = list()
        ## loading N samples 
        for N in range(1, nsamples + 1): 
            lmean = np.loadtxt("GDSC-sample%d-U1-latentmean.csv" % N, delimiter=",")
            link  = np.loadtxt("GDSC-sample%d-U1-link.csv" % N,       delimiter=",")
            V     = np.loadtxt("GDSC-sample%d-U2-latents.csv" % N,    delimiter=",")
    
            ## predicted latent vector for new cell lines from sample N
            uhat = drug_features_test.dot(link.transpose()) + lmean
            ## use predicted latent vector to predict reponse for new cell lines across drugs
            Yhat = uhat.dot(V) + meanvalue
            Ystore.append(Yhat) 
    
        dfs = Ystore
        Yhat_avg = pd.concat([each.stack() for each in dfs],axis=1)\
                 .apply(lambda x:x.mean(),axis=1)\
                 .unstack()
    
        
        ### compute pearson correlation by cell line ###
        Ypred_by_cell = Yhat_avg
        from scipy import stats
        Yobs_by_cell = RES.iloc[ Ypred_by_cell.index , : ]
        pp_by_cell = []
        for d in range(0 , len(Yobs_by_cell.iloc[0,:]) ):
           ii = np.flatnonzero(np.isfinite(Yobs_by_cell.iloc[ : , d ]))
           pcorr = stats.pearsonr( Yobs_by_cell.iloc[ ii , d ] , Ypred_by_cell.iloc[ ii , d ] ) [0]
           pp_by_cell.append(pcorr)
        
        pp_by_cell = pd.DataFrame(pp_by_cell) 
        pcorr_by_cell.append(pp_by_cell.iloc[:,0])  
            
        filelist = glob.glob("*.csv")
        for f in filelist:
           os.remove(f)
        os.rmdir(mydir) 

                     
os.chdir( path + 'DATA_RESULT_STORAGE/'+ result_folder ) 
pcorr_by_cell = pd.DataFrame(pcorr_by_cell) ;  pcorr_by_cell = np.mean(pcorr_by_cell)
pcorr_by_cell.to_csv('pcorr_'+result_folder+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+str(random.randint(1, 10000))+'.csv') 


    
#######################################################  AVERAGING RESULT FROM MULTIPLE RUNS #######################################################
prefixed = [filename for filename in os.listdir('.') if filename.startswith("pcorr_newDrug_by_cell")]

if len(prefixed) >= 20:      
    list = []
    for file in prefixed:
        df = pd.read_csv(file, header = None)  ; df = df.iloc[:,1]
        list.append(df)
    list = pd.DataFrame(list)
    mean = np.mean(list)
    mean.to_csv('pcorr_'+result_folder+'.csv')
####################################################################################################################################################
    


