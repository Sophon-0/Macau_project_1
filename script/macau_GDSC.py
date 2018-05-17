#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 16:50:34 2016
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
name_of_response = 'IC50_tissue'  ## IC50   IC50_tissue
name_of_cell_feature = 'progeny11'  ## GEX_tissue_OUT , GEX, progeny11 , progeny14 , NES_GDSC_VIPER , NES_GDSC_GSVA
name_of_drug_feature = 'target'  ## target , pathway , target_pathway
## loading data
RES = pd.read_csv(path + 'DATA/' + name_of_response, index_col = False) ; RES = RES.drop(RES.columns[0], axis=1) 
cell_feature = pd.read_csv(path + 'DATA/'  + name_of_cell_feature, index_col = False) ; cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False) ; drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
drug_feature = scipy.sparse.coo_matrix(drug_feature)

repetition = 1  ## should always be 1 or we have to add a new list
nsamples   = 600
nfolds = 10
num_latent = 10

#  result_folder = 'newCell_by_drug_'+name_of_cell_feature+'_'+name_of_drug_feature+'_'+name_of_response+'_'+'rep'+str(repetition)+'_fold'+str(nfolds)+'_sample'+str(nsamples)+'_latent'+str(num_latent)
result_folder = 'newCell_by_drug_'+name_of_cell_feature+'_'+name_of_response+'_'+'rep'+str(repetition)+'_fold'+str(nfolds)+'_sample'+str(nsamples)+'_latent'+str(num_latent)

os.chdir(path + 'DATA_RESULT_STORAGE') 
if os.path.exists(path + 'DATA_RESULT_STORAGE/'+result_folder ):
    os.chdir(result_folder)
else:
    os.makedirs(result_folder) ; os.chdir(result_folder)


## cross validation, split the data
from sklearn import cross_validation

for ite in range (0, repetition):
    kf = cross_validation.KFold(len(RES.iloc[1,:]), n_folds = nfolds, shuffle = True,random_state = random.randint(1, 1000) )

    pcorr_by_drug = [] ## to collect all the rsult from the cross validations
    for train_index, test_index in kf:
       
        cell_feature_train, cell_feature_test = cell_feature.iloc[train_index , : ], cell_feature.iloc[test_index , : ]
        RES_train = RES.iloc[ : , train_index ]
    
        RES_train = scipy.sparse.coo_matrix(RES_train)
        cell_feature_train = scipy.sparse.coo_matrix(cell_feature_train)
    
        ## setting directory 
        os.chdir(path + 'DATA_RESULT_STORAGE/'+result_folder)
        mydir = os.getcwd()+"/"+str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))+"-"+ str(random.randint(1, 1000))
        os.makedirs(mydir) ; os.chdir(mydir)
           
        ## running factorization (Macau)
        result = macau.macau(Y = RES_train,
                         Ytest      = None,
                         side       = [ None , cell_feature_train ],
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
            lmean = np.loadtxt("GDSC-sample%d-U2-latentmean.csv" % N, delimiter=",")
            link  = np.loadtxt("GDSC-sample%d-U2-link.csv" % N,       delimiter=",")
            V     = np.loadtxt("GDSC-sample%d-U1-latents.csv" % N,    delimiter=",")
    
            ## predicted latent vector for new cell lines from sample N
            uhat = cell_feature_test.dot(link.transpose()) + lmean
            ## use predicted latent vector to predict reponse for new cell lines
            Yhat = uhat.dot(V) + meanvalue
            Ystore.append(Yhat) 
    
        dfs = Ystore
        Yhat_avg = pd.concat([each.stack() for each in dfs],axis=1)\
                 .apply(lambda x:x.mean(),axis=1)\
                 .unstack()
    
        ### compute pearson correlation by drug ###
        Ypred_by_drug = Yhat_avg.transpose()
        from scipy import stats
        Yobs_by_drug = RES.iloc[ : , Ypred_by_drug.columns ]
        
        pp_by_drug = []
        for d in range(0 , len(Yobs_by_drug) ):
           ii = np.flatnonzero(np.isfinite( Yobs_by_drug.iloc[ d , : ] ))
           pcorr = stats.pearsonr( Yobs_by_drug.iloc[ d , ii ]  , Ypred_by_drug.iloc[ d , ii ]  ) [0]
           pp_by_drug.append(pcorr)
       
        pp_by_drug = pd.DataFrame(pp_by_drug) 
        pcorr_by_drug.append(pp_by_drug.iloc[:,0])  
        
        filelist = glob.glob("*.csv")
        for f in filelist:
           os.remove(f)
        os.rmdir(mydir)
        
                
##################################################################### save result ##################################################################
os.chdir( path + 'DATA_RESULT_STORAGE/'+result_folder )  
pcorr_by_drug = pd.DataFrame(pcorr_by_drug)  ;  pcorr_by_drug = np.mean(pcorr_by_drug)
pcorr_by_drug.to_csv('pcorr_'+result_folder+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+str(random.randint(1, 10000))+'.csv')  

     
#######################################################  AVERAGING RESULT FROM MULTIPLE RUNS #######################################################
prefixed = [filename for filename in os.listdir('.') if filename.startswith("pcorr_newCell_by_drug")]
if len(prefixed) >= 20:      
    list = []
    for file in prefixed:
        df = pd.read_csv(file, header = None)  ; df = df.iloc[:,1]
        list.append(df)
    list = pd.DataFrame(list)
    mean = np.mean(list)
    mean.to_csv('pcorr_'+result_folder+'.csv')







