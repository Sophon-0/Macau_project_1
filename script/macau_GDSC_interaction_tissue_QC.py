#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 18:17:49 2017
@author: miyang
"""

import scipy.io
import glob, os, datetime
import pandas as pd
import numpy as np
import macau
import random

#  origin_path = '/Users/miyang/Documents/RWTH_Aachen/'
origin_path = '/home/my871390/MI_YANG/RWTH_Aachen/'

path = origin_path+'macau_work_dir/macau_test_sanger/'
path_tissue = origin_path+'SANGER_DATA/TISSUE/'


os.chdir(path_tissue)
prefixed = [filename for filename in os.listdir(path_tissue) if not filename.startswith('.') ] 
prefixed = sorted(prefixed)
            

for i in range( 8 , len(prefixed) ):      ## len(prefixed)      
    ## SELECT THE TISSUE
    os.chdir(path_tissue + prefixed[i] + '/') 
    
    ## SELECT THE DATA
    name_of_response = 'RES'         ## RES   RES_target_Leiden  
    name_of_cell_feature = 'GEX_SLC_ABC' ## progeny11 progeny14 progeny_old11_new3 GEX GEX_by_gene GEX_SLC_ABC SNP_CNV TF_GSVA TF_VIPER 
    name_of_drug_feature = 'target'  ## target target_Leiden ECFP4
           
    RES = pd.read_csv(name_of_response, index_col = False); RES = RES.drop(RES.columns[0], axis=1) 
    cell_feature = pd.read_csv(name_of_cell_feature, index_col = False); cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
    drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False); drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
    # drug_feature = pd.read_csv(name_of_drug_feature, index_col = False); drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
    
    nsamples   = 600
    nfolds = 10
    num_latent = 30
    prefix = "GDSC"
    
    result_folder = origin_path+'MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/'+ name_of_drug_feature+'_'+name_of_cell_feature+'_QC'+"/"+prefixed[i]
    if os.path.exists(result_folder):
       os.chdir(result_folder)
    else:
       os.makedirs(result_folder) ; os.chdir(result_folder)
    
           
    ## cross validation, split the data
    from sklearn import cross_validation
    
    ############  2 simultaneous loops  ############
    kf_cell = cross_validation.KFold(len(RES.iloc[1,:]), n_folds = nfolds, shuffle = True,random_state = random.randint(1, 1000) )
    kf_drug = cross_validation.KFold(len(RES.iloc[:,1]), n_folds = nfolds, shuffle = True,random_state = random.randint(1, 1000) )
    
    pp = []
    for (drug_train_index, drug_test_index),(cell_train_index, cell_test_index) in zip(kf_drug, kf_cell):
         
        drug_features_train, drug_features_test = drug_feature.iloc[ drug_train_index  , : ], drug_feature.iloc[ drug_test_index , : ]
        drug_features_train = scipy.sparse.coo_matrix(drug_features_train)

        cell_features_train, cell_features_test = cell_feature.iloc[ cell_train_index  , : ], cell_feature.iloc[ cell_test_index , : ]
        cell_features_train = scipy.sparse.coo_matrix(cell_features_train)
        
        RES_train = RES.iloc[ drug_train_index , cell_train_index ]
        RES_train = scipy.sparse.coo_matrix(RES_train)
        
        ## setting directory 
        os.chdir(result_folder)
        mydir = os.getcwd()+"/"+str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))+"-"+ str( random.randint(1, 1000) )
        os.makedirs(mydir) ;  os.chdir(mydir)
        
        result = macau.macau(Y = RES_train,
                         Ytest      = None,
                         side       = [drug_features_train, cell_features_train],
                         num_latent = num_latent,
                         precision  = "adaptive",
                         burnin     = 400,
                         nsamples   = nsamples,
                         univariate = True,
                         save_prefix= prefix)
    
        meanvalue = np.loadtxt("GDSC-meanvalue.csv").tolist()
            
        Ystore = list()
        ## loading N samples 
        for N in range(1, nsamples + 1): 
            lmean_U1 = np.loadtxt(prefix+"-sample%d-U1-latentmean.csv" % N, delimiter=",")
            link_U1  = np.loadtxt(prefix+"-sample%d-U1-link.csv" % N,       delimiter=",")
            lmean_U2 = np.loadtxt(prefix+"-sample%d-U2-latentmean.csv" % N, delimiter=",")
            link_U2  = np.loadtxt(prefix+"-sample%d-U2-link.csv" % N,       delimiter=",")
    
            ## predicted latent vector for new cell lines and for new drugs
            uhat_U1 = drug_features_test.dot(link_U1.transpose()) + lmean_U1
            uhat_U2 = cell_features_test.dot(link_U2.transpose()) + lmean_U2

            Yhat = uhat_U1.dot(uhat_U2.transpose()) + meanvalue
            Ystore.append(Yhat) 
    
        dfs = Ystore
        Yhat_avg = pd.concat([each.stack() for each in dfs],axis=1)\
                 .apply(lambda x:x.mean(),axis=1)\
                 .unstack()
    
        ### compute pearson correlation  ###
        Ypred = Yhat_avg
        from scipy import stats
        Yobs = RES.iloc[ Ypred.index , Ypred.columns ]
        
        Ypred_vec = np.ravel(Ypred).T
        Yobs_vec = np.ravel(Yobs).T
      
        ii = np.flatnonzero(np.isfinite(Yobs_vec))
        pcorr = stats.pearsonr( Ypred_vec[ii]  , Yobs_vec[ii]  ) [0]        
        pp.append(pcorr)
                     
        filelist = glob.glob("*.csv")
        for f in filelist:
           os.remove(f)
        os.rmdir(mydir)
    
    os.chdir(result_folder)          
    pcorr_test_set = pd.DataFrame(pp)
    pcorr_test_set.to_csv('pcorr_'+prefixed[i]+'_'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+str(random.randint(1, 1000))+'.csv') 
      
    #######################################################  AVERAGING RESULT FROM MULTIPLE RUNS #######################################################
    prefixed_intermediate = [filename for filename in os.listdir('.') if filename.startswith("pcorr_"+prefixed[i] )]
    
    if len(prefixed_intermediate) >= 20:      
        final_result = []
        for file in prefixed_intermediate:
            df = pd.read_csv(file) ; df = df.iloc[:,1]
            final_result.append(df)
        
        mean = np.mean(final_result) ; os.makedirs('mean='+str(round(mean, 3))) 
        sd   = np.std(final_result)  ; os.makedirs('sd='  +str(round(sd  , 3))) 
       
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    