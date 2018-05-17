#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 14:10:39 2017
@author: miyang
"""


import scipy.io
import glob, os, datetime
import pandas as pd
import numpy as np
import macau

origin_path = '/Users/miyang/Documents/RWTH_Aachen/'
#  origin_path = '/home/my871390/MI_YANG/RWTH_Aachen/'

path = origin_path+'macau_work_dir/macau_test_sanger/'
path_tissue = origin_path+'SANGER_DATA/TISSUE/'


os.chdir(path_tissue)
prefixed = [filename for filename in os.listdir(path_tissue) if not filename.startswith('.') ] 
prefixed = sorted(prefixed)

for i in range( 0 , len(prefixed) ):      ## len(prefixed)      
    ## SELECT THE TISSUE
    os.chdir(path_tissue + prefixed[i] + '/') 
    
    ## SELECT THE DATA
    name_of_response = 'RES'   ##  RES   RES_GSVA_Reactome  RES_target_Leiden  RES_CRISPR
    name_of_cell_feature = 'progeny11'   ## progeny11  progeny14  GEX  SNP_CNV  TF_score 
    name_of_drug_feature = 'target'  ## target  target_Leiden
    result_folder = origin_path+'MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/'+ name_of_drug_feature+'_'+name_of_cell_feature
    
    RES = pd.read_csv(name_of_response, index_col = False) ; RES = RES.drop(RES.columns[0], axis=1) 
    cell_feature = pd.read_csv(name_of_cell_feature, index_col = False) ; cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
    drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False) ; drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
   
    side_1_name = drug_feature.columns
    side_2_name = cell_feature.columns
    
    nsamples   = 600
    num_latent = 30
    prefix = "GDSC"
    
    RES = scipy.sparse.coo_matrix(RES)
    drug_feature = scipy.sparse.coo_matrix(drug_feature)
    cell_feature = scipy.sparse.coo_matrix(cell_feature)
    
    ## running factorization (Macau)
    os.chdir(path + "MACAU_SAVE_TEST")
    mydir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    os.makedirs(mydir) ;             
    os.chdir(mydir)
    
    ## running factorization (Macau)
    result = macau.macau(Y = RES,
                     Ytest      = None,
                     side       = [ drug_feature, cell_feature ],
                     num_latent = num_latent,
                     precision  = "adaptive",
                     burnin     = 400,
                     nsamples   = nsamples,
                     univariate = True,
                     save_prefix = prefix)
    
    meanvalue = np.loadtxt(prefix + "-meanvalue.csv").tolist()
    
    
    interaction_store = list()
    for N in range(1, nsamples + 1 ): 
        link_1  = np.loadtxt(prefix + "-sample%d-U1-link.csv" % N, delimiter=",")
        link_1 = pd.DataFrame(link_1)
        
        link_2  = np.loadtxt(prefix + "-sample%d-U2-link.csv" % N, delimiter=",")
        link_2 = pd.DataFrame(link_2)
        link_2 = link_2.transpose()
        
        interaction = link_2.dot(link_1)   ;  interaction = interaction.transpose()
        interaction_store.append(interaction)
        
    interaction_mean = pd.concat([each.stack() for each in interaction_store],axis=1).apply(lambda x:x.mean(),axis=1).unstack()     
    
    
    filelist = glob.glob("*.csv")
    for f in filelist:
       os.remove(f)
    os.rmdir(mydir) 
    
    
    ##### give the label back and multiply by -1 
    interaction_mean = interaction_mean * -1
    interaction_mean.index = side_1_name  
    interaction_mean.columns = side_2_name           
                
    if os.path.exists(result_folder):
        os.chdir(result_folder)
    else:
        os.makedirs(result_folder) ; os.chdir(result_folder)

    interaction_mean.to_csv('interaction_'+prefixed[i]+'_'+name_of_drug_feature+'_'+name_of_cell_feature +'_sample'+str(nsamples)+'_latent'+str(num_latent)+'.csv') 


    
    
