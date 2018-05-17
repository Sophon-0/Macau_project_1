#!/usr/bin/env python2
#import sys
#tissue = sys.argv[1]
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 10:39:44 2018
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

      
## SELECT THE TISSUE
i = 2
os.chdir(path_tissue + prefixed[i] + '/') 

## SELECT THE DATA
name_of_response = 'RES'   ##  RES   RES_GSVA_Reactome  RES_target_Leiden  RES_CRISPR
name_of_cell_feature = 'progeny11'   ## progeny11  GEX  SNP_CNV  CRISPR
name_of_drug_feature = 'target'  ## target  target_Leiden
interaction_folder = origin_path+'MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/'+ name_of_drug_feature+'_'+name_of_cell_feature 
count_folder = origin_path+'MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/'+ name_of_drug_feature+'_'+name_of_cell_feature+'_PERMUTATION/'+prefixed[i]

nsamples   = 600
num_latent = 30
prefix = "GDSC"
number_permutation = 1000
reference_interaction = pd.read_csv(interaction_folder+'/interaction_'+prefixed[i]+'_'+name_of_drug_feature+'_'+name_of_cell_feature +'_sample'+str(nsamples)+'_latent'+str(num_latent)+'.csv')
reference_interaction = reference_interaction.drop(reference_interaction.columns[0], axis=1) 


count_store = pd.DataFrame()  ##  COUNT MATRIX   

## PERMUTATION
for permut in range(0,number_permutation):

    os.chdir(path_tissue + prefixed[i] + '/') 
    RES = pd.read_csv(name_of_response, index_col = False) ; RES = RES.drop(RES.columns[0], axis=1) 
    cell_feature = pd.read_csv(name_of_cell_feature, index_col = False) ; cell_feature = cell_feature.drop(cell_feature.columns[0], axis=1)
    drug_feature = pd.read_csv(path + 'DATA/'  + name_of_drug_feature, index_col = False) ; drug_feature = drug_feature.drop(drug_feature.columns[0], axis=1)
   
    side_1_name = drug_feature.columns
    side_2_name = cell_feature.columns
    
    ## PERMUTATION of cell line feature
    for sample in range(0,len(cell_feature.index)): 
        cell_feature.iloc[sample, : ] = np.random.permutation( cell_feature.iloc[sample, : ] )
    
    RES = scipy.sparse.coo_matrix(RES)
    drug_feature = scipy.sparse.coo_matrix(drug_feature)
    cell_feature = scipy.sparse.coo_matrix(cell_feature)
    
    ## running factorization (Macau)
    if os.path.exists(count_folder):
        os.chdir(count_folder)
    else:
        os.makedirs(count_folder) ; os.chdir(count_folder)
        
    mydir = os.getcwd()+"/"+str(datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))+"-"+ str( random.randint(1, 1000) )
    os.makedirs(mydir) ;  os.chdir(mydir)
    
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

    ##### UPDATE COUNT MATRIX
    reference_interaction.index = interaction_mean.index
    mat = interaction_mean.subtract(reference_interaction) 
    mat[mat > 0] = 1  ## when greater than reference value
    mat[mat < 0] = 0

    count_store = count_store.add(mat, fill_value=0)
    
                
os.chdir(count_folder)
count_store.to_csv('interaction_COUNT_'+prefixed[i]+'_'+name_of_drug_feature+'_'+name_of_cell_feature +'_sample'+str(nsamples)+'_latent'+str(num_latent)+'.csv') 

    
    
