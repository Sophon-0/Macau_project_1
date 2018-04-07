path <- "~/Documents/RWTH_Aachen"
source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))
load(paste(path,"/SANGER_DATA/v17a_IC50s.Rdata", sep=""))
DRUG_ANALYSIS_SET <- read.csv("~/Documents/RWTH_Aachen/SANGER_DATA/DRUG_ANALYSIS_SET",row.names = 1, check.names = F)
library(openxlsx)
DRUG_ANALYSIS_SET_update <- read.xlsx("/Users/miyang/Documents/RWTH_Aachen/SANGER_DATA/Screened_Compounds.xlsx",1)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

IC50 <- t(v17a_IC50s); IC50 <- convert_drugID(IC50)
v21.meta.per_compound <- read.delim("~/Documents/RWTH_Aachen/CTD2/Basal_Gene_Expression_and_Copy_Number_Correlation_analysis/CTRPv2.1_2016_pub_NatChemBiol_12_109/v21.meta.per_compound.txt")

print_target_GDSC <- function(protein_target, target_matrix, drug_names , prediction=F) {
  target_number  <- which(colnames(target_matrix)==protein_target)
  drug_concerned <- rownames(target_matrix)[which(target_matrix[ ,target_number] ==1)]
  drug_concerned <- drug_names[as.numeric(drug_concerned)]
  DRUG_ANALYSIS_SET$DRUG_NAME <- as.character(DRUG_ANALYSIS_SET$DRUG_NAME)
  DRUG_ANALYSIS_SET$PUTATIVE_TARGET <- DRUG_ANALYSIS_SET_update$Target
  subset <- DRUG_ANALYSIS_SET[ which(DRUG_ANALYSIS_SET$DRUG_NAME %in% drug_concerned) , c(3,6)  ]
  
  if(prediction==T) {
    drug_prediction <- read.csv("~/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA_RESULT_STORAGE/pcorr_newCell_by_drug_GEX_target_IC50_rep20_fold10_sample600_latent30.csv", header=FALSE)
    IC50 <- read.csv("~/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA/IC50")
    drug_prediction$V1 <- IC50$X
    prediction <- drug_prediction$V2[which(drug_prediction$V1 %in% subset[ ,1])]
    subset <- cbind(subset, prediction)
  }
  return(subset)
}


print_target_CTRP <- function(protein_target, target_matrix, drug_names ) {
  target_number  <- which(colnames(target_matrix)==protein_target)
  drug_concerned <- which(target_matrix[ ,target_number] ==1)
  drug_concerned <- drug_names[drug_concerned]
  v21.meta.per_compound$cpd_name <- as.character(v21.meta.per_compound$cpd_name)
  subset <- v21.meta.per_compound[ which(v21.meta.per_compound$cpd_name %in% drug_concerned) , c(2,6)  ]
  return(subset)
}

target_combo_gene <- function(mat, association,  corr_diff)  {
  t <- rowSums(mat) ; t <- t[order(-t)] ; print("Most efficient targets: ") ;print(names(t)[1:5] ) 
  # boxplot( rowSums(mat), main="Total target effect captured by RNA",medlwd=0.5)  
  mat2 <- mat[ rownames(mat)[which(rowSums(mat) > association)] , ] ## select strongest targets
  mat2 <- subset_col_up(mat2, up_limit = 0.9, obs=1)  ## remove useless genes
  
  cor_mat2 <- cor(t(mat2))
  coor <- which(cor_mat2 < corr_diff, arr.ind = T) 
  
  top_hits <- cbind(colnames(cor_mat2)[ coor[ ,1] ] , colnames(cor_mat2)[ coor[ ,2] ], cor_mat2[coor] ) 
  top_hits <- data.frame(top_hits) ; colnames(top_hits) <- c("target1", "target2" , "corr_diff")
  top_hits <- top_hits[!duplicated(top_hits[,c("corr_diff")]), ] ; 
  top_hits <- top_hits[ order(as.numeric(as.character(top_hits$corr_diff))) , ] 
  top_hits$corr_diff <- round(as.numeric(as.character(top_hits$corr_diff)), digits = 3)
  rownames(top_hits)<-1:length(rownames(top_hits))
  
  L <- list(cor_mat2,top_hits)
  return(L)
}

target_combo_pathway <- function(mat, top_association,  correlation)  {
  t <- rowSums(mat) ; t <- t[order(-t)] ; ## print("Most efficient targets: ") ; print(names(t)[1:5] ) 
  # boxplot( rowSums(mat), main="Total target effect captured by RNA",medlwd=0.5)  
  mat2 <- cbind(apply(mat,1,max,na.rm=TRUE) , mat)
  mat2 <- mat2[ order(-mat2$`apply(mat, 1, max, na.rm = TRUE)`) , ] ## add column max effect on single pathway
  mat2 <- mat2[ ,  -mat2$`apply(mat, 1, max, na.rm = TRUE)` ]
  mat2 <- mat2[ 1:top_association , ]  ## take the top single pathways effect
  
  cor_mat2 <- cor(t(mat2))
  coor <- which(cor_mat2 < correlation, arr.ind = T) 
  
  top_hits <- cbind(colnames(cor_mat2)[ coor[ ,1] ] , colnames(cor_mat2)[ coor[ ,2] ], cor_mat2[coor] ) 
  top_hits <- data.frame(top_hits) ; colnames(top_hits) <- c("target1", "target2" , "correlation")
  top_hits <- top_hits[!duplicated(top_hits[,c("correlation")]), ] ; 
  top_hits <- top_hits[ order(as.numeric(as.character(top_hits$correlation))) , ] 
  top_hits$correlation <- round(as.numeric(as.character(top_hits$correlation)), digits = 3)
  rownames(top_hits)<-1:length(rownames(top_hits))
  
  L <- list(cor_mat2,top_hits)
  return(L)
}

## drug_feature_name="target" ; cell_feature_name="GEX_SLC_ABC" ; selection="conservative" ; cut_off="single" ; N=9
save_interaction_plot <- function(drug_feature_name,cell_feature_name,selection,cut_off,N = 1:length(f) , limit_row=0.95,limit_col=0.95,significance=0.05) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,"/",sep = "")
  path_permutation <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,"_PERMUTATION/" )
  setwd(path) ; f <- list.files(path)
  for(i in N) {
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    
    ##################### decide about target selection: conservative or all
    if(selection == "conservative") { x <- x[ - target_to_remove ,  ] ; mat <- x }
    if(selection == "conservative_include_BIRC5") { mat <- x[ - target_to_remove ,  ] ; mat <- rbind(mat,x["BIRC5", ]) }
    if(selection == "all_target") {  mat <- x }
    ##################### decide about target cut off: single top hit of sum across all
    if(cut_off=="absolute") {
      v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
      if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
      mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:25] ]
      } 
    }
    if(cut_off=="single") {
      mat <- subset_row_abs(mat, abs_limit = limit_row, obs=1) ;  v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
      if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
      mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:25] ]
      } 
    }
    mat <- mat[complete.cases(mat),]
    tissue <- regmatches(f[i], regexpr('interaction.+?target', f[i])) ; tissue <- gsub('.{7}$', '', tissue) ; tissue <- substring(tissue, 13)
    
    interaction_pvalue_CORRECTED <- read.csv(paste0(path_permutation,tissue,"/interaction_pvalue_CORRECTED"), row.names = 1)
    interaction_pvalue_CORRECTED <- interaction_pvalue_CORRECTED[rownames(mat),colnames(mat)]
    
#     ######## For target Leiden
#     table_name_conversion <- read.csv("~/Documents/RWTH_Aachen/DIVERSE_AND_PROJECT/LEIDEN/DATA/table_name_conversion", row.names=1)
#     common <- intersect(table_name_conversion$pref_name, rownames(mat))
#     table_name_conversion <- table_name_conversion[table_name_conversion$pref_name  %in% common , ] 
#     mat <- mat[common, ] ; interaction_pvalue_CORRECTED <- interaction_pvalue_CORRECTED[common, ]
#     rownames(mat) <- table_name_conversion$HGNC ; rownames(interaction_pvalue_CORRECTED) <- table_name_conversion$HGNC
    
    pdf(paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/PLOTS/GDSC_",drug_feature_name,"_",cell_feature_name,"/",selection,'/',cut_off,'/',tissue,".pdf", sep = ""), width = 14.6 , height = 15.5, onefile = F ) # width = 14 , height = 15
    
    library(pheatmap) ; library(grid)
    plot_pheatmap <- function(mat, row_names, col_names , title ,cluster_rows=T,cluster_cols=T,fontsize=33,fontsize_row=33, fontsize_col=33, scale="none") {
      setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.95, name="vp", just=c("right","top"))), action="prepend")
      pheatmap(mat, main=title, fontsize=fontsize, fontsize_row=fontsize_row,fontsize_col=fontsize_col,cluster_rows = cluster_rows, cluster_cols = cluster_cols,scale=scale,
               display_numbers = matrix(ifelse(interaction_pvalue_CORRECTED < significance, "*", ""), nrow(interaction_pvalue_CORRECTED)))
      setHook("grid.newpage", NULL, "replace")
      grid.text(row_names, x=-0.03, rot=90, gp=gpar(fontsize=36))
      grid.text(col_names, y=0.01, gp=gpar(fontsize=36)  )
    }
    plot_pheatmap(mat, "Drug target" , cell_feature_name , tissue ) 
    dev.off()
  }
  return(mat)
}

significance_interaction_plot <- function(drug_feature_name,cell_feature_name,selection,cut_off,N = 1:length(f), database="GDSC",target_to_remove=target_to_remove, limit_row=0.95,limit_col=0.95) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_",database,"/",drug_feature_name,"_",cell_feature_name,"/",sep = "")
  path_permutation <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_",database,"/",drug_feature_name,"_",cell_feature_name,"_PERMUTATION/" )
  setwd(path) ; f <- list.files(path)
  for(i in N) {  #  i = 1
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    x <- x[ - target_to_remove ,  ] ; 
    interaction_mat <- x 
    tissue <- regmatches(f[i], regexpr('interaction.+?target', f[i])) ; tissue <- gsub('.{7}$', '', tissue) ; tissue <- substring(tissue, 13)
    
    filenames <- list.files(paste0(path_permutation,tissue,"/"), pattern="*COUNT")
    x <- read.csv(paste0(path_permutation,tissue,"/",filenames), row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    x <- x[ - target_to_remove ,  ] ; 
    COUNT_mat <- x 
    
    COUNT_mat_pos <- COUNT_mat ; COUNT_mat_pos[which(interaction_mat<0, arr.ind = T)] <- 0
    COUNT_mat_neg <- COUNT_mat ; COUNT_mat_neg[which(interaction_mat>0, arr.ind = T)] <- NA
    COUNT_mat_neg <- 1000 - COUNT_mat_neg ; COUNT_mat_neg[is.na(COUNT_mat_neg)] <- 0
  
    COUNT_mat <- COUNT_mat_pos + COUNT_mat_neg 
    pvalue_mat <- COUNT_mat/1000
    write.csv(pvalue_mat,paste0(path_permutation,tissue,"/interaction_pvalue"))
    pvalue_mat <- t( apply(pvalue_mat, 1, function(x) p.adjust(x, method = "BY" , n=length(colnames(pvalue_mat)) ) ) ) # for each target, correct for number of pathways, as this is what was randomized
    write.csv(pvalue_mat,paste0(path_permutation,tissue,"/interaction_pvalue_CORRECTED"))
  }
  return(pvalue_mat)
}


save_interaction_plot_quantile <- function(folder,drug_feature_name,cell_feature_name,selection,cut_off,N = 1:length(f) , limit_row=0.95,limit_col=0.95) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,"/",sep = "")
  setwd(path) ; f <- list.files(path)
  for(i in N) {
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    
    ##################### decide about target selection: conservative or all
    if(selection == "conservative") { x <- x[ - target_to_remove ,  ] ; mat <- x }
    if(selection == "conservative_include_BIRC5") { mat <- x[ - target_to_remove ,  ] ; mat <- rbind(mat,x["BIRC5", ]) }
    if(selection == "all_target") {  mat <- x }
    ##################### decide about target cut off: single top hit of sum across all
    if(cut_off=="absolute") {
      v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
      if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
        mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:25] ]
      } 
    }
    if(cut_off=="single") {
      mat <- subset_row_abs(mat, abs_limit = limit_row, obs=1) ;  v <- apply(mat,1,var) ; v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:20] , ]
      if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
        mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ; v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:20] ]
      } 
    }
    mat <- mat[complete.cases(mat),]
    mat <- quantile_normalisation(mat)  ####### NORMALIZATION
    tissue <- regmatches(f[i], regexpr('interaction.+?target', f[i])) ; tissue <- gsub('.{7}$', '', tissue) ; tissue <- substring(tissue, 13)
    pdf(paste("/Users/miyang/Documents/RWTH_Aachen/",folder,"/PLOTS/GDSC_",drug_feature_name,"_",cell_feature_name,"/",selection,'/',cut_off,'/',tissue,".pdf", sep = ""), width = 17 , height = 15, onefile = F ) # 13.5 15
    plot_pheatmap(mat,"Drug target",cell_feature_name,tissue) 
    dev.off()
  }
  return(mat)
}

retrieve_RANK_cell <- function(folder,drug_feature_name,cell_feature_name,cell_feature_selection,N = 1:length(f)) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,sep = "")
  setwd(path) ; f <- list.files(path)
  tissue_rank <- c()
  for(i in N) {
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    
    x_cell <- x[ ,cell_feature_selection] ; names(x_cell) <- rownames(x)
    x_cell <- x_cell[order(-abs(x_cell)) ]  
    tissue_rank <- cbind(tissue_rank, names(x_cell)[1:20] )
  }
  colnames(tissue_rank) <- tissue
  return(tissue_rank)
}

retrieve_RANK_drug <- function(folder,drug_feature_name,cell_feature_name,drug_feature_selection,N = 1:length(f)) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,sep = "")
  setwd(path) ; f <- list.files(path)
  tissue_rank <- c()
  for(i in N) {
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    
    x_drug <- x[ drug_feature_selection, ] ; names(x_drug) <- colnames(x)
    x_drug <- x_drug[order(-abs(x_drug)) ]  
    tissue_rank <- cbind(tissue_rank, names(x_drug)[1:20] )
  }
  colnames(tissue_rank) <- tissue
  return(tissue_rank)
}

retrieve_RANK_interaction <- function(folder,drug_feature_name,cell_feature_name,interaction_selection,N = 1:length(f)) {
  path <- paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/TISSUE_SPECIFIC_GDSC/",drug_feature_name,"_",cell_feature_name,sep = "")
  setwd(path) ; f <- list.files(path)
  tissue_rank <- c()
  for(i in N) {
    x <- read.csv(f[i], row.names=1)
    rownames(x)[grep("G9a and GLP methyltransferases",rownames(x))] <- "G9a and GLP"
    rownames(x)[grep("dsDNA break induction",rownames(x))] <- "dsDNA break"
    rownames(x)[grep("FNTA",rownames(x))] <- "FNTA"
    
    x_interaction <- x[ interaction_selection[1], interaction_selection[2] ]  
    x <- data.matrix(x) ; x <- as.numeric(x) ; x <- x[order(-abs(x))]
    tissue_rank <- rbind(tissue_rank, c(which(x==x_interaction), round(x_interaction,4) ) )
  }
  colnames(tissue_rank) <- c("rank","weight")
  rownames(tissue_rank) <- tissue
  return(tissue_rank)
}


retrieve_interaction_plot_single <- function(database="TISSUE_SPECIFIC_GDSC",drug_feature_name,cell_feature_name,selection,cut_off,N = 1:length(f), tissue_name , limit_row=0.95,limit_col=0.95) {
  path <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/",database,"/",drug_feature_name,"_",cell_feature_name,"/" )
  setwd(path) ; f <- list.files(path)
  mat_list <- list()
  pvalue <- list()
  for(t in N) {
    x <- read.csv(f[t], row.names=1)
    ##################### decide about target selection: conservative or all
    if(selection == "conservative") { x <- x[ - target_to_remove ,  ] ; mat <- x }
    if(selection == "all_target") {  mat <- x }
    ##################### decide about target cut off: single top hit of sum across all
    #     if(cut_off=="absolute") {
    #       v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
    #       if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
    #         mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:25] ]
    #       } 
    #     }
    #     if(cut_off=="single") {
    #       mat <- subset_row_abs(mat, abs_limit = limit_row, obs=1) ;  v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
    #       if(cell_feature_name %not in% c("progeny11","progeny14") ){ 
    #         mat <- subset_col_abs(mat, abs_limit = limit_col, obs=1) ;  v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ , names(v)[1:25] ]
    #       } 
    #     }
    mat_list[[t]] <- mat[complete.cases(mat),]
  }
  mat_list <- mat_list[lapply(mat_list,length)>0]
  
  return(mat_list )
}

retrieve_interaction_plot <- function(database="TISSUE_SPECIFIC_GDSC",drug_feature_name,cell_feature_name,selection,cut_off,N = 1:length(f), tissue_name) {
  path <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/",database,"/",drug_feature_name,"_",cell_feature_name,"/" )
  path_permutation <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/",database,"/",drug_feature_name,"_",cell_feature_name,"_PERMUTATION/" )
  setwd(path) ; f <- list.files(path)
  mat_list <- list()
  pvalue <- list()
  for(t in N) {
    x <- read.csv(f[t], row.names=1)
    ##################### decide about target selection: conservative or all
    if(selection == "conservative") { x <- x[ - target_to_remove ,  ] ; mat <- x }
    if(selection == "all_target") {  mat <- x }
    
    mat_list[[t]] <- mat[complete.cases(mat),]
    x <- read.csv(paste0(path_permutation,tissue_name[t],"/interaction_pvalue_CORRECTED"),row.names = 1)
    pvalue[[t]] <- x
  }
  mat_list <- mat_list[lapply(mat_list,length)>0]
  pvalue <- pvalue[lapply(pvalue,length)>0]
  return(list(mat_list,pvalue))
}



interaction_plot <- function(database="TISSUE_SPECIFIC_GDSC",drug_feature_name,cell_feature_name,selection,cut_off, t ) {
  path <- paste0("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/DATA_RESULT_STORAGE/",database,"/",drug_feature_name,"_",cell_feature_name,"/" )
  f <- list.files(path)
    x <- read.csv(paste0(path,f[t]), row.names=1)
    ##################### decide about target selection: conservative or all
    if(selection == "conservative") { x <- x[ - target_to_remove ,  ] ; mat <- x }
    if(selection == "all_target") {  mat <- x }
    ##################### decide about target cut off: single top hit of sum across all
    if(cut_off=="absolute") {
      v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
      v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; if(cell_feature_name!="progeny11"){ mat <- mat[ ,names(v)[1:14]] } 
    }
    if(cut_off=="single") {
      mat <- subset_row_abs(mat, abs_limit = 0.95, obs=1) ;  v <- apply(mat,1,var) ;  v <- v[ order(-abs(v)) ] ; mat <- mat[ names(v)[1:25] , ]
      v <- apply(mat,2,var) ;  v <- v[ order(-abs(v)) ] ; if(cell_feature_name %not in% c("progeny11","progeny14") ){ mat <- mat[ ,names(v)[1:14]] } 
    }
    mat <- mat[complete.cases(mat),]
    tissue <- regmatches(f[t], regexpr('interaction.+?target', f[t])) ; tissue <- gsub('.{7}$', '', tissue) ; tissue <- substring(tissue, 13)
    plot_pheatmap(mat, paste("Drug ",drug_feature_name,sep = "") , cell_feature_name , paste(drug_feature_name," - ",cell_feature_name," (",tissue,")",sep = "")  ) 
}


# Generate top hits drug target
# pdf(paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/PLOTS/GDSC_Target_Pathway/",group,'/',cut_off,'/',tissue,"_top_hits_corrplot" ,".pdf", sep = ""), width = 11 , height = 14, onefile = F )
# L <- target_combo_pathway (mat, top_association = 5 , corr_diff = 0.3 )
# library(corrplot) ; corrplot( L[[1]] , order = "hclust")
# dev.off()
# sga_hits <- L[[2]]
# require(gridExtra)
# pdf(paste("/Users/miyang/Documents/RWTH_Aachen/MACAU_PROJECT/PLOTS/GDSC_Target_Pathway/",group,'/',cut_off,'/',tissue,"_top_hits_TABLE" ,".pdf", sep = ""), width = 11 , height = 14, onefile = T )
# total_rows_per_page = 40
# start_row = 1 
# if(total_rows_per_page > nrow(sga_hits)){
#   end_row = nrow(sga_hits)
# }  else  {  end_row = total_rows_per_page  }    
# for(i in 1:ceiling(nrow(sga_hits)/total_rows_per_page)){
#   grid.newpage()   
#   grid.table(sga_hits[start_row:end_row, ] )
#   start_row = end_row + 1
#   if((total_rows_per_page + end_row) < nrow(sga_hits)){
#     end_row = total_rows_per_page + end_row
#   } else {  end_row = nrow(sga_hits)  }    
# }
# dev.off()