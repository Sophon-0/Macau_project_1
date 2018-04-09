path <- "~/Documents/RWTH_Aachen"
source(paste(path,"/FUNCTIONS/general_functions.R", sep=""))
load(paste(path,"/SANGER_DATA/MASTER_LIST_22112013.ro", sep="")) 
DRUG_ANALYSIS_SET <- read.csv("~/Documents/RWTH_Aachen/SANGER_DATA/DRUG_ANALYSIS_SET", row.names = 1)

load("~/Documents/RWTH_Aachen/SANGER_DATA/v17a_IC50s.Rdata") ; IC50 <- t(v17a_IC50s); IC50 <- convert_drugID(IC50) 
load(paste(path, "/SANGER_DATA/FEATURES/GEX_ALL.Rdata", sep=""))
tissue_label_gdsc <- read.csv("~/Documents/RWTH_Aachen/SANGER_DATA/tissue_label_gdsc_ID", row.names=1)
tissue_label <- tissue_label_gdsc ; table(tissue_label$tissue)

features <- GEX_ALL
response <- t(IC50)

load("~/Documents/RWTH_Aachen/SANGER_DATA/PANCAN_simple_MOBEM.rdata")
x <- MoBEM ; x <- x[ -grep("HypMET", rownames(x)) , ] ; mutation <- t(x)
mutation <- mutation[ order(rownames(mutation)), ] 
write.csv(mutation, "~/Documents/RWTH_Aachen/SANGER_DATA/FEATURES/SNP_CNV")

source( paste(path, "/TUTORIALS/PathActivities_Dorothea_with_methylation/src/TFactivities/code/lib_SLEA.r", sep = "" ) ) 
load(paste(path, "/TUTORIALS/PathActivities_Dorothea_with_methylation/src/TFactivities/regulons/merged_062016.rdata", sep = "" ) )

NES_GDSC_GSVA <- read.csv("/Users/miyang/Documents/RWTH_Aachen/SANGER_DATA/FEATURES/NES_GDSC_GSVA",row.names = 1, check.names = F)
NES_GDSC_VIPER <- read.csv("/Users/miyang/Documents/RWTH_Aachen/SANGER_DATA/FEATURES/NES_GDSC_VIPER",row.names = 1, check.names = F)

######### prepare subset of gene, selected by prior knowledge. 
target <- read.csv("/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA/target", check.names = F) ; drug_names <- target[ ,1] ; target <- target[ ,-1]
weights_rank_GEX_IC50_50ite <- read.csv("~/Documents/RWTH_Aachen/SANGER_DATA/DATA_RESULT_STORAGE_IC50/weights_rank_GEX_IC50_50ite", row.names=1)
weights_rank_GEX_IC50_50ite <- t(weights_rank_GEX_IC50_50ite)
v <- as.character ( weights_rank_GEX_IC50_50ite[ , 1 ] )  
relevant_gene <- unique(c(colnames(target), v ))

##################################### create GEX/progeny from tissue #####################################
###################################### only if more than 20 samples ######################################

result_folder <-  "/Users/miyang/Documents/RWTH_Aachen/SANGER_DATA/TISSUE"

retrieve_tissue (name="aero_dig_tract",tissue_label, features,mutation, response, tissue=c("aero_dig_tract") )
retrieve_tissue (name="bone",tissue_label, features,mutation, response, tissue=c("bone"))
retrieve_tissue (name="brain",tissue_label, features,mutation, response, tissue=c("nervous_system", "neuroblastoma"))
retrieve_tissue (name="breast",tissue_label, features,mutation, response, tissue=c("breast"))
retrieve_tissue (name="colon",tissue_label, features,mutation, response,tissue=c("large_intestine"))
retrieve_tissue (name="kidney",tissue_label, features,mutation, response, tissue=c("kidney"))
retrieve_tissue (name="leukemia",tissue_label, features,mutation, response, tissue=c("leukemia"))
retrieve_tissue (name="liver",tissue_label, features,mutation, response, tissue=c("liver"))  
retrieve_tissue (name="lung_NSCLC",tissue_label, features,mutation, response, tissue=c("lung_NSCLC"))
retrieve_tissue (name="lung_SCLC",tissue_label, features,mutation, response, tissue=c("lung_SCLC"))
retrieve_tissue (name="lymphoma",tissue_label, features,mutation, response, tissue=c("lymphoma"))
retrieve_tissue (name="ovary", tissue_label, features,mutation, response, tissue=c("ovary"))
retrieve_tissue (name="pancreas",tissue_label, features,mutation, response, tissue=c("pancreas"))
retrieve_tissue (name="skin",tissue_label, features,mutation, response, tissue=c("skin"))
retrieve_tissue (name="soft_tissue",tissue_label, features,mutation, response, tissue=c("soft_tissue"))
retrieve_tissue (name="stomach",tissue_label, features,mutation, response, tissue=c("stomach"))

##################################### create GEX R object from tissue ####################################
##########################################################################################################
x <- GEX_ALL

GEX_tissue <- list()
for(i in 1:length(table(tissue_label$tissue))) {
  tissue_ID <- tissue_label[ which(tissue_label$tissue == names(table(tissue_label$tissue))[i]) , 1]
  t <- x[ rownames(x) %in% tissue_ID , ] ; t <- t[order(rownames(t)) , ]
  GEX_tissue[[i]] <- t
}
names(GEX_tissue) <- names(table(tissue_label$tissue))
save(GEX_tissue,file = paste("~/Documents/RWTH_Aachen/SANGER_DATA/FEATURES/GEX_tissue.Rdata",sep=""))


################################## create SNP + CNV R object from tissue #################################
##########################################################################################################
load("~/Documents/RWTH_Aachen/SANGER_DATA/PANCAN_simple_MOBEM.rdata")
x <- MoBEM ; x <- x[ -grep("HypMET", rownames(x)) , ] ; x <- t(x)

MoBEM_tissue <- list()
for(i in 1:length(table(tissue_label$tissue))) {
  tissue_ID <- tissue_label[ which(tissue_label$tissue == names(table(tissue_label$tissue))[i]) , 1]
  t <- x[ rownames(x) %in% tissue_ID , ] ; t <- t[order(rownames(t)) , ]
  MoBEM_tissue[[i]] <- t
}
names(MoBEM_tissue) <- names(table(tissue_label$tissue))
save(MoBEM_tissue,file = paste("~/Documents/RWTH_Aachen/SANGER_DATA/FEATURES/MoBEM_tissue.Rdata",sep=""))


##################################### create features of tissue type #####################################
##########################################################################################################
tissue <- create_features_dummy (MASTER_LIST, "gdsc_desc_3") [[1]]
storedata(tissue, "~/Documents/RWTH_Aachen/SANGER_DATA/FEATURES")
storedata(tissue, "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA")

common <- intersect(rownames(GEX_ALL), rownames(tissue))
GEX_ALL <- GEX_ALL[common, ] ; tissue <- tissue[common, ]
GEX_tissue <- cbind(tissue , GEX_ALL)

common <- intersect(colnames(IC50) , rownames(GEX_tissue))
IC50 <- IC50[ , common ] ; GEX_tissue <- GEX_tissue[common, ]

GEX_tissue <- GEX_tissue[order(rownames(GEX_tissue)),  ]
IC50 <- IC50[ , order(colnames(IC50))   ]

IC50_tissue <- IC50

storedata(IC50_tissue, "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA")
storedata(GEX_tissue, "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA")

######################################### REGRESSING tissue out  #########################################
##########################################################################################################

common <- intersect(rownames(GEX_ALL), rownames(MASTER_LIST))
GEX_ALL <- GEX_ALL[common, ] ; tissue <- MASTER_LIST[common,"gdsc_desc_3"]
GEX_tissue <- cbind(tissue , GEX_ALL)

common <- intersect(colnames(IC50) , rownames(GEX_tissue))
IC50 <- IC50[ , common ] ; GEX_tissue <- GEX_tissue[common, ]

GEX_tissue_catego <- GEX_tissue[order(rownames(GEX_tissue)),  ]

options(contrasts = c("contr.treatment", "contr.poly"))
for(i in 2:length(colnames(GEX_tissue_catego))) {
  GEX_tissue_catego[,i] <- resid(m1 <- lm( I(GEX_tissue_catego[,i]) ~ tissue, data = GEX_tissue_catego))
}

GEX_tissue_OUT <- GEX_tissue_catego[ ,-1]
storedata(GEX_tissue_OUT, "/Users/miyang/Documents/RWTH_Aachen/macau_work_dir/macau_test_sanger/DATA")


