######################## Pulmmonary carcinoid #################
# This repository was script in  Macintosh HD⁩ ▸ ⁨Utilisateurs⁩ ▸ ⁨mathian⁩ ▸ ⁨Bureau⁩ ▸ ⁨INSA_4_2⁩ ▸ ⁨STGE_CIRC⁩ ▸ 
# Figure 1 : 
# LNEN MOFA MAP 
# Attributes to include 
# S1 : Sample Overview -> And ML predictions -> and Survival 
# S4 : Mutation 
# Immune check point with VST Table
# retinoid and xenobiotic metabolism pathways : CYP and UGT
# Cancer revelant core genes
# IHC Data : LAMP3 et CD1 ? NO
# Homeobox genes : HNF1A and HNFA4
# Expresison of DLL3 and ASCL1 NOTCH1-3 -MKI67 - EIF1AX  Cf:  fig 5B , fig 5C and fig 6 (very important)

#######################
#     LIBRAIRIES      #
#######################
library(openxlsx)
library(data.table)
library(dplyr)
 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("bumphunter")

library(bumphunter)
library(minfi)

#######################
# IMPORTATION OF DATA #
#######################
Sample_overview <- read.xlsx("../SupplementaryTables_R1_20190318.xlsx", sheet = 1, startRow = 41, colNames = TRUE,
          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
          skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
          namedRegion = NULL, na.strings = "NA")

Somatic_mutation  <- read.xlsx("../SupplementaryTables_R1_20190318.xlsx", sheet = 4, startRow = 51, colNames = TRUE,
                              rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                              skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
                              namedRegion = NULL, na.strings = "NA")


Data_vst_50 <- read.table("../data/VST_nosex_50pc_TCACLCNECSCLC.txt",  sep = " ", dec="." , header = TRUE,   quote="")
Data_vst_all <- read.table("../data/VST_nosex_TCACLCNECSCLC.txt",  sep = " ", dec="." , header = TRUE,   quote="")

Ref_gene <- read.table("../data/VST_nosex_50pc_TCACLCNECSCLC_annot.txt",  sep = " ", dec="." , header =FALSE,   quote="")
Ref_gene_all <- read.table("../data/VST_nosex_TCACLCNECSCLC_annot.txt",  sep = " ", dec="." , header =FALSE,   quote="")

Methylation_data <- load('../methylation_final_LM.RData')
###########################
# Coordinates MOFA fig 1  #
###########################

# Summary
# ========
table(Sample_overview$Histopathology)
sum(is.na(Sample_overview$LF1.LNEN) ==F)
summary(Sample_overview$LF1.LNEN)
table(Sample_overview$Histopathology[is.na(Sample_overview$LF1.LNEN) ==F])
summary(Sample_overview$LF1.LNEN_SCLC)
sum(is.na(Sample_overview$LF1.LNEN_SCLC) ==F)
table(Sample_overview$Histopathology[is.na(Sample_overview$LF1.LNEN_SCLC) ==F])
table(Sample_overview$Epic.850K)

# Coords
# ======

Coords_MOFA_fig1 <- data.frame("Sample_ID" = Sample_overview$Sample_ID , "Axis1" = Sample_overview$LF1.LNEN , "Axis2" = Sample_overview$LF2.LNEN *-1)
Coords_MOFA_fig1 <- Coords_MOFA_fig1[complete.cases(Coords_MOFA_fig1),]
plot(Sample_overview$LF1.LNEN, Sample_overview$LF2.LNEN , col = as.factor(Sample_overview$Histopathology ))

# MOFA Sample ID
# ==============

mofa_lnen_sample_id = data.frame("Sample_ID"=Coords_MOFA_fig1$Sample_ID) # Complete case of mofa coordinates
All_sample_id = data.frame("Sample_ID"=Sample_overview$Sample_ID) 
###########################
# Attributes               #
###########################

# Clinical attributes
# ====================
Attributes_from_overview <- data.frame("Sample_ID" = Sample_overview$Sample_ID  ,"Histopathology" = Sample_overview$Histopathology , "Stage_UICC" = Sample_overview$Stage_UICC , "Age"= Sample_overview$Age , "Age_class" = Sample_overview$Age_class , 
                                      "Sex" = Sample_overview$Sex , "Smoking_status" = Sample_overview$Smoking_status , "Professional_Asbestos_exposure" = Sample_overview$Professional_exposure , "Survival_months" = Sample_overview$Survival_months,
                                       "Neutrophil.to.Lymphocyte_ratio" = Sample_overview$Neutrophil.to.Lymphocyte_ratio  , "Cluster_LNEN" =Sample_overview$cluster_LNEN , "Cluster_LNET" =Sample_overview$cluster_LNET)

Attributes_from_overview <- cbind(Attributes_from_overview , Sample_overview[ , 42:52])

# Later
##  Attributes_from_overview <- merge(Attributes_from_overview , mofa_lnen_sample_id, by="Sample_ID")
  
# Genes of interest
# =================

# Immune checkpoint and ligand receptor Fig2E
# ____________________________________________

Embl_lag3 = as.character(Ref_gene$V1[Ref_gene$V7 == "LAG3"])
Embl_IGSF11 = as.character(Ref_gene$V1[Ref_gene$V7 == "IGSF11" ])
# VISTA VSIR : https://www.uniprot.org/uniprot/Q9H7M9
Embl_VISTA = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "C10orf54" ])
Embl_PDCD1LG2 = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "PDCD1LG2" ]) # Becarful not in 50pc
Embl_LGALS9  = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "LGALS9" ]) # Becarful not in 50pc
Embl_CD274 =  as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "CD274"]) # Becarful not in 50pc
Embl_PDCD1 = as.character(Ref_gene$V1[Ref_gene$V7 == "PDCD1" ])
Embl_HAVCR2 = as.character(Ref_gene$V1[Ref_gene$V7 == "HAVCR2" ])
Embl_CD86 = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "CD86"]) # Becarful not in 50pc
Embl_CD80 = as.character(Ref_gene$V1[Ref_gene$V7 == "CD80" ])
Embl_CTLA4 = as.character(Ref_gene$V1[Ref_gene$V7 == "CTLA4" ])



# HLA-D Sum of 20 genes or 20 attributes : (?)
index_HLA_D <- c( grep("HLA-D", Ref_gene$V7 ))
HLA_D <-Ref_gene$V1[c(index_HLA_D )] 
names(HLA_D)<-Ref_gene$V7[c(index_HLA_D )] 

# Clean data frame : 
t_Data_vst_all <-as.data.frame(t( Data_vst_all ))
t_Data_vst_all <- setDT(t_Data_vst_all , keep.rownames = TRUE)[]
colnames(t_Data_vst_all)[1] <- "Sample_ID"
t_Data_vst_all <- as.data.frame(t_Data_vst_all)
Data_vst_all_with_sample <- merge(t_Data_vst_all , All_sample_id, by='Sample_ID')
# Rmq : Only 158 samples with RNA seq data
Sample_id_rna_seq = Data_vst_all_with_sample$Sample_ID

# Name of genes of interest :
gene_interest_names <- c("LAG3" , "IGSF11" , "VISTA_VSIR" , "PDCD1LG2" , "LGALS9" , "CD274" ,
                   "PDCD1" , "HAVCR2" , "CD86" , "CD80" , "CTLA4", names(HLA_D) )
gene_interest_embl <- c(Embl_lag3 , Embl_IGSF11,Embl_VISTA, Embl_PDCD1LG2, Embl_LGALS9, Embl_CD274, Embl_PDCD1, 
                        Embl_HAVCR2, Embl_CD86, Embl_CD80, Embl_CTLA4 )
for (i in 1:12){
  gene_interest_embl <- append(as.character(HLA_D[[i]]) ,gene_interest_embl ) 
}

# Data frame :
gene_interest_fig2E <- data.frame("Sample_ID" = Sample_id_rna_seq)
for (i in 1:23){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl[i]))
  #print(as.numeric(n_col))
  gene_name <- as.character(gene_interest_names[i])
  #print(gene_name)
  gene_interest_fig2E[gene_name] <- Data_vst_all_with_sample[,n_col ]
}

HLA_D_sum = rowSums(gene_interest_fig2E[,13:24], na.rm = TRUE)
HLA_D_mean= apply(gene_interest_fig2E[,13:24], 1 , mean)
# Mutation Fig3A:
length(unique(Somatic_mutation$Gene.Symbol))
mutation_matrix =matrix(ncol = length(unique(Somatic_mutation$Gene.Symbol)) , nrow= length(mofa_lnen_sample_id$Sample_ID))
colnames(mutation_matrix)<- unique(Somatic_mutation$Gene.Symbol)
rownames(mutation_matrix) <- mofa_lnen_sample_id$Sample_ID

for (i in 1:dim(Somatic_mutation)[1]){
  if(identical(which(rownames(mutation_matrix)== Somatic_mutation$Sample_ID[i]),  integer(0))  == F){
    n_row = which(rownames(mutation_matrix)== Somatic_mutation$Sample_ID[i])
    n_col = which(colnames(mutation_matrix)== Somatic_mutation$Gene.Symbol[i])
    mutation_matrix[n_row,n_col] = 1
  }
}


sort( colSums(mutation_matrix, na.rm = TRUE, dims = 1) , decreasing = TRUE  )[1:20]

fig3A_index_col_mutation_matrix = c(which(colnames(mutation_matrix)== "MEN1"),which(colnames(mutation_matrix)== "ARID1A"),
                                    which(colnames(mutation_matrix)== "ARID2"),   which(colnames(mutation_matrix)== "ATM"),
                                    which(colnames(mutation_matrix)== "DNAH17"), which(colnames(mutation_matrix)== "DOT1L"),
                                    which(colnames(mutation_matrix)== "HERC2"),  which(colnames(mutation_matrix)== "JMJD1C"),
                                    which(colnames(mutation_matrix)== "KDM5C"),which(colnames(mutation_matrix)== "NIPBL"),
                                    which(colnames(mutation_matrix)== "PALMD"), which(colnames(mutation_matrix)== "PSIP1"),
                                    which(colnames(mutation_matrix)== "RLIM"), which(colnames(mutation_matrix)== "ROBO1"),
                                    which(colnames(mutation_matrix)== "SEC31A"), which(colnames(mutation_matrix)== "SELP"),
                                    which(colnames(mutation_matrix)== "SMARCA1"),which(colnames(mutation_matrix)== "SMARCA2"),
                                    which(colnames(mutation_matrix)== "TP53"),which(colnames(mutation_matrix)== "WDR26"),
                                    which(colnames(mutation_matrix)== "BAP1"), which(colnames(mutation_matrix)== "RB1"))

mutation_matrix = mutation_matrix[,c(fig3A_index_col_mutation_matrix )]
New_col_name = c()
for (i in 1:dim(mutation_matrix)[2]){
  New_col_name[i] = paste("Mutation_" , colnames(mutation_matrix)[i], sep="")
}
colnames(mutation_matrix) <- New_col_name
# Doit-on tous les inclure ?

# Genes Fig5A :
Embl_ASCL1 =as.character(Ref_gene$V1[Ref_gene$V7 == "ASCL1"])
Embl_DLL3= as.character(Ref_gene$V1[Ref_gene$V7 == "DLL3"])
Embl_SLIT1 = as.character(Ref_gene$V1[Ref_gene$V7 == "SLIT1"])
Embl_ROBO1 = as.character(Ref_gene$V1[Ref_gene$V7 == "ROBO1"])
Embl_ANGPTL3 = as.character(Ref_gene$V1[Ref_gene$V7 == "ANGPTL3"])
Embl_ERBB4 = as.character(Ref_gene$V1[Ref_gene$V7 == "ERBB4"])
Embl_OTP = as.character(Ref_gene$V1[Ref_gene$V7 == "OTP"])
Embl_NKX2_1 = as.character(Ref_gene$V1[Ref_gene$V7 == "NKX2-1"])

gene_interest_names_fi5A <- c("ASCL1","DLL3","SLIT1" , "ROBO1","ANGPTL3","ERBB4","OTP", "NKX2-1" )
gene_interest_embl_fig5A <- c(Embl_ASCL1,Embl_DLL3,Embl_SLIT1,Embl_ROBO1,Embl_ANGPTL3,Embl_ERBB4,Embl_OTP,Embl_NKX2_1)
gene_interest_fig5A <- data.frame("Sample_ID" = Sample_id_rna_seq)
for (i in 1:8){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_fig5A[i]))
  gene_name <- as.character(gene_interest_names_fi5A[i])
  gene_interest_fig5A[gene_name] <- Data_vst_all_with_sample[,n_col ]
}

# Machine Learning Fig1 :

ML_prediction <- function(data){
  res_ml = c()
  c_methyl = 0
  for (i in 1:dim(data)[1]){
    val <- data[i,]
    names(val)<- c("Atypical", "Typical", "LCNEC")
    if (sum(is.na(val)) ==3){
       res_ml[i] = NA
    }
    else {
      s_val = (sort(val, decreasing = TRUE))
      #print(s_val)
      c_methyl = c_methyl + 1
      if (s_val[1]/s_val[2] < 1.5){
        res_ml[i] = "Unclassified"
      }
      else{
        #print(names(s_val[1]))
        res_ml[i] = names(s_val[1])
      }
    }
  }
  return(res_ml)
}
ML_methyl <- ML_prediction(Sample_overview[,53:55])
ML_expr <- ML_prediction(Sample_overview[,56:58])
ML_MKI67 <- ML_prediction(Sample_overview[,59:61]) # Fig S11
ML_Mofa <- ML_prediction(Sample_overview[,62:64])
ML_expr_methyl <- ML_prediction( Sample_overview[,65:67])  # Fig S9

table(ML_Mofa)
# Cobfusion matrix of ML and Mofa sample :

Histpopathology_4_classes = Sample_overview$Histopathology
table(Histpopathology_4_classes)
for (i in 1:length(Histpopathology_4_classes )){
  if ( Histpopathology_4_classes[i]=="LCNEC combined SCLC" ){
    Histpopathology_4_classes[i] = "LCNEC"
  }
  else if (  Histpopathology_4_classes[i]== "LCNEC combined SCLC"){
    Histpopathology_4_classes[i] = "LCNEC"
  }
  else if (Histpopathology_4_classes[i]== "LCNEC combined SqCC"){
    Histpopathology_4_classes[i] = "LCNEC"
  }
  else if (Histpopathology_4_classes[i]== "LCNEC combined ADC"){
    Histpopathology_4_classes[i] = "LCNEC"
  }
}

table( Histpopathology_4_classes)


# Confusion Matrix
# -------------------

# For MOFA Laten factor

MOFA_ML = data.frame("Pred_MOFA"= ML_Mofa, "Sample_ID"  = Sample_overview$Sample_ID)
merge_pred_mofa = merge(HISTO_df, MOFA_ML, by="Sample_ID")

tab_conf =  table(merge_pred_mofa$Histpopathology_4_classes , merge_pred_mofa$Pred_MOFA)
prop.table(tab_conf, 2)
#Methyl <- load("../methylation_final_LM.RData")


# For Expr 
# --------

HISTO_df  = data.frame("Histpopathology_4_classes"= Histpopathology_4_classes, "Sample_ID"  = Sample_overview$Sample_ID)
EXPR_ML = data.frame("Expr_ML"= ML_expr, "Sample_ID"  = Sample_overview$Sample_ID)
merge_pred_expr = merge(HISTO_df, EXPR_ML, by="Sample_ID")




tab_conf =  table(merge_pred_mofa$Histpopathology_4_classes , merge_pred_expr$Expr_ML)
prop.table(tab_conf, 1)  # OK 


# ML Type Fig
# ------------

ML_histopatholical_type <- data.frame("ML_methyl" = ML_methyl  ,"ML_expr" = ML_expr , stringsAsFactors=F) 
Res_type_ml = c()

for (i in 1:dim(ML_histopatholical_type)[1])
  {
  # Aucun NA  sur la ligne et les deux types predits sont différents et ne sont pas de la catégorie unclassified -> Précdiction Unclassified
  if(sum(is.na( ML_histopatholical_type[i,]) == F) == 2 & 
     ML_histopatholical_type$ML_methyl[i] != ML_histopatholical_type$ML_expr[i] &
     ML_histopatholical_type$ML_methyl[i] != "Unclassified"  &  
     ML_histopatholical_type$ML_expr[i] != "Unclassified" )
  {
    Res_type_ml[i] = "Unclassified"
  }
  # Si 1 NA le type prédit est celui restant
  else if (sum(is.na( ML_histopatholical_type[i,]) == F) == 1 ) 
  {
    Res_type_ml[i] =as.character( ML_histopatholical_type[i, which(is.na(ML_histopatholical_type[i,])==F)] )
  }
  # Si Aucun  NA le type prédit est au  un Unclassified -> Type restant si deux unclassied -> Type Un classifies
  else if (sum(is.na( ML_histopatholical_type[i,]) == F) == 2  )
    { 
    if (ML_histopatholical_type$ML_methyl[i] == "Unclassified" &
        ML_histopatholical_type$ML_expr[i] == "Unclassified")
      {
      Res_type_ml[i] = "Unclassified"
      }
    else
    {
      Res_type_ml[i] =as.character( ML_histopatholical_type[i, which(ML_histopatholical_type[i,]!="Unclassified") ] )
    }
  }
  # Si Aucun NA et deux et Type identiques -> On conserve le mm
  else if (sum(is.na( ML_histopatholical_type[i,]) == F) == 2 &  
           ML_histopatholical_type$ML_methyl[i]== ML_histopatholical_type$ML_expr[i] )
    {
    Res_type_ml[i] = as.character( ML_histopatholical_type[i, 1 ]  )
  }
  # Si deux NA -> NA
  else if(is.na(ML_histopatholical_type$ML_methyl[i]) == T & is.na(ML_histopatholical_type$ML_expr[i])==T)
    {
    Res_type_ml[i] =NA
  }
  
}

Res_type_ml

# Confusion matrix  fig1B
# ------------------------

table_conf_f1B = table(merge_pred_mofa$Histpopathology_4_classes ,Res_type_ml)
prop.table(table_conf_f1B,1)


# Figure 4B 
# ----------
Embl_CD1A =as.character(Ref_gene$V1[Ref_gene$V7 == "CD1A"])
Embl_LAMP3 =as.character(Ref_gene$V1[Ref_gene$V7 == "LAMP3"])

gene_interest_names_fi4B <- c("CD1A",  "LAMP3" )
gene_interest_embl_fig4B <- c(Embl_CD1A , Embl_LAMP3)
gene_interest_fig4B <- data.frame("Sample_ID" =Sample_id_rna_seq)
for (i in 1:2){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_fig4B[i]))
  gene_name <- as.character(gene_interest_names_fi4B[i])
  gene_interest_fig4B[gene_interest_names_fi4B ] <- Data_vst_all_with_sample[,n_col ]
}



# Retinoid and Xenobiotic
# -----------------------
Embl_CYP2C8=as.character(Ref_gene$V1[Ref_gene$V7 == "CYP2C8"])
Embl_CYP2C9=as.character(Ref_gene$V1[Ref_gene$V7 == "CYP2C9"])
Embl_CYP2C19=as.character(Ref_gene$V1[Ref_gene$V7 == "CYP2C19"])
Embl_CYP3A5=as.character(Ref_gene$V1[Ref_gene$V7 == "CYP3A5"])
Embl_CYP39A1=as.character(Ref_gene$V1[Ref_gene$V7 == "CYP39A1"])
Embl_UGT2A3 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2A3"])
Embl_UGT2B4 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B4"])
Embl_UGT2B7 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B7"])
Embl_UGT2B11 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B11"])
Embl_UGT2B15 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B15"])
Embl_UGT2B17 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B17"])
Embl_UGT2B10 = as.character(Ref_gene$V1[Ref_gene$V7 == "UGT2B10"])


gene_interest_names_fi4D <- c("CYP2C8", "CYP2C9", "CYP2C19", "CYP3A5", "UGT2A3","UGT2B4", "UGT2B7", "UGT2B11","UGT2B15", "UGT2B17" )
gene_interest_embl_fig4D <- c(Embl_CYP2C8,Embl_CYP2C9, Embl_CYP2C19,Embl_CYP3A5, Embl_UGT2A3, Embl_UGT2B4, Embl_UGT2B7, Embl_UGT2B11, Embl_UGT2B15 , Embl_UGT2B17 )
gene_interest_fig4D <- data.frame("Sample_ID" =Sample_id_rna_seq)
for (i in 1:10){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_fig4D[i]))
  gene_name <- as.character(gene_interest_names_fi4D[i])
  gene_interest_fig4D[gene_interest_names_fi4D ] <- Data_vst_all_with_sample[,n_col ]
}

# Gene expression Fig 6
# ----------------------
Embl_MKI67 = as.character(Ref_gene$V1[Ref_gene$V7 == "MKI67"])
gene_interest_names_fig6 <- c("MKI67" )
gene_interest_embl_fig6 <- c(Embl_MKI67 )
gene_interest_fig6 <- data.frame("Sample_ID" =Sample_id_rna_seq)
for (i in 1:1){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_fig6[i]))
  gene_name <- as.character(gene_interest_names_fig6[i])
  gene_interest_fig6[gene_interest_names_fig6 ] <- Data_vst_all_with_sample[,n_col ]
}


#############################
# CHEMOkINES               #
#############################
Embl_CCL2=as.character(Ref_gene$V1[Ref_gene$V7 == "CCL2"])
Embl_CCL7=as.character(Ref_gene$V1[Ref_gene$V7 == "CCL7"])
Embl_CCL19=as.character(Ref_gene$V1[Ref_gene$V7 == "CCL19"])
Embl_CCL21=as.character(Ref_gene$V1[Ref_gene$V7 == "CCL21"])
Embl_CCL22=as.character(Ref_gene$V1[Ref_gene$V7 == "CCL22"])
Embl_IL8=as.character(Ref_gene$V1[Ref_gene$V7 == "IL8"])
Embl_CXCL1=as.character(Ref_gene$V1[Ref_gene$V7 == "CXCL1"])
Embl_CXCL3=as.character(Ref_gene$V1[Ref_gene$V7 == "CXCL3"])
Embl_CXCL5=as.character(Ref_gene$V1[Ref_gene$V7 == "CXCL5"])

gene_interest_names_chemokines <- c("CCL2","CCL7" , "CCL19",  "CCL21","CCL22", "IL8","CXCL1", "CXCL3","CXCL5")
gene_interest_embl_chemokines<- c(Embl_CCL2, Embl_CCL7, Embl_CCL19 ,Embl_CCL21, Embl_CCL22, Embl_IL8, Embl_CXCL1, Embl_CXCL3, Embl_CXCL5)
gene_interest_chemokines <- data.frame("Sample_ID" =Sample_id_rna_seq)
for (i in 1:9){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_chemokines[i]))
  gene_name <- as.character(gene_interest_names_chemokines[i])
  gene_interest_chemokines[gene_interest_names_chemokines ] <- Data_vst_all_with_sample[,n_col ]
}

New_col_name = c()
for (i in 2:10){
  New_col_name[i-1] = paste("Chemokine_" , colnames(gene_interest_chemokines)[i], sep="")
}
colnames(gene_interest_chemokines)[2:10] <- New_col_name

# Fig S24
# __________

Embl_NOTCH1=as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "NOTCH1"])
Embl_NOTCH2=as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "NOTCH2"])
Embl_NOTCH3=as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "NOTCH3"])
Embl_NOTCH4=as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "NOTCH4"])
gene_interest_names_notch <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
gene_interest_embl_notch<- c(Embl_NOTCH1, Embl_NOTCH2, Embl_NOTCH3, Embl_NOTCH4)
gene_interest_notch<- data.frame("Sample_ID" =Sample_id_rna_seq)
for (i in 1:4){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl_notch[i]))
  gene_name <- as.character(gene_interest_names_notch[i])
  gene_interest_notch[gene_interest_names_notch ] <- Data_vst_all_with_sample[,n_col ]
}


# Fig 5C
# Fig S23 -> Methylation

# Metylation 
# ----------
metadata=pData(funnometa) # metadata
#Mdata = minfi::getM(funnometa)
#colnames(Mdata)=sapply(colnames(Mdata),function(i) metadata$Ms_id[which(metadata$barcode==i)])
#minfi::getMeth(funnometa)


#############################
# MERGE ALL ATTRIBUTES      #
#############################

dim(Attributes_from_overview)
dim(mutation_matrix)
Mutation_df <- as.data.frame(mutation_matrix, rownames = T , colnames = T)
Mutation_df  <- setDT(Mutation_df  , keep.rownames = TRUE)[]
colnames(Mutation_df )[1] <- "Sample_ID"
head(mutation_matrix,3)

length(ML_methyl)
length(ML_expr)
length(ML_MKI67)
length(ML_Mofa )
length(ML_expr_methyl)
length(Res_type_ml) # Fig 1A pred

dim(gene_interest_fig5A)
dim(gene_interest_fig2E)
length(HLA_D_mean)
dim(gene_interest_fig4B)
dim(gene_interest_fig4D)
dim(gene_interest_fig6)
dim(gene_interest_chemokines)
dim(gene_interest_notch)
# Merge Genes Expr
# ----------------
gene_interest = merge(gene_interest_fig5A ,gene_interest_fig2E , by= "Sample_ID")
gene_interest = merge(gene_interest , gene_interest_fig4B , by= "Sample_ID")
gene_interest = merge(gene_interest , gene_interest_fig4D , by= "Sample_ID")
gene_interest = merge(gene_interest , gene_interest_fig6 , by= "Sample_ID")
gene_interest = merge(gene_interest , gene_interest_notch , by= "Sample_ID")
gene_interest = cbind(gene_interest , HLA_D_mean)
gene_interest = merge(gene_interest ,gene_interest_chemokines, by= "Sample_ID")


# Data frame ML
# -------------

ML_prediction_df = data.frame("Sample_ID" = Sample_overview$Sample_ID ,"ML_predictions_Methylation_Data" =ML_methyl , "ML_predictions_Expression_Data" = ML_expr , "ML_prediction_MKI67_Data"= ML_MKI67,
                              "ML_predictions_Mofa_Data"= ML_Mofa, "ML_predictions_Expression_Methylation_Data"= ML_expr_methyl, "ML_predictions_fig1"= Res_type_ml )

# Merging of iffrent attributes :
# -------------------------------

# Add Supra carcinoid class

HISTO_df  = data.frame("Histpopathology_4_classes"= as.character(Histpopathology_4_classes), "Sample_ID"  = as.character(Sample_overview$Sample_ID), stringsAsFactors = F)
supra_carcinoid_sample_id = c("LNEN005", "LNEN022", "LNEN012" , "S01522" , "S01513" , "LNEN021")
for (i in 1:258){
  if (as.character(HISTO_df$Sample_ID[i] ) %in% supra_carcinoid_sample_id ){
    HISTO_df$Histpopathology_4_classes[i] = "Supra_carcinoids"
  }
}
HISTO_df$Histpopathology_4_classes[HISTO_df$Sample_ID %in% supra_carcinoid_sample_id] ="Supra_carcinoid"
table(HISTO_df$Histpopathology_4_classes)

Histpopathology_simplified<- HISTO_df$Histpopathology_4_classes 
Attributes_from_overview <- cbind(Attributes_from_overview,Histpopathology_simplified  )
Attributes2 <- merge(Attributes_from_overview, gene_interest , by='Sample_ID' , all = TRUE )
Attributes2 <- merge(Attributes2, ML_prediction_df , by='Sample_ID' , all = TRUE )
Attributes2 <- merge(Attributes2, Mutation_df , by='Sample_ID' , all = TRUE )

#######################
# WRITE TABLE         #
#######################

# Fig 1A 
# _________
Fig1_LNEN_sample_id = data.frame("Sample_ID"=Coords_MOFA_fig1$Sample_ID)
Attributes_fig1A = merge(Attributes2 , Fig1_LNEN_sample_id , by="Sample_ID")
write.table(Coords_MOFA_fig1, file='Coords_MOFA_fig1.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F)
write.table(Attributes_fig1A, file='Attributes_fig1A.tsv', quote=FALSE, sep='\t', row.names = F)

# Fig 4A :
# --------

Coords_MOFA_fig4A <- data.frame("Sample_ID" = Sample_overview$Sample_ID , "Axis1" = Sample_overview$LF1.LNET *-1 , "Axis2" = Sample_overview$LF2.LNET )
Coords_MOFA_fig4A <- Coords_MOFA_fig4A[complete.cases(Coords_MOFA_fig4A),]
Sample_id_fig4A = data.frame("Sample_ID"=Coords_MOFA_fig4A$Sample_ID)
Attributes_fig4A = merge(Attributes2 , Sample_id_fig4A  , by="Sample_ID")
write.table(Coords_MOFA_fig4A, file='Coords_MOFA_fig4A.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig4A, file='Attributes_fig4A.tsv', quote=FALSE, sep='\t', row.names = F)

# Fig S6 :
# --------

# FIG 6A
# -------


PCA_RNA_seq  <- read.xlsx("../SupplementaryTables_R1_20190318.xlsx", sheet = 2, startRow = 11, colNames = TRUE,
                               rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                               skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
                               namedRegion = NULL, na.strings = "NA")

Coords_PCA_S6A <- data.frame("Sample_ID" = PCA_RNA_seq$Sample_ID , "Axis1" = PCA_RNA_seq$PC1.LNEN_SCLC , "Axis2" = PCA_RNA_seq$PC2.LNEN_SCLC)
Coords_PCA_S6A<- Coords_PCA_S6A[complete.cases(Coords_PCA_S6A),]
Sample_id_fig6A = data.frame("Sample_ID"=Coords_PCA_S6A$Sample_ID)
Attributes_fig6A = merge(Attributes2 , Sample_id_fig6A  , by="Sample_ID")
write.table(Coords_PCA_S6A,  file='Coords_PCA_S6A.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig6A, file='Attributes_fig6A.tsv', quote=FALSE, sep='\t', row.names = F)

# FIG 6B
# -------
Coords_PCA_S6B <- data.frame("Sample_ID" = PCA_RNA_seq$Sample_ID , "Axis1" = PCA_RNA_seq$PC1.LNEN , "Axis2" = PCA_RNA_seq$PC2.LNEN)
Coords_PCA_S6B<- Coords_PCA_S6B[complete.cases(Coords_PCA_S6B),]
Sample_id_fig6B = data.frame("Sample_ID"=Coords_PCA_S6B$Sample_ID)
Attributes_fig6B = merge(Attributes2 , Sample_id_fig6B  , by="Sample_ID")
write.table(Coords_PCA_S6B,  file='Coords_PCA_S6B.tsv', quote=FALSE, sep='\t', row.names = F,  col.names = F)
write.table(Attributes_fig6B, file='Attributes_fig6B.tsv', quote=FALSE, sep='\t', row.names = F)

# FIG 6C
# -------

Coords_PCA_S6C <- data.frame("Sample_ID" = PCA_RNA_seq$Sample_ID , "Axis1" = PCA_RNA_seq$PC1.LNET_SCLC , "Axis2" = PCA_RNA_seq$PC2.LNET_SCLC)
Coords_PCA_S6C<- Coords_PCA_S6C[complete.cases(Coords_PCA_S6C),]
Sample_id_fig6C = data.frame("Sample_ID"=Coords_PCA_S6C$Sample_ID)
Attributes_fig6C = merge(Attributes2 , Sample_id_fig6C  , by="Sample_ID")
write.table(Coords_PCA_S6C,  file='Coords_PCA_S6C.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig6C, file='Attributes_fig6C.tsv', quote=FALSE, sep='\t', row.names = F)



# FIG 6D
# -------

Coords_PCA_S6D <- data.frame("Sample_ID" = PCA_RNA_seq$Sample_ID , "Axis1" = PCA_RNA_seq$PC1.LNET , "Axis2" = PCA_RNA_seq$PC2.LNET)
Coords_PCA_S6D<- Coords_PCA_S6D[complete.cases(Coords_PCA_S6D),]
Sample_id_fig6D = data.frame("Sample_ID"=Coords_PCA_S6D$Sample_ID)
Attributes_fig6D = merge(Attributes2 , Sample_id_fig6D  , by="Sample_ID")
write.table(Coords_PCA_S6D,  file='Coords_PCA_S6D.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig6D, file='Attributes_fig6D.tsv', quote=FALSE, sep='\t', row.names = F)

# Fig S7 :
# --------
PCA_methylation  <- read.xlsx("../SupplementaryTables_R1_20190318.xlsx", sheet = 3, startRow = 9, colNames = TRUE,
                          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                          skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
                          namedRegion = NULL, na.strings = "NA")

# FIG 7A
# -------

Coords_PCA_S7A <- data.frame("Sample_ID" = PCA_methylation$Sample_ID , "Axis1" = PCA_methylation$PC1.LNEN, "Axis2" = PCA_methylation$PC2.LNEN)
Coords_PCA_S7A<- Coords_PCA_S7A[complete.cases(Coords_PCA_S7A),]
Sample_id_fig7A = data.frame("Sample_ID"=Coords_PCA_S7A$Sample_ID)
Attributes_fig7A = merge(Attributes2 , Sample_id_fig7A  , by="Sample_ID")
write.table(Coords_PCA_S7A,  file='Coords_PCA_S7A.tsv', quote=FALSE, sep='\t', row.names = F,  col.names = F)
write.table(Attributes_fig7A, file='Attributes_fig7A.tsv', quote=FALSE, sep='\t', row.names = F)


# FIG 7B
# -------

Coords_PCA_S7B <- data.frame("Sample_ID" = PCA_methylation$Sample_ID , "Axis1" = PCA_methylation$PC1.LNET, "Axis2" = PCA_methylation$PC2.LNET)
Coords_PCA_S7B<- Coords_PCA_S7B[complete.cases(Coords_PCA_S7B),]
Sample_id_fig7B = data.frame("Sample_ID"=Coords_PCA_S7B$Sample_ID)
Attributes_fig7B = merge(Attributes2 , Sample_id_fig7B  , by="Sample_ID")
write.table(Coords_PCA_S7B,  file='Coords_PCA_S7B.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig7B, file='Attributes_fig7B.tsv', quote=FALSE, sep='\t', row.names = F)

# Fig S13 :
# --------

# Fig S13A:
# --------
Coords_MOFA_S13A<- data.frame("Sample_ID" = Sample_overview$Sample_ID , "Axis1" = Sample_overview$LF1.LNEN_SCLC, "Axis2" = Sample_overview$LF2.LNEN_SCLC)
Coords_MOFA_S13A <- Coords_MOFA_S13A[complete.cases(Coords_MOFA_S13A),]
Sample_id_fig13A = data.frame("Sample_ID"=Coords_MOFA_S13A$Sample_ID)
Attributes_fig13A = merge(Attributes2 , Sample_id_fig13A  , by="Sample_ID")
write.table(Coords_MOFA_S13A,  file='Coords_MOFA_S13A.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig13A, file='Attributes_fig13A.tsv', quote=FALSE, sep='\t', row.names = F)


# Fig S13C:
# --------
Coords_MOFA_S13B<- data.frame("Sample_ID" = Sample_overview$Sample_ID , "Axis1" = Sample_overview$LF1.LNET_SCLC, "Axis2" = Sample_overview$LF2.LNET_SCLC)
Coords_MOFA_S13B <- Coords_MOFA_S13B[complete.cases(Coords_MOFA_S13B),]
Sample_id_fig13B = data.frame("Sample_ID"=Coords_MOFA_S13B$Sample_ID)
Attributes_fig13B = merge(Attributes2 , Sample_id_fig13B  , by="Sample_ID")
write.table(Coords_MOFA_S13B,  file='Coords_MOFA_S13B.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(Attributes_fig13B, file='Attributes_fig13B.tsv', quote=FALSE, sep='\t', row.names = F)
