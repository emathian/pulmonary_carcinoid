######################## Pulmmonary carcinoid #################

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

#######################
# IMPORTATION OF DATA #
#######################
Sample_overview <- read.xlsx("SupplementaryTables_R1_20190318.xlsx", sheet = 1, startRow = 41, colNames = TRUE,
          rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
          skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
          namedRegion = NULL, na.strings = "NA")

Somatic_mutation  <- read.xlsx("SupplementaryTables_R1_20190318.xlsx", sheet = 4, startRow = 51, colNames = TRUE,
                              rowNames = FALSE, detectDates = FALSE, skipEmptyRows = TRUE,
                              skipEmptyCols = TRUE, rows = NULL, cols = NULL, check.names = TRUE,
                              namedRegion = NULL, na.strings = "NA")


Data_vst_50 <- read.table("data/VST_nosex_50pc_TCACLCNECSCLC.txt",  sep = " ", dec="." , header = TRUE,   quote="")
Data_vst_all <- read.table("data/VST_nosex_TCACLCNECSCLC.txt",  sep = " ", dec="." , header = TRUE,   quote="")

Ref_gene <- read.table("data/VST_nosex_50pc_TCACLCNECSCLC_annot.txt",  sep = " ", dec="." , header =FALSE,   quote="")
Ref_gene_all <- read.table("data/VST_nosex_TCACLCNECSCLC_annot.txt",  sep = " ", dec="." , header =FALSE,   quote="")


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


# Coords
# ======

Coords_MOFA_fig1 <- data.frame("Sample_ID" = Sample_overview$Sample_ID , "Axis1" = Sample_overview$LF1.LNEN , "Axis2" = Sample_overview$LF2.LNEN)
Coords_MOFA_fig1 <- Coords_MOFA_fig1[complete.cases(Sample_id_rna_seq),]
plot(Sample_overview$LF1.LNEN, Sample_overview$LF2.LNEN , col = as.factor(Sample_overview$Histopathology ))

# MOFA Sample ID
# ==============

mofa_lnen_sample_id = data.frame("Sample_ID"=Coords_MOFA_fig1$Sample_ID) # Complete case of mofa coordinates


###########################
# Attributes               #
###########################

# Clinical attributes
# ====================
Attributes_from_overview <- data.frame("Histopathology" = Sample_overview$Histopathology , "Stage_UICC" = Sample_overview$Stage_UICC , "Age"= Sample_overview$Age , "Age_class" = Sample_overview$Age_class , 
                                      "Sex" = Sample_overview$Sex , "Smoking_status" = Sample_overview$Smoking_status , "Professional_Asbestos_exposure" = Sample_overview$Professional_exposure , "Survival_months" = Sample_overview$Survival_months,
                                      "cluster_LNEN" = Sample_overview$cluster_LNEN , "Neutrophil.to.Lymphocyte_ratio" = Sample_overview$Neutrophil.to.Lymphocyte_ratio )
Attributes_from_overview <- cbind(Attributes_from_overview , Sample_overview[ , 42:52])
Attributes_from_overview <- merge(Attributes_from_overview , mofa_lnen_sample_id, by="Sample_ID")
  
# Genes of interest
# =================

# Immune checkpoint and ligand receptor Fig2E
# ____________________________________________

Embl_lag3 = as.character(Ref_gene$V1[Ref_gene$V7 == "LAG3"])
Embl_IGSF11 = as.character(Ref_gene$V1[Ref_gene$V7 == "IGSF11" ])
# VISTA VSIR : https://www.uniprot.org/uniprot/Q9H7M9
Embl_VISTA = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "C10orf54" ])
Embl_PDCD1LG2 = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "PDCD1LG2" ]) # Becarful not in 50
Embl_LGALS9  = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "LGALS9" ]) # Becarful not in 50
Embl_CD274 =  as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "CD274"]) # Becarful not in 50
Embl_PDCD1 = as.character(Ref_gene$V1[Ref_gene$V7 == "PDCD1" ])
Embl_HAVCR2 = as.character(Ref_gene$V1[Ref_gene$V7 == "HAVCR2" ])
Embl_CD86 = as.character(Ref_gene_all$V1[Ref_gene_all$V7 == "CD86"]) # Becarful not in 50
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
Data_vst_all_with_sample <- merge(t_Data_vst_all , mofa_lnen_sample_id , by='Sample_ID')
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
  print(as.numeric(n_col))
  gene_name <- as.character(gene_interest_names[i])
  print(gene_name)
  gene_interest_fig2E[gene_name] <- Data_vst_all_with_sample[,n_col ]
}

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


