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
plot(Sample_overview$LF1.LNEN, Sample_overview$LF2.LNEN , col = as.factor(Sample_overview$Histopathology ))

###########################
# Attributes               #
###########################

# Clinical attributes
# ====================
Attributes_from_overview <- data.frame("Histopathology" = Sample_overview$Histopathology , "Stage_UICC" = Sample_overview$Stage_UICC , "Age"= Sample_overview$Age , "Age_class" = Sample_overview$Age_class , 
                                      "Sex" = Sample_overview$Sex , "Smoking_status" = Sample_overview$Smoking_status , "Professional_Asbestos_exposure" = Sample_overview$Professional_exposure , "Survival_months" = Sample_overview$Survival_months,
                                      "cluster_LNEN" = Sample_overview$cluster_LNEN , "Neutrophil.to.Lymphocyte_ratio" = Sample_overview$Neutrophil.to.Lymphocyte_ratio )


Attributes_from_overview <- cbind(Attributes_from_overview , Sample_overview[ , 42:52])

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
Sample_id <- data.frame ( "Sample_ID" =Sample_overview$Sample_ID ) 
t_Data_vst_all <-as.data.frame(t( Data_vst_all ))
t_Data_vst_all <- setDT(t_Data_vst_all , keep.rownames = TRUE)[]
colnames(t_Data_vst_all)[1] <- "Sample_ID"
t_Data_vst_all <- as.data.frame(t_Data_vst_all)
Data_vst_all_with_sample <- merge(t_Data_vst_all , Sample_id , by='Sample_ID')

# Name of genes of interest :
gene_interest_names <- c("LAG3" , "IGSF11" , "VISTA_VSIR" , "PDCD1LG2" , "LGALS9" , "CD274" ,
                   "PDCD1" , "HAVCR2" , "CD86" , "CD80" , "CTLA4", names(HLA_D) )
gene_interest_embl <- c(Embl_lag3 , Embl_IGSF11,Embl_VISTA, Embl_PDCD1LG2, Embl_LGALS9, Embl_CD274, Embl_PDCD1, 
                        Embl_HAVCR2, Embl_CD86, Embl_CD80, Embl_CTLA4 )
for (i in 1:12){
  gene_interest_embl <- append(as.character(HLA_D[[i]]) ,gene_interest_embl ) 
}

# Data frame :
gene_interest_fig2E <- data.frame("Sample_ID" = Sample_overview$Sample_ID)
for (i in 1:23){
  n_col = which(colnames(t_Data_vst_all) == as.name(gene_interest_embl[i]))
  print(as.numeric(n_col))
  gene_name <- as.character(gene_interest_names[i])
  print(gene_name)
  print(t_Data_vst_all$colname)
  gene_interest_fig2E[gene_name] <- t_Data_vst_all[,n_col ]
}
