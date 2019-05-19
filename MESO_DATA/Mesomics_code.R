## SCRIPT CONSTRUIT A PARTIR DU FICHIER MESO DATA

#######################
# LIBRAIRIES          #
#######################

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(data.table)
library(dplyr)
library(ade4)
library(umap)

#########################
# IMPORT DATA           #
########################
metadata <- read.table("metaBuenoTCGAred_FGS.txt", sep = "\t", dec="." , header = TRUE,   quote="")
attribute <- load("Rdata8nov2018.RData")
suptable1<- read.table("SupTables_R1_20190418/TableS1_Characteristics_of_gene-Tableau 1.csv", sep = ";", dec="," , header = TRUE,   quote="")
head(suptable1,4)

######################
#   COORDINATES     #
#####################

Coordinates <- data.frame( sample = suptable1$Samples , x = suptable1$Dimension.1 , y = suptable1$Dimension.2)
Coordinates2 <- Coordinates 
Coordinates2$y <- Coordinates2$y * -1

######################
#   CORRECT IDs     #
####################

Diff  = setdiff(groupsBuenoTCGAred$Sample , suptable1$Samples)


groupsBuenoTCGAred_2= groupsBuenoTCGAred
for (i in  1:284){
  if (groupsBuenoTCGAred$Sample[i] %in% Diff){
    groupsBuenoTCGAred_2$Sample[i] = paste("M", groupsBuenoTCGAred$Sample[i], sep="" )
  }
}
setdiff(groupsBuenoTCGAred_2$Sample, suptable1$Samples)

########################
#     ATTRIBUTES      #
#######################


TILS <- data.frame( sample = suptable1$Samples , B = suptable1$B , Macrophages.M1 = suptable1$Macrophages.M1 , 
                    Macrophages.M2 = suptable1$Macrophages.M2 , Monocytes = suptable1$Monocytes, 
                    Neutrophils = suptable1$Neutrophils , NK.cells = suptable1$Macrophages.M2 ,
                    T.cells.CD4 = suptable1$T.cells.CD4  , T.cells.CD8 = suptable1$T.cells.CD8 , 
                    Tregs = suptable1$Tregs , Dendritic.cells = suptable1$Dendritic.cells )


ClinicalAttributes <- data.frame( sample = groupsBuenoTCGAred_2$Sample, Asbestos =  groupsBuenoTCGAred_2$Asbestos , 
                                  Sex= groupsBuenoTCGAred_2$Sex , Type = groupsBuenoTCGAred_2$Type , 
                                  Age = groupsBuenoTCGAred_2$Age , Survival = groupsBuenoTCGAred_2$Survival , cohort=metadata$Study ) 

# Finaly not include 
#Subtype_Type_revised =  data.frame( sample = as.character(metadata$Sample) , Subtype_FGS = metadata$Subtype.FGS , Type_FGS = metadata$Type.FGS , Necrosis = metadata$Necrosis, grade = metadata$Grade, stringsAsFactors=FALSE)
#str(Subtype_Type_revised)
#Diff= setdiff(Subtype_Type_revised$sample, suptable1$Samples); Diff
#for (i in  1:284){
#  if (Subtype_Type_revised$sample[i] %in% Diff){
#    Subtype_Type_revised$sample[i] = paste("M",as.character( Subtype_Type_revised$sample[i]), sep="" )
#  }
# }



#####################
# GENE OF INTEREST  #
##################### 


VISTA = c(dinomesomics["ENSG00000107738.19",])
CD8A = c(dinomesomics["ENSG00000153563.15",])
VEGFR1 = c(dinomesomics["ENSG00000102755.10",])
VEGFR2 = c(dinomesomics["ENSG00000128052.8",])
VEGFR3 = c(dinomesomics["ENSG00000037280.15",])
VEGFC = c(dinomesomics["ENSG00000150630.3",])
PDGFRB = c(dinomesomics["ENSG00000113721.13",])
PD1  = c(dinomesomics["ENSG00000188389.10",])
PDL1 = c(dinomesomics["ENSG00000120217.13",])
CTLA4 = c(dinomesomics["ENSG00000163599.14",])
TIM3 = c(dinomesomics["ENSG00000135077.8",])
LAG3 = c(dinomesomics["ENSG00000089692.8",])

gene_interet <- data.frame(VISTA,CD8A,VEGFR1, VEGFR2,VEGFR3, VEGFC , PDGFRB,PD1,PDL1, CTLA4 , TIM3 , LAG3 )
# Rmq quant au nom des ID
ID = sampleFilesred
New_ID = c()
for (i in 1:length(ID)){
  New_ID[i] =  sub("_count.txt","", ID[i])  
}

# Création d'un tableau de référence entre les IDs et le noms des éch

metadataID = data.frame(ID = New_ID, sample = suptable1$Sample)
gene_interet = setDT(gene_interet, keep.rownames = TRUE)[]
colnames(metadataID)[1] <- "ID"
colnames(gene_interet)[1] <- "ID"

gene_interet_samples_names <- merge(metadataID, gene_interet, by="ID")
gene_interet2 <- gene_interet_samples_names[,-1]

# PROFILES
#---------

max(suptable1$Dimension.1)
min(suptable1$Dimension.1)
max(suptable1$Dimension.2)
min(suptable1$Dimension.2)


profiles <- c()
for (i in 1:dim(suptable1)[1]){
  if (suptable1$Dimension.1[i] <= min(suptable1$Dimension.1) + (max(suptable1$Dimension.1)   + -1*min(suptable1$Dimension.1)) *1/3  & 
      suptable1$Dimension.2[i] >= min(suptable1$Dimension.2) + (max(suptable1$Dimension.2)   + -1*min(suptable1$Dimension.2)) *2/3 ){
    profiles[i] <- "Hot/IC+/Angio+"
  }
  else if (suptable1$Dimension.1[i] >=  min(suptable1$Dimension.1) + (max(suptable1$Dimension.1)   + -1*min(suptable1$Dimension.1)) *2/3
           & suptable1$Dimension.2[i] >= min(suptable1$Dimension.2) + (max(suptable1$Dimension.2)   + -1*min(suptable1$Dimension.2)) *2/3  ){
    profiles[i] <- "VEGFR2+/VISTA+"
  }
  else if (suptable1$Dimension.1[i] <=  min(suptable1$Dimension.1) + (max(suptable1$Dimension.1)   + -1*min(suptable1$Dimension.1)) *1/3   
           & suptable1$Dimension.2[i] <= min(suptable1$Dimension.2) + (max(suptable1$Dimension.2)   + -1*min(suptable1$Dimension.2)) *1/3  ){
    profiles[i] <- "Cold/Angio+"
  }
  else{
    profiles[i]<- NA
  }
}


plot(suptable1$Dimension.1 , suptable1$Dimension.2 , col= as.factor(profiles))

profile_df <- data.frame("sample" = suptable1$Sample , "Profiles" =profiles )
#################################
#   MERGE ATTRIBUTES            #
#################################



Attributes1 <- merge(TILS, ClinicalAttributes, by="sample")

Attributes2 <- merge(Attributes1, gene_interet2, by="sample")
Attributes2 <- merge(Attributes2 , profile_df , by="sample")
Attributes3 <- Attributes2 # FOR SIMPLIFY THE FOLLOWING SCRIPT

# Convert factor to char
Attributes3$Type = as.character(Attributes3$Type)
table(Attributes3$Type)


# Survie en mois
Attributes3$Survival<- Attributes3$Survival*12


#############################
#   ADD IHC DATA            #
#############################
IHC_333.2 <- read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_333_2.csv", sep = ";", dec="." , header = TRUE, na.strings=c("","NA"))
IHC_334.2 <- read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_334_2.csv", sep = ";", dec="." , header = TRUE, na.strings=c("","NA"))
IHC_351.1 <- read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_351_1.csv", sep = ";", dec="." , header = TRUE, na.strings=c("","NA"))
IHC_352.1 <-  read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_352_1.csv", sep = ";", dec="." , header = TRUE, na.strings=c("","NA"))
IHC_353.1 <-  read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_353_1.csv", sep = ";", dec="." , header = TRUE, na.strings=c("","NA"))
IHC_354.1 <-   read.csv("../IHC_Lecture_LM_20181121/IHC_Lecture_LM_20181116_354_1.csv", sep = ";", dec="." , header = TRUE,na.strings=c("","NA"))

change_sample_name_IHC_file <- function(data){
  TGCA_ID <- c()
  for (i in 1:dim(data)[1]){
    if (is.na(data$IARC_ID[i])==F){
      if (sum(is.na(data[i,4:8])) != 5){
        withM <- test_ID <- sub("GNT", "M", data$IARC_ID[i])[1]    
        g_expr = gregexpr(pattern ='_',withM)[[1]][1]
        if  (g_expr !=  -1){
          TGCA_ID[i] <- paste(substr(withM, 1, g_expr-1), "PT", sep="")[1]
        }
        else{
          TGCA_ID[i] <- paste(withM, "PT", sep="")[1]
        }
      }
      else{
        TGCA_ID[i] <- NA
      }
    }
    else{
      TGCA_ID[i] <- NA
    }
  }
  return(TGCA_ID)
}

IHC_333.2["Sample"] =  change_sample_name_IHC_file(IHC_333.2)
IHC_334.2["Sample"]  = change_sample_name_IHC_file(IHC_334.2)
IHC_351.1["Sample"]  = change_sample_name_IHC_file(IHC_351.1)
colnames(IHC_352.1)[2] <- "IARC_ID"
IHC_352.1["Sample"]  = change_sample_name_IHC_file(IHC_352.1)
colnames(IHC_353.1)[2] <- "IARC_ID"
IHC_353.1["Sample"]  = change_sample_name_IHC_file(IHC_353.1)
colnames(IHC_354.1)[2] <- "IARC_ID"
IHC_354.1["Sample"]  = change_sample_name_IHC_file(IHC_354.1)


# Pour CD8 et VEGFR3M passage en pourcentage 
# ------------------------------------------

CD8_VEGFR3_Pourcent <- function(data){
  CD8_p <- c()
  VEGFR3M_p <- c()
  for (i in 1:dim(data)[1]){
    if (is.na(data$CD8[i]) == F){
      if (data$CD8[i]==0){
        CD8_p[i]<- 0
      }
      else if(data$CD8[i]==1){
        CD8_p[i]<- 25
      }
      else if(data$CD8[i]==2){
        CD8_p[i]<- 50
      }
      else if(data$CD8[i]==3){
        CD8_p[i]<- 75
      }
      else if(data$CD8[i]==4){
        CD8_p[i]<- 100
      }
      else {
        CD8_p[i]<- NA
      }
    }
    else{
      CD8_p[i]<- NA
    }
    if (is.na(data$VEGFR3M[i])==F){
      if (data$VEGFR3M[i]==0){
        VEGFR3M_p[i]<- 0
      }
      else if(data$VEGFR3M[i]==1){
        VEGFR3M_p[i]<- 25
      }
      else if(data$VEGFR3M[i]==2){
        VEGFR3M_p[i]<- 50
      }
      else if(data$VEGFR3M[i]==3){
        VEGFR3M_p[i]<- 75
      }
      else if(data$VEGFR3M[i]==4){
        VEGFR3M_p[i]<- 100
      }
      else {
        VEGFR3M_p[i]<- NA
      }
    }
    else{
      VEGFR3M_p[i]<- NA
    }
  }
  CD8_VEGFR3_Pourcent_df <- data.frame(CD8_p, VEGFR3M_p) 
  return (CD8_VEGFR3_Pourcent_df)
}


IHC_333.2["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_333.2)[,1]
IHC_333.2["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_333.2)[,2]

IHC_334.2["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_334.2)[,1]
IHC_334.2["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_334.2)[,2]   

IHC_351.1["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_351.1)[,1]
IHC_351.1["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_351.1)[,2]    

IHC_352.1["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_352.1)[,1]
IHC_352.1["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_352.1)[,2] 


IHC_353.1["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_353.1)[,1]
IHC_353.1["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_353.1)[,2] 

IHC_354.1["CD8_p"]<-CD8_VEGFR3_Pourcent(IHC_354.1)[,1]
IHC_354.1["VEGFR3M_p"]<-CD8_VEGFR3_Pourcent(IHC_354.1)[,2] 


# Création d'une table finale 
#  __________________________


IHC_all <- data.frame()

IHC_333.2_2 <- data.frame( Sample = IHC_333.2$Sample , VISTA = IHC_333.2$VISTA, PDL1 = IHC_333.2$PDL1 , VEGFR3 = IHC_333.2$VEGFR3M_p , CD8 = IHC_333.2$CD8_p , VEGFR2 = IHC_333.2$VEGFR2)

IHC_all <- rbind(IHC_333.2_2 ,IHC_all)


IHC_334.2_2<- data.frame( Sample = IHC_334.2$Sample , VISTA = IHC_334.2$VISTA, PDL1 = IHC_334.2$PDL1 , VEGFR3 = IHC_334.2$VEGFR3M_p , CD8 = IHC_334.2$CD8_p , VEGFR2 = IHC_334.2$VEGFR2)

IHC_all <- rbind(IHC_334.2_2 ,IHC_all)

IHC_351.1_2<- data.frame( Sample = IHC_351.1$Sample , VISTA = IHC_351.1$VISTA, PDL1 = IHC_351.1$PDL1 , VEGFR3 = IHC_351.1$VEGFR3M_p , CD8 = IHC_351.1$CD8_p , VEGFR2 = IHC_351.1$VEGFR2)

IHC_all <- rbind(IHC_351.1_2 ,IHC_all)

IHC_352.1_2<- data.frame( Sample = IHC_352.1$Sample , VISTA = IHC_352.1$VISTA, PDL1 = IHC_352.1$PDL1 , VEGFR3 = IHC_352.1$VEGFR3M_p , CD8 = IHC_352.1$CD8_p , VEGFR2 = IHC_352.1$VEGFR2)
IHC_all <- rbind(IHC_352.1_2 ,IHC_all)


IHC_354.1_2<- data.frame( Sample = IHC_354.1$Sample , VISTA = IHC_354.1$VISTA, PDL1 = IHC_354.1$PDL1 , VEGFR3 = IHC_354.1$VEGFR3M_p , CD8 = IHC_354.1$CD8_p , VEGFR2 = IHC_354.1$VEGFR2)
IHC_all <- rbind(IHC_354.1_2 ,IHC_all)

tot_sample <- c(IHC_333.2$Sample , IHC_334.2$Sample, IHC_351.1$Sample, IHC_352.1$Sample,  IHC_353.1$Sample , IHC_354.1$Sample     ) #IHC_411.1$Sample , IHC_412.1$Sample 
uniq_tot_sample= unique(tot_sample)[-1]  


# Moyenne par échantillon
# -----------------------

IHC_moy_df <- data.frame(sample =uniq_tot_sample  ,  VEGFR2.IHC = rep(0,length(uniq_tot_sample)), VEGFR3M.IHC = rep(0,length(uniq_tot_sample)), 
                         PDL1.IHC = rep(0,length(uniq_tot_sample)), CD8.IHC = rep(0,length(uniq_tot_sample)), 
                         VISTA.IHC = rep(0,length(uniq_tot_sample)),  stringsAsFactors = F)

IHC_moy_sample <- function(data, res_df){
  for (i in 1:length(uniq_tot_sample)){
    c_sample_id = IHC_moy_df$sample[i]
    c_CD8_m =  mean(as.numeric(IHC_all$CD8[which(IHC_all$Sample== c_sample_id)]), na.rm=TRUE)
    res_df[i,5] <-c_CD8_m # Ok
    c_VISTA =  mean(as.numeric(IHC_all$VISTA[which(IHC_all$Sample== c_sample_id)]), na.rm=TRUE)
    res_df[i,6] <-c_VISTA 
    c_PDL1 =  mean(as.numeric(IHC_all$PDL1[which(IHC_all$Sample== c_sample_id)]), na.rm=TRUE)
    res_df[i,4] <-c_PDL1 #OK
    c_VEGFR3 =  mean(as.numeric(IHC_all$VEGFR3[which(IHC_all$Sample== c_sample_id)]), na.rm=TRUE)
    res_df[i,3] <-c_VEGFR3
    c_VEGFR2 =  mean(as.numeric(IHC_all$VEGFR2[which(IHC_all$Sample== c_sample_id)]), na.rm=TRUE)
    res_df[i,2] <- c_VEGFR2 
  }
  return(res_df)
  
}

IHC_moy_df = IHC_moy_sample(IHC_all , IHC_moy_df)
# Supression des données pour lesquels les données génomiques sont manquantes
v1 = which(IHC_moy_df$sample=="M654PT")
IHC_moy_df = IHC_moy_df [-v1,]
v2 = which(IHC_moy_df$sample=="M78PT")
IHC_moy_df = IHC_moy_df [-v2,]
v3 = which(IHC_moy_df$sample=="M601PT")
IHC_moy_df = IHC_moy_df [-v3,]



# Attributes with IHC 
# --------------------

Attributes4 <- merge(Attributes3 , IHC_moy_df , by= "sample" , all= TRUE) 
colnames(Attributes4)[which(colnames(Attributes4)=="B")] <- "B.Cells"


############################################
# Feature_DATA with the largest variance   #
############################################

data_lv  = merge(data , dataNOS , by=0) # Merge by row names
f_colnames = data_lv$Row.names
data_lv = as.data.frame(t(data_lv), col.names =FALSE )
colnames(data_lv) <- f_colnames 
data_lv = setDT(data_lv , keep.rownames = TRUE)[]
data_lv = data_lv[-1,]
colnames(data_lv)[1] <- "ID"
data_lv_t_sample = merge(metadataID,data_lv, by ='ID')
data_lv_t_sample = data_lv_t_sample[,-1]
data_lv_t_sample = t(data_lv_t_sample)


#######################################@
#             ACP Fig 1                #
########################################
df_acp_fig1 <-data_lv
df_acp_fig1.2 = as.data.frame( df_acp_fig1[,2:7146] )
rownames(df_acp_fig1.2) <- df_acp_fig1$ID
for (i in 1:7145){
  df_acp_fig1.2[,i] <- as.numeric(as.character(df_acp_fig1.2[,i]) )
}
str(df_acp_fig1.2)
acp_fig1 <- dudi.pca( df_acp_fig1.2 ,center = T , scale = F , scannf = F , nf=2) #
s.label(acp_fig1$li, xax = 1, yax = 2)
acp_fig1_li_df =  as.data.frame(acp_fig1$li)
acp_fig1_li_df  = setDT(acp_fig1_li_df , keep.rownames = TRUE)[]
colnames(acp_fig1_li_df )[1]<-'ID'
acp_fig1_li_df = merge(acp_fig1_li_df , metadataID , by ='ID')
acp_fig1_li_df = acp_fig1_li_df[,-1]

acp_fig1_li_df = merge(Attributes4 , acp_fig1_li_df , by="sample")


gcol = c("blue", "red", "orange", "darkgreen")
Type_acp_fig1  = as.factor(Attributes4$Type)
#s.class(acp_fig1$li,  Type_acp_fig1  , col = gcol, xax=1, yax=2)
colnames(suptable1)[1]<-"Sample"
suptable1.2 <- merge(suptable1, metadata, by='Sample')
plot(acp_fig1_li_df$Axis1 , acp_fig1_li_df$Axis2 , col=as.factor(acp_fig1_li_df$Type) )  # ;)   Victoire
plot(suptable1.2$Dimension.1 , suptable1.2$Dimension.2 , col=suptable1.2$Type ) #OK





######################
# Reproduction fig 3 #
######################

# Importation de la cohorte de replication 
# =========================================

Replication_data <-  read.csv("IHC_discovery_+_replication_cohorts/dataIHC.csv", sep = ",", dec="." , header = TRUE, na.strings=c("","NA"))

# Add IHC data
#_-----------

Vista_r_mb = Replication_data$VISTA.membranous
CD8_r = Replication_data$CD8
Replication_data$VEGFR2[26] = "1M"
Vegfr2_r_cyto = Replication_data$VEGFR2
Vegfr2_r_cyto.f = c()
for (i in 1:77){
  if (is.na(Vegfr2_r_cyto[i])==F){
    if (Vegfr2_r_cyto[i] =='0'){
      Vegfr2_r_cyto.f[i] = 0
    }
    else if (Vegfr2_r_cyto[i]=='1M'){
      Vegfr2_r_cyto.f[i] = 0.25
    }
    else if (Vegfr2_r_cyto[i]=='2M'){
      Vegfr2_r_cyto.f[i] = 0.5
    }
    else if (Vegfr2_r_cyto[i]=='3M'){
      Vegfr2_r_cyto.f[i] = 0.75
    }
    else if (Vegfr2_r_cyto[i]=='4M'){
      Vegfr2_r_cyto.f[i] = 1
    }
  }
  else{
    Vegfr2_r_cyto.f[i] = NA
  }
}


Vegfr2_r_mb = Replication_data$VEGFR2.membranous

PDL1_r = Replication_data$PDL1.Tumor
Vegfr3_r = Replication_data$VEGFR3.membranous
Vegfr3_r.f = c()
for (i in 1:77){
  if (is.na(Vegfr3_r[i])==F){
    if (Vegfr3_r[i] ==0){
      Vegfr3_r.f[i] = 0
    }
    else if (Vegfr3_r[i]==1){
      Vegfr3_r.f[i] = 0.25
    }
    else if (Vegfr3_r[i]==2){
      Vegfr3_r.f[i] = 0.5
    }
    else if (Vegfr3_r[i]==3){
      Vegfr3_r.f[i] = 0.75
    }
    else if (Vegfr3_r[i]==4){
      Vegfr3_r.f[i] = 1
    }
  }
  else{
    Vegfr3_r.f[i] = NA
  }
}

Vegfr3_r_cyto =Replication_data$VEGFR3.cytoplasmic
Vegfr3_r_cyto.f =c()
for (i in 1:77){
  if (is.na(Vegfr3_r_cyto[i])==F){
    if (Vegfr3_r_cyto[i] ==0){
      Vegfr3_r_cyto.f[i] = 0
    }
    else if (Vegfr3_r_cyto[i]==1){
      Vegfr3_r_cyto.f[i] = 0.25
    }
    else if (Vegfr3_r_cyto[i]==2){
      Vegfr3_r_cyto.f[i] = 0.5
    }
    else if (Vegfr3_r_cyto[i]==3){
      Vegfr3_r_cyto.f[i] = 0.75
    }
    else if (Vegfr3_r_cyto[i]==4){
      Vegfr3_r_cyto.f[i] = 1
    }
  }
  else{
    Vegfr3_r_cyto.f[i] = NA
  }
}

PDL1.TILS_r = Replication_data$PDL1.TILS
PDL1.TILS_r.f = c()
for (i in 1:77){
  if (is.na(PDL1.TILS_r[i])==F){
    if (PDL1.TILS_r[i] ==0){
      PDL1.TILS_r.f[i] = 0
    }
    else if (PDL1.TILS_r[i]==1){
      PDL1.TILS_r.f[i] = 0.25
    }
    else if (PDL1.TILS_r[i]==2){
      PDL1.TILS_r.f[i] = 0.5
    }
    else if ( PDL1.TILS_r[i]==3){
      PDL1.TILS_r.f[i] = 0.75
    }
    
  }
  else{
    PDL1.TILS_r.f[i] = NA
  }
}


# ACP avec les marqueurs IHC
# --------------------------

replication_acp = data.frame(Vista_r_mb, Vegfr2_r_mb , PDL1_r, PDL1.TILS_r.f,  CD8_r,Vegfr3_r.f  ) #     , ,Vegfr2_r_cyto.f Vegfr3_r_cyto.f
rownames(replication_acp)=Replication_data$X
replication_acp = replication_acp[complete.cases(replication_acp),]
acp_rep <- dudi.pca(replication_acp , center = T  , scale = F , scannf = F , nf =2)
acp_rep_li = as.data.frame(acp_rep$li)
acp_rep_li = setDT(acp_rep_li, keep.rownames = TRUE)[]
colnames(acp_rep_li)[1] <- "X"
Attributes_replication = merge(Replication_data, acp_rep_li , by="X", all=TRUE)
replication_acp = setDT(replication_acp, keep.rownames = TRUE)[]
colnames(replication_acp)[1] = "X"
Attributes_replication = merge(Attributes_replication , replication_acp , by="X", all=FALSE)
plot(Attributes_replication$Axis1*-1, Attributes_replication$Axis2 , col = as.factor(Attributes_replication$GRP)) # Victoire

Coord_fig3_b <- data.frame("X"=Attributes_replication$X ,"Axis1"= Attributes_replication$Axis1 , "Axis2" = Attributes_replication$Axis2)
Coord_fig3_b$X = as.character(Coord_fig3_b$X)

sample.2 =c()
for (i in 1:63){
  sample.2[i]= paste("M", Coord_fig3_b$X[i], sep="" )
}
Coord_fig3_b$X = sample.2
Coord_fig3_b$Axis1 = Coord_fig3_b$Axis1*-1
Coord_fig3_b$Axis2 = Coord_fig3_b$Axis2*-1




Attribute_fig3a_b <-data.frame( "sample"=Attributes_replication$X, "Sex" = Attributes_replication$SEXE ,"Age"  = Attributes_replication$Age, "Asbestos"= Attributes_replication$Expo.A.codée , "Smoking"= Attributes_replication$Tabagisme , "Type"= Attributes_replication$GRP , "Survival"= Attributes_replication$Survie, "CD8.IHC" =Attributes_replication$CD8_r  , "Vegfr2.mb.IHC"= Attributes_replication$Vegfr2_r_mb , "Vegfr3.mb.IHC"= Attributes_replication$Vegfr3_r.f , "PDL1.IHC" = Attributes_replication$PDL1_r , "PDL1.TILS" = Attributes_replication$PDL1.TILS_r.f , "VISTA.IHC" = Attributes_replication$Vista_r_mb)

Attribute_fig3a_b$Asbestos <- as.character(Attribute_fig3a_b$Asbestos)
Asbestos.2 =c()
for (i in 1:63){
  if (is.na(Attribute_fig3a_b$Asbestos[i])==F){
    if (Attribute_fig3a_b$Asbestos[i] == "non"){
      
      Asbestos.2[i] = "No"
    }
    else if (Attribute_fig3a_b$Asbestos[i] == "oui"){
      Asbestos.2[i] = "Yes"
    }
    else  if (Attribute_fig3a_b$Asbestos[i] == "poss"){
      Asbestos.2[i] = "Possible"
    }
    else  if (Attribute_fig3a_b$Asbestos[i]== "prob"){
      Asbestos.2[i] = "Probable"
    }
  }
  else{
    Asbestos.2[i] = NA
  }
}


Attribute_fig3a_b$Asbestos =  Asbestos.2
Attribute_fig3a_b$Smoking <- as.character(Attribute_fig3a_b$Smoking)
Smoking.2 = c()
for (i in 1:63){
  if (is.na(Attribute_fig3a_b$Smoking[i])==F){
    if (Attribute_fig3a_b$Smoking[i] == "Non fumeur"){
      
      Smoking.2[i] = "No.Smoker"
    }
    else if (Attribute_fig3a_b$Smoking[i] == "Fumeur"){
      Smoking.2[i] = "Smoker"
    }
    else  if (Attribute_fig3a_b$Smoking[i] == "Ex-Fumeur"){
      Smoking.2[i] = "Former.smoker"
    }
  }
  else{
    Smoking.2[i] = NA
  }
}
Attribute_fig3a_b$Smoking <- Smoking.2

Attribute_fig3a_b$Type <- as.character(Attribute_fig3a_b$Type)
Type.2 <-c()

for (i in 1:63){
  if (is.na(Attribute_fig3a_b$Type[i])==F){
    if (Attribute_fig3a_b$Type[i] == "MME.long.surv"){
      
      Type.2[i] = "Epithelioid.long.survival"
    }
    else if (Attribute_fig3a_b$Type[i] == "MME.short.surv"){
      Type.2[i] = "Epithelioid.short.survival"
    }
    else  if (Attribute_fig3a_b$Type[i] == "MMS"){
      Type.2[i] = "Sarcomatoid"
    }
  }
  else{
    Type.2[i] = NA
  }
}

Attribute_fig3a_b$Type = Type.2

Attribute_fig3a_b$sample = as.character(Attribute_fig3a_b$sample)
sample.2 =c()
for (i in 1:63){
  sample.2[i]= paste("M", Attribute_fig3a_b$sample [i], sep="" )
}

Attribute_fig3a_b$sample = sample.2


# Fig 3A Top
#___________
# Fig 3a top Panel (29/04)

# Echantillons 
Surv_table <- read.table("surv_data.txt", sep = "\t", dec="." , header = TRUE,   quote="")
# Df avec échantillons et type_+_survie
Surv_table_2 <- data_frame("sample" = Surv_table$sample , "Type_and_survival" = Surv_table$Group)
# Joint Surv_table_2 avec les 284 élmts de Attributes4
Attributes5 <- merge(Attributes4 , Surv_table_2 , by="sample")
# Ajout des ID au tableau
Attributes5 <- merge(metadataID  , Attributes5 , by="sample")
# Supression des Type_+_survival non def
To_delete = which(is.na(Attributes5$Type_and_survival)==T )
data_attributes_fi3a = Attributes5[-c(To_delete),]
# Taille du df 113 ech data_attributes_fi3a

# Liste des IDs conservés
# Création de la matrice de comptage pour les 113 échantillons
ID_selected = as.character(data_attributes_fi3a$ID)
ID_selected_df = data.frame("ID"= ID_selected)
data_lv = setDT(data_lv, keep.rownames = TRUE)[]
data_lv_fig3 = merge(data_lv ,ID_selected_df , by="ID")
rownames(data_lv_fig3) = data_lv_fig3$ID
data_lv_fig3 = data_lv_fig3[,-1]
data_lv_fig3 = as.data.frame(data_lv_fig3 )

for ( i in 1:7145){
  data_lv_fig3[,i]= as.numeric(as.character(data_lv_fig3[,i]))
}
rownames(data_lv_fig3) = ID_selected_df$ID

# ACP avec tous les 7145 gènes
#------------------------------
acp_fig3a <- dudi.pca(data_lv_fig3, center = T , scale = F , scannf = F, nf=2)
summary(acp_fig3a)
acp_fig3a_li <- acp_fig3a$li
acp_fig3a_li <- setDT(acp_fig3a_li, keep.rownames = TRUE)[]
colnames(acp_fig3a_li)[1] = "ID"

data_attributes_fi3a = merge (data_attributes_fi3a ,acp_fig3a_li , by='ID' )
plot(data_attributes_fi3a$Axis1 , data_attributes_fi3a$Axis2 , col = data_attributes_fi3a$Type_and_survival)

# Fig 3a ACP avec les 5 gènes d'intérêt
# --------------------------------------


VISTA =  dinomesomics["ENSG00000107796.12",]  #ENSG00000107738.19
ID_dinomesomics = colnames(dinomesomics)
VISTA.f = c()
Sample_for_vista= c()
data_attributes_fi3a_ID =as.list (as.character(data_attributes_fi3a$ID))
c= 0                          
for (i in 1:284){
  if (ID_dinomesomics[i] %in% data_attributes_fi3a_ID){
    c = c+1
    VISTA.f[i]= VISTA[i]
    Sample_for_vista[i]=as.character(Attributes5$sample[which(Attributes5$ID == ID_dinomesomics[i])])
  }
}

Vista_sample_df  = data.frame("sample" = Sample_for_vista , "VISTA.f"= VISTA.f)
Vista_sample_df = Vista_sample_df[complete.cases(Vista_sample_df),]


CD8A = data_lv_fig3$ENSG00000153563.15
CD8A.f = c()
for (i in 1:113){
  CD8A.f[i]=CD8A[[i]]
}
VEGFR2 = data_lv_fig3$ENSG00000128052.8
VEGFR2.f = c()
for (i in 1:113){
  VEGFR2.f[i]=VEGFR2[[i]]
}

VEGFR3 = data_lv_fig3$ENSG00000037280.15
VEGFR3.f = c()
for (i in 1:113){
  VEGFR3.f[i]=VEGFR3[[i]]
}

PDL1 = data_lv_fig3$ENSG00000120217.13
PDL1.f = c()
for (i in 1:113){
  PDL1.f[i]=PDL1[[i]]
}


dat_fig3a_5_genes = data.frame("sample" = data_attributes_fi3a$sample  , "PDL1"= PDL1.f, "VEGFR2"  = VEGFR2.f , "VEGFR3"=VEGFR3.f , "CD8A" = CD8A.f )
dat_fig3a_5_genes = merge(Vista_sample_df,dat_fig3a_5_genes , by="sample")
rownames(dat_fig3a_5_genes) <-dat_fig3a_5_genes$sample
dat_fig3a_5_genes <- dat_fig3a_5_genes[,-1]

acp_5_genes <- dudi.pca(dat_fig3a_5_genes, center = T , scale = F, scannf = F, nf=2)
summary(acp_5_genes)
acp_5_genes_li <- acp_5_genes$li
acp_5_genes_li <- setDT(acp_5_genes_li, keep.rownames = TRUE)[]
colnames(acp_5_genes_li)[1]="sample"
data_attributes_fi3a= merge(data_attributes_fi3a, acp_5_genes_li , by="sample")
to_write_data_attributes_fi3a_top = data_attributes_fi3a[ , -2]
plot(data_attributes_fi3a$Axis1.y , data_attributes_fi3a$Axis2.y , col = data_attributes_fi3a$Type_and_survival)


# Test du 30/05 Reproduction de l'ACP 
Keep_sample = as.character(Surv_table$sample[is.na(Surv_table$Group)==F ])
Keep_ID = as.character(metadataID$ID[as.character(metadataID$sample)%in%Keep_sample])
Keep_ID_Sample = data.frame("sample" = as.character(Keep_sample), "ID" = as.character(Keep_ID))
Surv_table_complete = Surv_table[complete.cases(Surv_table$Group),]
Surv_table_complete = merge(Surv_table_complete , Keep_ID_Sample , by="sample")
data_lv_with_surv_table =merge(data_lv ,Keep_ID_Sample , by='ID' )
data_lv_with_surv_table.ID = as.character(data_lv_with_surv_table$ID )
data_lv_with_surv_table.sample = as.character(data_lv_with_surv_table$sample)
# 
which(colnames(data_lv_with_surv_table)=="sample")
data_lv_with_surv_table= data_lv_with_surv_table[,-7147]
data_lv_with_surv_table= data_lv_with_surv_table[,-1]
data_lv_with_surv_table[] <- lapply(data_lv_with_surv_table[,1:7145], function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
#sapply(data_lv_with_surv_table, class)
rownames(data_lv_with_surv_table) = data_lv_with_surv_table.sample
acp_data_lv_surv <- dudi.pca(data_lv_with_surv_table, center = T , scale = F , scannf = F , nf=2)
acp_data_lv_surv_li <- acp_data_lv_surv$li
acp_data_lv_surv_li["sample"] = data_lv_with_surv_table.sample
Surv_table_complete = merge(Surv_table_complete, acp_data_lv_surv_li , by="sample")
plot(Surv_table_complete$Axis1 , Surv_table_complete$Axis2*-1 , col = Surv_table_complete$Group )


##################################
# COORDONNEES Figure3A top Panel #
##################################

#head(Keep_ID_Sample ,3)
colnames(Keep_ID_Sample)[1]<- "Sample"
#head(suptable1 ,3)
Coordinates_fig3a_top <- merge( Keep_ID_Sample , suptable1 ,by = "Sample")
Coordinates_fig3a_top$Dimension.2 =Coordinates_fig3a_top$Dimension.2*-1
Coordinates_fig3a_top <- Coordinates_fig3a_top[,-2] # Remove IDs
Coordinates_fig3a_top <- Coordinates_fig3a_top[,1:3]








##################################
#   Checking before writing      # 
##################################

head(Coordinates2,3)
head(Attributes4, 3)
head(data_lv_t_sample, 3)
head(Attribute_fig3a_b,3)
head(Coord_fig3_b,3)
head(Coordinates_fig3a_top,3)
head(to_write_data_attributes_fi3a_top,3)

dim(to_write_data_attributes_fi3a_top)
to_write_data_attributes_fi3a_top <- to_write_data_attributes_fi3a_top[,-c(37,38,39,40)]
head(to_write_data_attributes_fi3a_top,3)
############################
#   WRITE FILES           #
###########################

# Coordinates
# ____________

# Fig 1  with ACP coordinates
#____________________________


write.table(Coordinates2, file='Coordinates_PCA_f2_invY.tsv', quote=FALSE, sep='\t', row.names = F) # Coordinates2 = coordinates1 with y*-1
write.table(Attributes4, file='Attributes_PCA_f1.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F)

# Fig 1 with Features Data
#_________________________

write.table(Attributes4, file='Attributes_4_with_IHC.tsv', quote=FALSE, sep='\t', row.names = F)
write.table(data_lv_t_sample, file='feature_data_with_lv_2.tsv', quote=FALSE, sep='\t', row.names = T, col.names = F)



# Fig 3 bottom  panel
#_________________________
write.table(Attribute_fig3a_b, file='Attribute_fig3a_b.tsv', quote=FALSE, sep='\t', row.names = F)
write.table(Coord_fig3_b, file='Coordinates_PCA_f3a_B.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F)


# Fig 3 top panel
#_________________________

write.table(Coordinates_fig3a_top, file='Coordinates_fig3a_top.tsv', quote=FALSE, sep='\t', row.names = F, col.names = F)
write.table(to_write_data_attributes_fi3a_top, file='to_write_data_attributes_fi3a_top.tsv', quote=FALSE, sep='\t', row.names = F)




#################################
#             UMAP              #
#################################

umap.defaults

data_lv_sample <- t(data_lv_t_sample)  

data_lv.L <- data.frame('sample'= data_lv_sample[,1] )
Type_df <-data.frame( 'sample'= Attributes4$sample , 'Type'=Attributes4$Type )
data_lv.L <- merge(data_lv.L , Type_df , by = "sample")
data_lv.L <- data_lv.L[order(data_lv.L$sample),] 
  

data_lv.D <-as.data.frame( data_lv_sample )
data_lv.D <- data_lv.D[order(data_lv.D$sample),]
data_lv.D <- data_lv.D[,-1] 
data_lv.D <- apply(data_lv.D, 2, as.numeric)
test1.umap = umap(data_lv.D)

coord_PCA1 <- merge(Coordinates2, Type_df , by = "sample")

par(mfrow=c(1,2))
plot(test1.umap$layout[,1] ,test1.umap$layout[,2], col= c("orange", "black", "forestgreen", "brown4")[as.factor(data_lv.L$Type)], xlab = "x" , ylab ="y" , pch=20  )
plot(coord_PCA1$x , coord_PCA1$y , col= c("orange", "black", "forestgreen", "brown4")[as.factor(coord_PCA1$Type)], xlab = "x" , ylab ="y" , pch = 20)

# Test Min Dist 
# -------------

min_dist_c = c(0.002,0.2,0.4,0.6,0.8,0.9)
for (i in 1:length(min_dist_c)){
  Meso.umap = umap(data_lv.D, random_state = 123, min_dist = min_dist_c[i])
  #par(mfrow=c(1,2))
  plot(Meso.umap$layout[,1] ,Meso.umap$layout[,2], main = paste("min_dist =" ,as.character(min_dist_c[i]) ) , col= c("orange", "black", "forestgreen", "brown4")[as.factor(data_lv.L$Type)], xlab = "x" , ylab ="y" , pch=20  )
 # plot(coord_PCA1$x , coord_PCA1$y , col= c("orange", "black", "forestgreen", "brown4")[as.factor(coord_PCA1$Type)], xlab = "x" , ylab ="y" , pch = 20)
}


# Test N Neighbors
# -----------------

n_neighbors_c = c(2,10,20,30,50,80,100,120,130,150,160,170,200,230,250,260)
for (i in 1:length(n_neighbors_c )){
  Meso.umap = umap(data_lv.D, random_state = 123, n_neighbors =  n_neighbors_c[i] )
  #par(mfrow=c(1,2))
  plot(Meso.umap$layout[,1] ,Meso.umap$layout[,2], main = paste("n_neighbors =" ,as.character(n_neighbors_c[i]) ) , col= c("orange", "black", "forestgreen", "brown4")[as.factor(data_lv.L$Type)], xlab = "x" , ylab ="y" , pch=20  )
  #plot(coord_PCA1$x , coord_PCA1$y , col= c("orange", "black", "forestgreen", "brown4")[as.factor(coord_PCA1$Type)], xlab = "x" , ylab ="y" , pch = 20)
}


# Test with N neighbors and min dist
# -----------------------------------

min_dist_c = c(0.2,0.4,0.6,0.8,0.9)
n_neighbors_c = c(10,20,50,80,100,120,170,200,230,260)

for (i in 1:length(min_dist_c)){
  for (j in 1:length(n_neighbors_c )){
    Meso.umap = umap(data_lv.D, random_state = 123, n_neighbors =  n_neighbors_c[j], min_dist = min_dist_c[i] )
   # par(mfrow=c(1,2))
    plot(Meso.umap$layout[,1] ,Meso.umap$layout[,2], main = paste("n_neighbors =" ,as.character(n_neighbors_c[j]), "min_dist = " ,  min_dist_c[i] ) , col= c("orange", "black", "forestgreen", "brown4")[as.factor(data_lv.L$Type)], xlab = "x" , ylab ="y" , pch=20  )
   # plot(coord_PCA1$x , coord_PCA1$y , col= c("orange", "black", "forestgreen", "brown4")[as.factor(coord_PCA1$Type)], xlab = "x" , ylab ="y" , pch = 20)

    }
}


# Write for TumorMap analysis
#----------------------------

# Min dist 

Meso_MD02.umap = umap(data_lv.D, random_state = 123, min_dist =0.2 )
Meso_MD02_df_coord = data.frame("sample"= data_lv.L$sample , "x"=Meso_MD02.umap$layout[,1] , "y"=Meso_MD02.umap$layout[,2] )
#write.table(Meso_MD02_df_coord, file='Meso_MD02_df_coord.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F) # Coordinates2 = coordinates1 with y*-1

Meso_MD09.umap = umap(data_lv.D, random_state = 123, min_dist =0.9 )
Meso_MD09_df_coord = data.frame("sample"= data_lv.L$sample , "x"=Meso_MD09.umap$layout[,1] , "y"=Meso_MD09.umap$layout[,2] )
#write.table(Meso_MD09_df_coord, file='Meso_MD09_df_coord.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F) # Coordinates2 = coordinates1 with y*-1

# Nearest nighbors

Meso_MD02.umap = umap(data_lv.D, random_state = 123, min_dist =0.2 )
Meso_MD02_df_coord = data.frame("sample"= data_lv.L$sample , "x"=Meso_MD02.umap$layout[,1] , "y"=Meso_MD02.umap$layout[,2] )
#write.table(Meso_MD02_df_coord, file='Meso_MD02_df_coord.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F) # Coordinates2 = coordinates1 with y*-1

Meso_MD09.umap = umap(data_lv.D, random_state = 123, min_dist =0.9 )
Meso_MD09_df_coord = data.frame("sample"= data_lv.L$sample , "x"=Meso_MD09.umap$layout[,1] , "y"=Meso_MD09.umap$layout[,2] )
#write.table(Meso_MD09_df_coord, file='Meso_MD09_df_coord.tsv', quote=FALSE, sep='\t', row.names = F , col.names = F) # Coordinates2 = coordinates1 with y*-1



