###################################
#       LIBRAIRIES                #
###################################

library(plotly)
library(reticulate)
use_python('/Users/mathian/miniconda2/bin/python2.7')

#############################################################
#     COMPARAISON BETWEEN REAL NEIGHBORHOOD AND OTHERS      #
#############################################################
source("Graph_comp_function.R")

# Coords
# ------


TM_coords  <- read.table("Meso_tm_coords_v2.tab", sep = "\t", dec="." , header = TRUE,   quote="")
colnames(TM_coords)[1] <- "sample"
PCA_coords  <- read.table("Meso_pca_coords.tab", sep = "\t", dec="." , header = TRUE,   quote="")
colnames(PCA_coords)[1] <- "sample"
UMAP_coords_NN230 <- read.table( "umap_coords_nn230.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN230)[1]<- "sample"
UMAP_coords_NN20 <-  read.table( "umap_coords_nn20.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN20)[1]<- "sample"
UMAP_coords_MD09 <-  read.table( "umap_coords_md09.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_MD09)[1]<- "sample"
UMAP_coords_MD02 <-  read.table( "umap_coord_md_02.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_MD02 )[1]<- "sample"
UMAP_coords_NN150_MD05 <- read.table( "umap_coords_nn150_md05.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN150_MD05 )[1]<- "sample"

# Centrality preservation 
# -----------------------

CP_PCA_TM <- read.table("CP_MesosomicsV2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_PCA <- read.table("CP_meso_PCA_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_TM <- read.table("CP_meso_TM_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_NN150_MD05 <- read.table("CP_meso_UMAP_NN150_MD05_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")

CP_R_UMAP_NN230 <- read.table("CP_meso_UMAP_NN230_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_NN20 <- read.table("CP_meso_UMAP_NN20_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD09 <- read.table("CP_meso_UMAP_MD09_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD02 <- read.table("CP_meso_UMAP_MD02_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
par(mfrow=c(1,1))


ku_stack =c(60,170,260,280)
plot_centrality_preservation(CP_PCA_TM, TM_coords, ku_stack  , "2", "CP2 drawn on the 2D projection : Openord & PCA" )
plot_centrality_preservation(CP_PCA_TM, TM_coords, ku_stack  , "N", "CPN drawn on the 2D projection : Openord & PCA" )
par(mfrow=c(2,2))
plot(CP_PCA_TM$CP2[which(CP_PCA_TM$K == 60)], CP_PCA_TM$CPN[which(CP_PCA_TM$K == 60)], xlab="CP_TM",  ylab="CP_PCA")
plot(CP_PCA_TM$CP2[which(CP_PCA_TM$K == 170)], CP_PCA_TM$CPN[which(CP_PCA_TM$K == 170)], xlab="CP_TM",  ylab="CP_PCA")
plot(CP_PCA_TM$CP2[which(CP_PCA_TM$K == 260)], CP_PCA_TM$CPN[which(CP_PCA_TM$K == 260)], xlab="CP_TM",  ylab="CP_PCA")
plot(CP_PCA_TM$CP2[which(CP_PCA_TM$K == 280)], CP_PCA_TM$CPN[which(CP_PCA_TM$K == 280)], xlab="CP_TM",  ylab="CP_PCA")
par(mfrow=c(1,1))


plot_centrality_preservation(CP_R_PCA, PCA_coords, ku_stack  , "2", "CP2 drawn on 2D projection : 'Real' & PCA" )
plot_centrality_preservation(CP_R_PCA, PCA_coords, ku_stack  , "N", "CPN drawn on 2D projection : 'Real' & PCA" )
par(mfrow=c(2,2))
plot(CP_R_PCA$CP2[which(CP_R_PCA$K == 60)], CP_R_PCA$CPN[which(CP_R_PCA$K == 60)], xlab="CP_PCA",  ylab="CPN")
plot(CP_R_PCA$CP2[which(CP_R_PCA$K == 170)], CP_R_PCA$CPN[which(CP_R_PCA$K == 170)], xlab="CP_PCA",  ylab="CPN")
plot(CP_R_PCA$CP2[which(CP_R_PCA$K == 280)], CP_R_PCA$CPN[which(CP_R_PCA$K == 280)], xlab="CP_PCA",  ylab="CPN")
par(mfrow=c(1,1))

plot_centrality_preservation(CP_R_TM, TM_coords, ku_stack  , "2", "CP2 drawn on 2D projection : Real & OpenOrd ")
plot_centrality_preservation(CP_R_TM, TM_coords, ku_stack  , "N", "CPN preservation drawn on 2D projection : Real & OpenOrd ")
par(mfrow=c(2,2))
plot(CP_R_TM$CP2[which(CP_R_TM$K == 60)], CP_R_TM$CPN[which(CP_R_TM$K == 60)], xlab="CP_TM",  ylab="CPN")
plot(CP_R_TM$CP2[which(CP_R_TM$K == 170)], CP_R_TM$CPN[which(CP_R_TM$K == 170)], xlab="CP_TM",  ylab="CPN")
plot(CP_R_TM$CP2[which(CP_R_TM$K == 280)], CP_R_TM$CPN[which(CP_R_TM$K == 280)], xlab="CP_TM",  ylab="CPN")
par(mfrow=c(1,1))

plot_centrality_preservation(CP_R_UMAP_NN150_MD05, UMAP_coords_NN150_MD05, ku_stack  , "2", "CP2 drawn on 2D projection : Real & Umap nn=150 md=0.5  ")
plot_centrality_preservation(CP_R_UMAP_NN150_MD05, UMAP_coords_NN150_MD05, ku_stack  , "N", "CPN drawn on 2D projection : Real & Umap nn=150 md=0.5  ")
par(mfrow=c(2,2))
plot(CP_R_UMAP_NN150_MD05$CP2[which(CP_R_UMAP_NN150_MD05$K == 60)], CP_R_UMAP_NN150_MD05$CPN[which(CP_R_UMAP_NN150_MD05$K == 60)], xlab="CP_UMAP_NN120_MD_05",  ylab="CPN")
plot(CP_R_UMAP_NN150_MD05$CP2[which(CP_R_UMAP_NN150_MD05$K == 170)], CP_R_UMAP_NN150_MD05$CPN[which(CP_R_UMAP_NN150_MD05$K == 170)], xlab="CP_UMAP_NN120_MD_05",  ylab="CPN")
plot(CP_R_UMAP_NN150_MD05$CP2[which(CP_R_UMAP_NN150_MD05$K == 280)], CP_R_UMAP_NN150_MD05$CPN[which(CP_R_UMAP_NN150_MD05$K == 280)], xlab="CP_UMAP_NN120_MD_05",  ylab="CPN")
par(mfrow=c(1,1))

plot_centrality_preservation(CP_R_UMAP_NN230, UMAP_coords_NN230, ku_stack  , "2", "CP2 drawn on 2D projection : Real & Umap nn=230  ")
plot_centrality_preservation(CP_R_UMAP_NN230, UMAP_coords_NN230, ku_stack  , "N", "CPN drawn on 2D projection : Real & Umap nn=230  ")
par(mfrow=c(2,2))
plot(CP_R_UMAP_NN230$CP2[which(CP_R_UMAP_NN230$K == 60)], CP_R_UMAP_NN230$CPN[which(CP_R_UMAP_NN230$K == 60)], xlab="CP_UMAP_NN230",  ylab="CPN")
plot(CP_R_UMAP_NN230$CP2[which(CP_R_UMAP_NN230$K == 170)], CP_R_UMAP_NN230$CPN[which(CP_R_UMAP_NN230$K == 170)], xlab="CP_UMAP_NN230",  ylab="CPN")
plot(CP_R_UMAP_NN230$CP2[which(CP_R_UMAP_NN230$K == 280)], CP_R_UMAP_NN230$CPN[which(CP_R_UMAP_NN230$K == 280)], xlab="CP_UMAP_NN230",  ylab="CPN")
par(mfrow=c(1,1))

plot_centrality_preservation(CP_R_UMAP_NN20, UMAP_coords_NN20, ku_stack  , "2", "CP2 drawn on 2D projection : Real & Umap nn=20  ")
plot_centrality_preservation(CP_R_UMAP_NN20, UMAP_coords_NN20, ku_stack  , "N", "CPN drawn on 2D projection : Real & Umap nn=20  ")
par(mfrow=c(2,2))
plot(CP_R_UMAP_NN20$CP2[which(CP_R_UMAP_NN20$K == 60)], CP_R_UMAP_NN20$CPN[which(CP_R_UMAP_NN20$K == 60)], xlab="CP_UMAP_NN20",  ylab="CPN")
plot(CP_R_UMAP_NN20$CP2[which(CP_R_UMAP_NN20$K == 170)], CP_R_UMAP_NN20$CPN[which(CP_R_UMAP_NN20$K == 170)], xlab="CP_UMAP_NN20",  ylab="CPN")
plot(CP_R_UMAP_NN20$CP2[which(CP_R_UMAP_NN20$K == 280)], CP_R_UMAP_NN20$CPN[which(CP_R_UMAP_NN20$K == 280)], xlab="CP_UMAP_NN20",  ylab="CPN")
par(mfrow=c(1,1))

plot_centrality_preservation(CP_R_UMAP_MD09, UMAP_coords_MD09, ku_stack  , "2", "CP2 drawn on 2D projection : Real & Umap md=0.9  ")
plot_centrality_preservation(CP_R_UMAP_MD09, UMAP_coords_MD09, ku_stack  , "N", "CPN drawn on 2D projection : Real & Umap md=0.9  ")
par(mfrow=c(2,2))
plot(CP_R_UMAP_MD09$CP2[which(CP_R_UMAP_MD09$K == 60)], CP_R_UMAP_MD09$CPN[which(CP_R_UMAP_MD09$K == 60)], xlab="CP_UMAP_MD09",  ylab="CPN")
plot(CP_R_UMAP_MD09$CP2[which(CP_R_UMAP_MD09$K == 170)], CP_R_UMAP_MD09$CPN[which(CP_R_UMAP_MD09$K == 170)], xlab="CP_UMAP_MD09",  ylab="CPN")
plot(CP_R_UMAP_MD09$CP2[which(CP_R_UMAP_MD09$K == 280)], CP_R_UMAP_MD09$CPN[which(CP_R_UMAP_MD09$K == 280)], xlab="CP_UMAP_MD09",  ylab="CPN")
par(mfrow=c(1,1))


plot_centrality_preservation(CP_R_UMAP_MD02, UMAP_coords_MD02, ku_stack  , "2", "CP2 drawn on 2D projection : Real & Umap md=0.2  ")
plot_centrality_preservation(CP_R_UMAP_MD02, UMAP_coords_MD02, ku_stack  , "N", "CPN drawn on 2D projection : Real & Umap md=0.2  ")

par(mfrow=c(2,2))
plot(CP_R_UMAP_MD02$CP2[which(CP_R_UMAP_MD02$K == 60)], CP_R_UMAP_MD02$CPN[which(CP_R_UMAP_MD02$K == 60)], xlab="CP_UMAP_MD02",  ylab="CPN")
plot(CP_R_UMAP_MD02$CP2[which(CP_R_UMAP_MD02$K == 170)], CP_R_UMAP_MD02$CPN[which(CP_R_UMAP_MD02$K == 170)], xlab="CP_UMAP_MD02",  ylab="CPN")
plot(CP_R_UMAP_MD02$CP2[which(CP_R_UMAP_MD02$K == 280)], CP_R_UMAP_MD02$CPN[which(CP_R_UMAP_MD02$K == 280)], xlab="CP_UMAP_MD02",  ylab="CPN")
par(mfrow=c(1,1))


list_CP_df = list(CP_PCA_TM , CP_R_PCA , CP_R_TM , CP_R_UMAP_NN150_MD05 , CP_R_UMAP_NN230 , CP_R_UMAP_NN20 , CP_R_UMAP_MD09, CP_R_UMAP_MD02)
Name = c("CP_PCA_TM" , "CP_R_PCA" , "CP_R_TM", "CP_R_UMAP_NN150_MD05" , "CP_R_UMAP_NN230" , "CP_R_UMAP_NN20" , "CP_R_UMAP_MD09", "CP_R_UMAP_MD02")
CP_mean_by_k(list_CP_df,Name)

# Set Diff
# ---------

Set_diff_PCA_TM <- read.table("set_diff_meso.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_PCA <-  read.table("set_diff_meso_PCA_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_TM <-  read.table("set_diff_meso_TM_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_UMAP_NN150_MD_05 <-  read.table("set_diff_meso_UMAP_NN150_MD05_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_UMAP_NN230 <-  read.table("set_diff_meso_UMAP_NN230_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_UMAP_NN20 <-  read.table("set_diff_meso_UMAP_NN20_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_UMAP_MD09 <-  read.table("set_diff_meso_UMAP_MD09_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Set_diff_R_UMAP_MD02 <-  read.table("set_diff_meso_UMAP_MD02_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")

plot_set_diff(Set_diff_PCA_TM, TM_coords , c(5,10,15,20,50,150,200,250) , "Set diff PCA TM")
plot_set_diff(Set_diff_R_PCA, PCA_coords , c(5,10,15,20,50,150,200,250) , "Set diff R PCA")
plot_set_diff(Set_diff_R_TM, TM_coords , c(5,10,15,20,50,150,200,250) , "Set diff R TM")
plot_set_diff(Set_diff_R_UMAP_NN150_MD_05, UMAP_coords_NN150_MD05 , c(5,10,15,20,50,150,200,250) , "Set diff R UMAP_NN150_MD05")
plot_set_diff(Set_diff_R_UMAP_NN230, UMAP_coords_NN230 , c(5,10,15,20,50,150,200,250) , "Set diff R UMAP_NN230")
plot_set_diff(Set_diff_R_UMAP_NN20, UMAP_coords_NN20 , c(5,10,15,20,50,150,200,250) , "Set diff R UMAP_NN20")
plot_set_diff(Set_diff_R_UMAP_MD09, UMAP_coords_MD09 , c(5,10,15,20,50,150,200,250) , "Set diff R UMAP_MD09")
plot_set_diff(Set_diff_R_UMAP_MD02, UMAP_coords_MD02 , c(5,10,15,20,50,150,200,250) , "Set diff R UMAP_MD02")


library(ggplot2)

#boxplot(Set_diff_PCA_TM$set_diff[Set_diff_PCA_TM$k==20], Set_diff_R_PCA$set_diff[Set_diff_R_PCA$k==20],
#        Set_diff_R_TM$set_diff[Set_diff_R_TM$k==20], Set_diff_R_UMAP_NN150_MD_05$set_diff[Set_diff_R_UMAP_NN150_MD_05$k==20],
#        Set_diff_R_UMAP_NN20$set_diff[Set_diff_R_UMAP_NN20$k==20])

dist_set <- function(K , list_data_frame , Name ){
  
  par(mfrow= c(length(list_data_frame),length(K))) 
  c = 1
   for (j in 1:length(list_data_frame )){
    for (i in 1:length(K)) {
        c_data_frame = as.data.frame(list_data_frame[j])
        plot(c_data_frame$set_diff[c_data_frame$k==K[i]], col = c, main = paste("k = ", K[i] , Name[j]),  xlab="Index", ylab="set diff value")
      }
    c =c+1
  }  
}

#list_df_set= list(Set_diff_PCA_TM, Set_diff_R_PCA, Set_diff_R_TM, Set_diff_R_UMAP_NN150_MD_05)
list_df_set= list(Set_diff_R_UMAP_NN230, Set_diff_R_UMAP_NN20,Set_diff_R_UMAP_MD09,Set_diff_R_UMAP_MD02)
#Name_comp = c("PCA_TM", "R_PCA", "R_TM", "R_UMAP_NN150_MD_05")
Name_comp = c("R_UMAP_NN230", "R_UMAP_NN20","R_UMAP_MD09","R_UMAP_MD02")
dist_set(c(20,40,60,100),list_df_set, Name_comp )

plot(x=unique(Set_diff_PCA_TM$k) ,y= tapply(Set_diff_PCA_TM$set_diff ,Set_diff_PCA_TM$k, mean))

Set_diff_mean_by_k  <-function (list_set_diff , Name){
  par(mfrow=c(1,2))
  c= 1 
  c_data_frame = as.data.frame(list_set_diff[1])
  Set_diff_mean_df = data.frame("k" =unique(c_data_frame$k ))
  for (i in 1:length(list_set_diff)){
    c_data_frame = as.data.frame(list_set_diff[i])
    Mean_by_k =tapply(c_data_frame$set_diff ,c_data_frame$k, mean)
    Set_diff_mean_df <- cbind(Set_diff_mean_df,Mean_by_k )
    colnames(Set_diff_mean_df)[dim(Set_diff_mean_df)[2]] <- Name[i]
   # print(head(Set_diff_mean_df))
    
  }
  Col = c()
  C=1
  for (i in 2:dim(Set_diff_mean_df)[2]){
   
    Col = c(Col,C)
    if (i ==2 ){
        plot(Set_diff_mean_df$k, Set_diff_mean_df[,i] , col=C  , type ='l' )
      C =C+1
      Col = c(Col,C)
    }
    else{
      lines(Set_diff_mean_df$k, Set_diff_mean_df[,i] , col=C )
     
     
      C =C+1
      Col = c(Col,C)
    }
  
  }
  print(unique(Col))
  plot.new()
  legend("bottomleft", 
         legend = Name, 
         col = c(unique(Col)),
         pch = 15)


}
    
    
 

list_df_set= list(Set_diff_R_PCA, Set_diff_R_TM, Set_diff_R_UMAP_NN150_MD_05 ,Set_diff_PCA_TM, 
                  Set_diff_R_UMAP_NN230, Set_diff_R_UMAP_NN20,Set_diff_R_UMAP_MD09, Set_diff_R_UMAP_MD02)
Name = c( "R_PCA", "R_TM", "R_UMAP_NN150_MD_05","PCA_TM",  "R_UMAP_NN230", "R_UMAP_NN20","R_UMAP_MD09","R_UMAP_MD02")
testP <- Set_diff_mean_by_k(list_df_set, Name)
testP 

# CONCLUSION :  TRES PEU DE DIFFERENCE ENTRE LES MOYENNES

# Seq Diff 
# --------

Seq_diff_PCA_TM <- read.table("seq_diff_meso.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_PCA <- read.table("seq_diff_meso_PCA_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_TM <- read.table("seq_diff_meso_TM_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_UMAP_NN230 <- read.table("seq_diff_NN230_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_UMAP_NN20 <- read.table("seq_diff_NN20_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_UMAP_MD09 <- read.table("seq_diff_MD09_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_UMAP_MD02 <- read.table("seq_diff_MD02_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff_R_UMAP_MD05_NN150 <- read.table("seq_diff_NN150_MD05_real.txt", sep = "\t", dec="." , header = TRUE,   quote="")


plot_seq_diff(Seq_diff_PCA_TM, TM_coords , c(20,40,80,120,160,200,240,280) , "Seq diff PCA TM")
plot_seq_diff(Seq_diff_R_PCA, PCA_coords , c(20,40,80,120,160,200,240,280) , "Seq diff R PCA")
plot_seq_diff(Seq_diff_R_TM , TM_coords , c(20,40,80,120,160,200,240,280) , "Seq diff R TM")
plot_seq_diff(Seq_diff_R_UMAP_NN230 , UMAP_coords_NN230 , c(20,40,80,120,160,200,240,280) , "Seq diff R UMAP NN 230")
plot_seq_diff(Seq_diff_R_UMAP_NN20 , UMAP_coords_NN20 , c(20,40,80,120,160,200,240,280) , "Seq diff R UMAP NN 20")
plot_seq_diff(Seq_diff_R_UMAP_MD09 , UMAP_coords_MD09 , c(20,40,80,120,160,200,240,280) , "Seq diff R UMAP MD 09")
plot_seq_diff(Seq_diff_R_UMAP_MD02 , UMAP_coords_MD02 , c(20,40,80,120,160,200,240,280) , "Seq diff R UMAP MD 02")
plot_seq_diff(Seq_diff_R_UMAP_MD05_NN150 , UMAP_coords_NN150_MD05 , c(20,40,80,120,160,200,240,280) , "Seq diff R UMAP NN 150 MD 05")


plot(Seq_diff_R_UMAP_MD05_NN150$seq_diff)
plot(Seq_diff_R_TM$seq_diff)



Seq_diff_mean_by_k  <-function (list_Seq_diff , Name){
  par(mfrow=c(1,2))
  c= 1 
  c_data_frame = as.data.frame(list_Seq_diff[1])
  Seq_diff_mean_df = data.frame("k" =unique(c_data_frame$k ))
  for (i in 1:length(list_Seq_diff)){
    c_data_frame = as.data.frame(list_Seq_diff[i])
    Mean_by_k =tapply(c_data_frame$seq_diff ,c_data_frame$k, mean)
    Seq_diff_mean_df <- cbind(Seq_diff_mean_df,Mean_by_k )
    colnames(Seq_diff_mean_df)[dim(Seq_diff_mean_df)[2]] <- Name[i]
    # print(head(Seq_diff_mean_df))
    
  }
  Col = c()
  C=1
  for (i in 2:dim(Seq_diff_mean_df)[2]){
    
    Col = c(Col,C)
    if (i ==2 ){
      plot(Seq_diff_mean_df$k, Seq_diff_mean_df[,i] , col=C  , type ='l' , )
      C =C+1
      Col = c(Col,C)
    }
    else{
      lines(Seq_diff_mean_df$k, Seq_diff_mean_df[,i] , col=C )
      C =C+1
      Col = c(Col,C)
    }
    
  }
  print(unique(Col))
  plot.new()
  legend("bottomleft", 
         legend = Name, 
         col = c(unique(Col)),
         pch = 15)
}

list_df_seq= list( Seq_diff_PCA_TM ,Seq_diff_R_PCA ,Seq_diff_R_TM , Seq_diff_R_UMAP_NN230 ,
                   Seq_diff_R_UMAP_NN20, Seq_diff_R_UMAP_MD09 ,Seq_diff_R_UMAP_MD02,Seq_diff_R_UMAP_MD05_NN150 )

Name= c( 'R_PCA' , 'R_TM' , 'R_UMAP_NN230' ,
          'R_UMAP_NN20', 'R_UMAP_MD09' , 'R_UMAP_MD02','R_UMAP_MD05_NN150')
#'PCA_TM' , 
Seq_diff_mean_by_k(list_df_seq,Name)
        