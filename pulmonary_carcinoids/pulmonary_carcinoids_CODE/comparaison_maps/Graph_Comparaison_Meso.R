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

CP_PCA_TM <- read.table("CP_PCA_TM_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")  # HERE
CP_R_PCA <- read.table("CP_PCA_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_TM <- read.table("CP_meso_TM_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
#CP_R_UMAP_NN150_MD05_nearly_all <- read.table("CP_meso_UMAP_NN150_MD05_real_nearly_all_k.txt", sep = "\t", dec="." , header = TRUE,   quote="")
#CP_R_UMAP_NN150_MD05_Next <- read.table("CP_meso_UMAP_NN150_MD05_real_Next.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_NN150_MD05_all <- rbind(CP_R_UMAP_NN150_MD05_nearly_all ,CP_R_UMAP_NN150_MD05_Next  )
CP_R_UMAP_NN150_MD05 <- CP_R_UMAP_NN150_MD05[CP_R_UMAP_NN150_MD05$K %in% unique(CP_R_TM$K) ,  ]
CP_R_UMAP_NN230 <- read.table("CP_meso_UMAP_NN230_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_NN20 <- read.table("CP_meso_UMAP_NN20_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD09 <- read.table("CP_meso_UMAP_MD09_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD02 <- read.table("CP_meso_UMAP_MD02_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
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


list_CP_df = list( CP_R_PCA, CP_R_TM ,CP_R_UMAP_MD02)#, CP_R_UMAP_NN150_MD05 , CP_R_UMAP_NN230 , CP_R_UMAP_NN20 , CP_R_UMAP_MD09, CP_PCA_TM ,
Name = c("CP_R_PCA" , "CP_R_TM" ,"CP_R_UMAP_MD02")#, "CP_R_UMAP_NN150_MD05", "CP_R_UMAP_NN230" , "CP_R_UMAP_NN20" , "CP_R_UMAP_MD09", "CP_PCA_TM" 
CP_mean_by_k(list_CP_df,Name)

# Set Diff
# ---------
Set_diff_R_PCA_V2 <-  read.table("set_diff_PCA_real_V2.txt", sep = "\t", dec="." , header = TRUE,   quote="")

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



list_df_set= list(Set_diff_R_UMAP_NN230, Set_diff_R_UMAP_NN20,Set_diff_R_UMAP_MD09,Set_diff_R_UMAP_MD02)
#Name_comp = c("PCA_TM", "R_PCA", "R_TM", "R_UMAP_NN150_MD_05")
Name_comp = c("R_UMAP_NN230", "R_UMAP_NN20","R_UMAP_MD09","R_UMAP_MD02")
dist_set(c(20,40,60,100),list_df_set, Name_comp )

list_df_set= list(Set_diff_R_PCA, Set_diff_R_TM,  Set_diff_R_UMAP_NN150_MD_05 )#Set_diff_PCA_TM, Set_diff_R_PCA_V2,  ,
#Set_diff_R_UMAP_NN230, Set_diff_R_UMAP_NN20,Set_diff_R_UMAP_MD09, Set_diff_R_UMAP_MD02
Name = c( "R_PCA", "R_TM","R_UMAP_NN150_MD_05")#"PCA_TM","R_PCA_V2", ,  "R_UMAP_NN230", "R_UMAP_NN20","R_UMAP_MD09","R_UMAP_MD02"
Set_diff_mean_by_k(list_df_set, Name)


# CONCLUSION :  TRES PEU DE DIFFERENCE ENTRE LES MOYENNES

# Seq Diff 
# --------
Seq_diff_R_PCA_V2 <-read.table("seq_diff_PCA_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")


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

list_df_seq= list(Seq_diff_R_PCA  , Seq_diff_R_TM ,  Seq_diff_R_UMAP_MD05_NN150 ) #Seq_diff_PCA_TM ,Seq_diff_R_UMAP_NN230 ,
#Seq_diff_R_UMAP_NN20, Seq_diff_R_UMAP_MD09 ,Seq_diff_R_UMAP_MD02,, Seq_diff_R_PCA_V2

Name= c('R_PCA'  , 'R_TM' ,  'R_UMAP_MD05_NN150') #'R_UMAP_NN230' ,'R_PCA_V2',
#'R_UMAP_NN20', 'R_UMAP_MD09' , 'R_UMAP_MD02',
 
Seq_diff_mean_by_k(list_df_seq,Name)


Overall_df_seq_diff <- data.frame("Sample" = Seq_diff_PCA_TM$sample , "k"= Seq_diff_PCA_TM$k ,'R_PCA'= Seq_diff_R_PCA$seq_diff ,
                                  'R_TM' = Seq_diff_R_TM$seq_diff , 'R_UMAP_NN150_MD_05' =Seq_diff_R_UMAP_MD05_NN150$seq_diff , 'R_UMAP_MD02' = Seq_diff_R_UMAP_MD02$seq_diff,
                                  'R_UMAP_MD09'= Seq_diff_R_UMAP_MD09$seq_diff , 'R_UMAP_NN230' = Seq_diff_R_UMAP_NN230$seq_diff , 'R_UMAP_NN20' = Seq_diff_R_UMAP_NN20$seq_diff)
#  "PCA_TM" = Seq_diff_PCA_TM$seq_diff ,
Overall_df_seq_diff$Sample = as.character(Overall_df_seq_diff$Sample)
Select_sample = unique(Overall_df_seq_diff$Sample)
Select_sample = sample( Select_sample , 8)
Overall_df_seq_diff[Overall_df_seq_diff$Sample %in% Select_sample & Overall_df_seq_diff$k ==20 ,3:dim(Overall_df_seq_diff)[2]]
test1 = data.frame(Overall_df_seq_diff[Overall_df_seq_diff$Sample %in% Select_sample & Overall_df_seq_diff$k ==20 ,3:10])
rownames(test1) = Select_sample

test2 = data.frame(Overall_df_seq_diff[Overall_df_seq_diff$Sample %in% Select_sample & Overall_df_seq_diff$k ==20 ,4])
rownames(test2) = Select_sample

plot( test1$Val  , xaxt="n")  #Select_sample , 
axis(1, at=1:8, labels=rownames(test1),las=2)

View_sampling  <-function (overall_df , K , N, myylab){
  # This function allow to compare methods according their set difference value or sequence differnce values 
  # for |N| different sample, and for 3 K different levels.
  
  Select_sample = unique(overall_df$Sample)
  Select_sample = sample( Select_sample ,N)
  par(mfrow=c(2,2))
  for (i in 1:length(K)){
    to_draw = data.frame(overall_df[overall_df$Sample %in% Select_sample & overall_df$k ==K[i],3:dim(Overall_df_seq_diff)[2]])
    rownames(to_draw)= Select_sample
    Mycol_l = c(1)
    mycol = 1
    for ( j in 1:dim(to_draw)[2]){
      if(j == 1){
        plot(to_draw[,j], xaxt="n", col=mycol, ylab=myylab , main = paste("k = ", K[i]))
        axis(1, at=1:N, labels=rownames(to_draw),las=2)
        mycol =mycol+ 1
        Mycol_l = c(Mycol_l , mycol)
      }
      else{
        points(to_draw[,j], xaxt="n", col=mycol) 
        axis(1, at=1:N, labels=rownames(to_draw),las=2)
        mycol =mycol+ 1
        Mycol_l = c(Mycol_l , mycol)
      }
    }
  }
  Name = colnames(to_draw)
  print(Mycol_l)
  plot.new()
  legend("bottomleft",  legend = Name,  col = c(unique(Mycol_l)), pch = 15)
}

View_sampling(Overall_df_seq_diff, c(20,50,100), 6, "Seq Diff")

