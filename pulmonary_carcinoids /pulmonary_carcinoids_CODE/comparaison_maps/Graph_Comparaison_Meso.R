###################################
#       LIBRAIRIES                #
###################################

library(plotly)
library(reticulate)
use_python('/Users/mathian/miniconda2/bin/python2.7')

###################################
#       DATA                      #
###################################
TM_coords  <- read.table("Meso_tm_coords_v2.tab", sep = "\t", dec="." , header = TRUE,   quote="")
colnames(TM_coords)[1] <- "sample"
PCA_coords  <- read.table("Meso_pca_coords.tab", sep = "\t", dec="." , header = TRUE,   quote="")

Set_diff <- read.table("set_diff_meso.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff <- read.table("seq_diff_meso.txt", sep = "\t", dec="." , header = TRUE,   quote="")
Seq_diff$seq_diff <- as.numeric(Seq_diff$seq_diff)

###################################
#       SCRIPT PYTHON             #
###################################
source_python("jaccard_set_distance.py")


###################################
# SET DIFFERENCE VIEW             #
###################################

# K = 30
# ------
set_diff_k30 = Set_diff[8237:8520,]
Set_diff_with_tm <- merge(set_diff_k30,TM_coords, by="sample" )
colnames(Set_diff_with_tm)[2] = "SET_DIFFERENCE"
ATitle =as.character(  paste("k = ", as.character(30)) )
p<- plot_ly(Set_diff_with_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
        marker=list( size=20 , opacity=0.5), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample)) %>%
  layout(annotations=  list(x = 0.2 , y = 1.05, text = ATitle, showarrow = F, xref='paper', yref='paper'))
#colorbar(p, limits = c(0, 1))

p1 <- plot_ly(economics, x = ~date, y = ~unemploy) %>%
  add_lines(name = ~"unemploy")
p2 <- plot_ly(economics, x = ~date, y = ~uempmed) %>%
  add_lines(name = ~"uempmed")
#p <- subplot(p1, p2)


# Fonction for set difference visualisation 
# __________________________________________

ku_stack = c(unique(Set_diff$k))
ku_stack =c(6,20,30,50,100, 110 , 130 , 120 ,150, 180 ,200,250)

#ku_stack =c(10,20,50,100)
while (length(ku_stack)!=0) { # after k_list size
  n_node = length(which(Set_diff$k ==1))
  set_diff_k1st = Set_diff[(n_node*ku_stack[1] - n_node  + 1):(n_node *ku_stack[1]),]
  Set_diff_k1st_tm <- merge(set_diff_k1st,TM_coords, by="sample" )
  colnames(Set_diff_k1st_tm)[2] = "SET_DIFFERENCE"
  Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
  p1<- plot_ly(Set_diff_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
          marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample)) %>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper') , showlegend = FALSE)
  p1<- colorbar(p1, limits = c(0, 1))

  
  set_diff_k2nd = Set_diff[(n_node*ku_stack[2] - n_node  + 1):(n_node *ku_stack[2]),]
  Set_diff_k2nd_tm <- merge(set_diff_k2nd,TM_coords, by="sample" )
  colnames(Set_diff_k2nd_tm)[2] = "SET_DIFFERENCE"
  Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
  p2<- plot_ly(Set_diff_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper') ,  showlegend = FALSE)
  p2<- colorbar(p2, limits = c(0, 1))

  set_diff_k3rd = Set_diff[(n_node*ku_stack[3] - n_node  + 1):(n_node *ku_stack[3]),]
  Set_diff_k3rd_tm <- merge(set_diff_k3rd ,TM_coords, by="sample" )
  colnames(Set_diff_k3rd_tm)[2] = "SET_DIFFERENCE"
  Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
  p3<- plot_ly(Set_diff_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper') ,  showlegend = FALSE)
  p3<- colorbar(p3, limits = c(0, 1)) 
 
  
  set_diff_k4th = Set_diff[(n_node*ku_stack[4] - n_node  + 1):(n_node *ku_stack[4]),]
  
  Set_diff_k4th_tm <- merge(set_diff_k4th, TM_coords ,by="sample" )

  colnames(Set_diff_k4th_tm)[2] = "SET_DIFFERENCE"
  Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
  p4<- plot_ly(Set_diff_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  p4<- colorbar(p4, limits = c(0, 1)) 
  
  
  ku_stack <- ku_stack[-c(1:4)]
  p <- subplot(p1, p2, p3, p4)
  
  print(p)
}



Set_diff$k = as.factor(Set_diff$k)
Mean_set_diff =tapply(Set_diff$set_diff, Set_diff$k, mean)
Mean_set_diff = data.frame("k"= seq(1:284), "Mean_set_diff" =Mean_set_diff )
p_mean <- plot_ly(Mean_set_diff, x = ~k, y = ~Mean_set_diff, type = 'scatter', mode = 'lines')
p_mean



Set_diff_M100PT <- Set_diff[c(which(Set_diff$sample == "M100PT")), ]
plot(Set_diff_M100PT$k, Set_diff_M100PT$set_diff , xlab="k" , ylab='set_diff' , title='Set diff for M1OOPT for all k')



# Fonction for seq difference visualisation 
# __________________________________________


min(Seq_diff$seq_diff)
max(Seq_diff$seq_diff)
summary(Seq_diff$seq_diff)[2]
summary(Seq_diff$seq_diff)[3]
ku_stack = c(unique(Seq_diff$k)) ; ku_stack
ku_stack =c(6,20,30,50,100, 110 , 130 , 120 ,150, 180 ,200,250)




#min(Seq_diff$seq_diff)
#max(Seq_diff$seq_diff)
#plot(Seq_diff$seq_diff)
binf <- min(Seq_diff$seq_diff)
bsup <- max(Seq_diff$seq_diff)

#ku_stack =c(10,20,50,100)
while (length(ku_stack)!=0) { # after k_list size
  n_node = length(which(Seq_diff$k ==1))
  
  seq_diff_k1st = Seq_diff[(n_node*ku_stack[1] - n_node  + 1):(n_node *ku_stack[1]),]
  print(head(seq_diff_k1st))
  Seq_diff_k1st_tm <- merge(seq_diff_k1st, TM_coords ,by="sample" )
  colnames(Seq_diff_k1st_tm)[2] = "SEQ_DIFFERENCE"
  Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
  p1_seq<- plot_ly(Seq_diff_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p1_seq<- colorbar(p1_seq, limits = c(binf , bsup) )
  
  seq_diff_k2nd = Seq_diff[(n_node*ku_stack[2] - n_node  + 1):(n_node *ku_stack[2]),]
  Seq_diff_k2nd_tm <- merge(seq_diff_k2nd, TM_coords ,by="sample" )
  colnames(Seq_diff_k2nd_tm)[2] = "SEQ_DIFFERENCE"
  Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
  p2_seq<- plot_ly(Seq_diff_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p2_seq<- colorbar(p2_seq, limits = c(binf, bsup) )
  
  seq_diff_k3rd = Seq_diff[(n_node*ku_stack[3] - n_node  + 1):(n_node *ku_stack[3]),]
  Seq_diff_k3rd_tm <- merge(seq_diff_k3rd, TM_coords ,by="sample" )
  colnames(Seq_diff_k3rd_tm)[2] = "SEQ_DIFFERENCE"
  Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
  p3_seq<- plot_ly(Seq_diff_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #3_seq <- colorbar(p3_seq, limits = c(binf, bsup ) )
  
  
  seq_diff_k4th = Seq_diff[(n_node*ku_stack[4] - n_node  + 1):(n_node *ku_stack[4]),]
  Seq_diff_k4th_tm <- merge(seq_diff_k4th, TM_coords ,by="sample" )
  colnames(Seq_diff_k4th_tm)[2] = "SEQ_DIFFERENCE"
  Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
  p4_seq<- plot_ly(Seq_diff_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p4_seq<- colorbar(p4_seq, limits = c(binf,bsup) )
  
  ku_stack <- ku_stack[-c(1:4)]
  p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)
  
  print(p_seq)

}


##############################################
#     CENTRALITY PRESERVATION                #
##############################################

CP_all <- read.table("CP_MesosomicsV2.txt", sep = "\t", dec="." , header = TRUE,   quote="")

#CP_all <- read.table("CP_Mesosomics.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_2 <- data.frame("sample" = CP_all$sample ,"CP2"= CP_all$CP2 ,"K"= CP_all$K)
CP_2$CP2 <- as.numeric(CP_2$CP2)
ku_stack = unique(CP_all$K) ; ku_stack
#ku_stack =c(20,50,80,100 ,120,150,200,230)


binf <- min(CP_2$CP2)
bsup <- max(CP_2$CP2)

#ku_stack =c(10,20,50,100)
while (length(ku_stack)!=0) { # after k_list size
 
  CP_2_k1st = CP_2[min(which(CP_2$K == ku_stack[1])):max(which(CP_2$K == ku_stack[1])),]
  
  CP_2_k1st_tm <- merge(CP_2_k1st, TM_coords ,by="sample" )
  colnames(CP_2_k1st_tm)[2] = "CP_2"
  Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
  p1_seq<- plot_ly(CP_2_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p1_seq<- colorbar(p1_seq, limits = c(binf , bsup) )
  
  CP_2_k2nd = CP_2[min(which(CP_2$K == ku_stack[2])):max(which(CP_2$K == ku_stack[2])),]
  CP_2_k2nd_tm <- merge(CP_2_k2nd, TM_coords ,by="sample" )
  colnames(CP_2_k2nd_tm)[2] = "CP_2"
  Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
  p2_seq<- plot_ly(CP_2_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p2_seq<- colorbar(p2_seq, limits = c(binf, bsup) )
  
  CP_2_k3rd = CP_2[min(which(CP_2$K == ku_stack[3])):max(which(CP_2$K == ku_stack[3])),]
  CP_2_k3rd_tm <- merge(CP_2_k3rd, TM_coords ,by="sample" )
  colnames(CP_2_k3rd_tm)[2] = "CP_2"
  Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
  p3_seq<- plot_ly(CP_2_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
 # p3_seq <- colorbar(p3_seq, limits = c(binf, bsup ) )
  
  
  CP_2_k4th = CP_2[min(which(CP_2$K == ku_stack[4])):max(which(CP_2$K == ku_stack[4])),]
  CP_2_k4th_tm <- merge(CP_2_k4th, TM_coords ,by="sample" )
  colnames(CP_2_k4th_tm)[2] = "CP_2"
  Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
  p4_seq<- plot_ly(CP_2_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE)
  #p4_seq<- colorbar(p4_seq, limits = c(binf,bsup) )
  
  ku_stack <- ku_stack[-c(1:4)]
  p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)%>%
    layout(title = "centrality preservation",  margin = 0.04)
  
  print(p_seq)
  
}

