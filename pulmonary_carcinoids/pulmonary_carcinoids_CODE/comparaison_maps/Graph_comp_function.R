library(plotly)
library(reticulate)
use_python('/Users/mathian/miniconda2/bin/python2.7')

plot_set_diff <- function(Set_diff, Coords_df ,ku_stack , mytitle){
#""" Print set_diif plot according  k value 4 by 4"""
  while (length(ku_stack)!=0) { # after k_list size
    n_node = length(which(Set_diff$k ==Set_diff$k[1]))
    set_diff_k1st = Set_diff[(n_node*ku_stack[1] - n_node  + 1):(n_node *ku_stack[1]),]
    Set_diff_k1st_tm <- merge(set_diff_k1st,Coords_df, by="sample" )
    colnames(Set_diff_k1st_tm)[2] = "SET_DIFFERENCE"
    Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
  
    p1<- plot_ly(Set_diff_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample)) %>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper') , showlegend = FALSE)
    p1<- colorbar(p1, limits = c(0, 1))
  
  
    set_diff_k2nd = Set_diff[(n_node*ku_stack[2] - n_node  + 1):(n_node *ku_stack[2]),]
    Set_diff_k2nd_tm <- merge(set_diff_k2nd,Coords_df, by="sample" )
    colnames(Set_diff_k2nd_tm)[2] = "SET_DIFFERENCE"
    Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
    p2<- plot_ly(Set_diff_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper') ,  showlegend = FALSE)
    p2<- colorbar(p2, limits = c(0, 1))
  
    set_diff_k3rd = Set_diff[(n_node*ku_stack[3] - n_node  + 1):(n_node *ku_stack[3]),]
    Set_diff_k3rd_tm <- merge(set_diff_k3rd ,Coords_df, by="sample" )
    colnames(Set_diff_k3rd_tm)[2] = "SET_DIFFERENCE"
    Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
    p3<- plot_ly(Set_diff_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper') ,  showlegend = FALSE)
    p3<- colorbar(p3, limits = c(0, 1)) 
  
  
    set_diff_k4th = Set_diff[(n_node*ku_stack[4] - n_node  + 1):(n_node *ku_stack[4]),]
    Set_diff_k4th_tm <- merge(set_diff_k4th, Coords_df ,by="sample" )
    colnames(Set_diff_k4th_tm)[2] = "SET_DIFFERENCE"
    Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
    p4<- plot_ly(Set_diff_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
               marker=list( size=10 , opacity=1), color = ~SET_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
    p4<- colorbar(p4, limits = c(0, 1)) 
  
  
    ku_stack <- ku_stack[-c(1:4)]
    p <- subplot(p1, p2, p3, p4)%>%
    layout(title = mytitle,  margin = 0.04)
  
   print(p)
  }
  
} 




plot_seq_diff <- function(Seq_diff, Coords_df ,ku_stack, mytitle ){
# Representation of centrality preservtion calculated on "2" or "N"
  while (length(ku_stack)!=0) { # after k_list size
    n_node = length(which(Seq_diff$k ==Seq_diff$k[1]))
    seq_diff_k1st = Seq_diff[(n_node*ku_stack[1] - n_node  + 1):(n_node *ku_stack[1]),]
    Seq_diff_k1st_tm <- merge(seq_diff_k1st, Coords_df ,by="sample" )
    colnames(Seq_diff_k1st_tm)[2] = "SEQ_DIFFERENCE"
    Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
    p1_seq<- plot_ly(Seq_diff_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers",
                     marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p1_seq<- colorbar(p1_seq, limits = c(binf , bsup) )
  
  seq_diff_k2nd = Seq_diff[(n_node*ku_stack[2] - n_node  + 1):(n_node *ku_stack[2]),]
  Seq_diff_k2nd_tm <- merge(seq_diff_k2nd, Coords_df ,by="sample" )
  colnames(Seq_diff_k2nd_tm)[2] = "SEQ_DIFFERENCE"
  Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
  p2_seq<- plot_ly(Seq_diff_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p2_seq<- colorbar(p2_seq, limits = c(binf, bsup) )
  
  seq_diff_k3rd = Seq_diff[(n_node*ku_stack[3] - n_node  + 1):(n_node *ku_stack[3]),]
  Seq_diff_k3rd_tm <- merge(seq_diff_k3rd, Coords_df ,by="sample" )
  colnames(Seq_diff_k3rd_tm)[2] = "SEQ_DIFFERENCE"
  Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
  p3_seq<- plot_ly(Seq_diff_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #3_seq <- colorbar(p3_seq, limits = c(binf, bsup ) )
  
  
  seq_diff_k4th = Seq_diff[(n_node*ku_stack[4] - n_node  + 1):(n_node *ku_stack[4]),]
  Seq_diff_k4th_tm <- merge(seq_diff_k4th, Coords_df ,by="sample" )
  colnames(Seq_diff_k4th_tm)[2] = "SEQ_DIFFERENCE"
  Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
  p4_seq<- plot_ly(Seq_diff_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                   marker=list( size=10 , opacity=1), color = ~SEQ_DIFFERENCE, text = ~paste('Sample: ', sample))%>%
    layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
  #p4_seq<- colorbar(p4_seq, limits = c(binf,bsup) )
  
  ku_stack <- ku_stack[-c(1:4)]
  p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)%>%
    layout(title = mytitle,  margin = 0.04)
  
  print(p_seq)
  
}
  
}



plot_centrality_preservation <- function(CP,  Coords_df ,ku_stack, projection , mytitle){
  if (projection == "2"){
    CP_2 <- data.frame("sample" = CP$sample ,"CP2"= CP$CP2 ,"K"= CP$K)
    CP_2$CP2 <- as.numeric(CP_2$CP2)
    print(CP_2$CP2[1:10])
    while (length(ku_stack)!=0) { # after k_list size
    
      CP_2_k1st = CP_2[min(which(CP_2$K == ku_stack[1])):max(which(CP_2$K == ku_stack[1])),]
      print(122)
      CP_2_k1st_tm <- merge(CP_2_k1st, Coords_df ,by="sample" )
      print(124)
      colnames(CP_2_k1st_tm)[2] = "CP_2"
      Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
      print(126)
      p1_seq<- plot_ly(CP_2_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                     marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
      layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      #p1_seq<- colorbar(p1_seq, limits = c(binf , bsup) )
      print(132)
      CP_2_k2nd = CP_2[min(which(CP_2$K == ku_stack[2])):max(which(CP_2$K == ku_stack[2])),]
      CP_2_k2nd_tm <- merge(CP_2_k2nd, Coords_df ,by="sample" )
      colnames(CP_2_k2nd_tm)[2] = "CP_2"
      Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
      p2_seq<- plot_ly(CP_2_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                     marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
      layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      #p2_seq<- colorbar(p2_seq, limits = c(binf, bsup) )
    
      CP_2_k3rd = CP_2[min(which(CP_2$K == ku_stack[3])):max(which(CP_2$K == ku_stack[3])),]
      CP_2_k3rd_tm <- merge(CP_2_k3rd, Coords_df ,by="sample" )
      colnames(CP_2_k3rd_tm)[2] = "CP_2"
      Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
      p3_seq<- plot_ly(CP_2_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                     marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
      layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      # p3_seq <- colorbar(p3_seq, limits = c(binf, bsup ) )
    
    
      CP_2_k4th = CP_2[min(which(CP_2$K == ku_stack[4])):max(which(CP_2$K == ku_stack[4])),]
      CP_2_k4th_tm <- merge(CP_2_k4th, Coords_df ,by="sample" )
      colnames(CP_2_k4th_tm)[2] = "CP_2"
      Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
      p4_seq<- plot_ly(CP_2_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                     marker=list( size=10 , opacity=1), color = ~CP_2, text = ~paste('Sample: ', sample))%>%
      layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE)
      #p4_seq<- colorbar(p4_seq, limits = c(binf,bsup) )
    
      ku_stack <- ku_stack[-c(1:4)]
      p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)%>%
      layout(title = mytitle,  margin = 0.04)
    
      print(p_seq)
   }
    
  }
  else{
    
    CP_N <- data.frame("sample" = CP$sample ,"CPN"= CP$CPN ,"K"= CP$K)
    CP_N$CPN <- as.numeric(CP_N$CPN)
    

    while (length(ku_stack)!=0) { # after k_list size
      
      CP_N_k1st = CP_N[min(which(CP_N$K == ku_stack[1])):max(which(CP_N$K == ku_stack[1])),]
      CP_N_k1st_tm <- merge(CP_N_k1st, Coords_df ,by="sample" )
      colnames(CP_N_k1st_tm)[2] = "CP_N"
      Title1 = as.character(paste("k = ", as.character(ku_stack[1])))
      p1_seq<- plot_ly(CP_N_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                       marker=list( size=10 , opacity=1), color = ~CP_N, text = ~paste('Sample: ', sample))%>%
        layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      #p1_seq<- colorbar(p1_seq, limits = c(binf , bsup) )
      
      CP_N_k2nd = CP_N[min(which(CP_N$K == ku_stack[2])):max(which(CP_N$K == ku_stack[2])),]
      CP_N_k2nd_tm <- merge(CP_N_k2nd, Coords_df ,by="sample" )
      colnames(CP_N_k2nd_tm)[2] = "CP_N"
      Title2 = as.character(paste("k = ", as.character(ku_stack[2])))
      p2_seq<- plot_ly(CP_N_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                       marker=list( size=10 , opacity=1), color = ~CP_N, text = ~paste('Sample: ', sample))%>%
        layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      #p2_seq<- colorbar(p2_seq, limits = c(binf, bsup) )
      
      CP_N_k3rd = CP_N[min(which(CP_N$K == ku_stack[3])):max(which(CP_N$K == ku_stack[3])),]
      CP_N_k3rd_tm <- merge(CP_N_k3rd, Coords_df ,by="sample" )
      colnames(CP_N_k3rd_tm)[2] = "CP_N"
      Title3 = as.character(paste("k = ", as.character(ku_stack[3])))
      p3_seq<- plot_ly(CP_N_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                       marker=list( size=10 , opacity=1), color = ~CP_N, text = ~paste('Sample: ', sample))%>%
        layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      # p3_seq <- colorbar(p3_seq, limits = c(binf, bsup ) )
      
      
      CP_N_k4th = CP_N[min(which(CP_N$K == ku_stack[4])):max(which(CP_N$K == ku_stack[4])),]
      CP_N_k4th_tm <- merge(CP_N_k4th, Coords_df ,by="sample" )
      colnames(CP_N_k4th_tm)[2] = "CP_N"
      Title4 = as.character(paste("k = ", as.character(ku_stack[4])))
      p4_seq<- plot_ly(CP_N_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers", 
                       marker=list( size=10 , opacity=1), color = ~CP_N, text = ~paste('Sample: ', sample))%>%
        layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE)
      #p4_seq<- colorbar(p4_seq, limits = c(binf,bsup) )
      
      ku_stack <- ku_stack[-c(1:4)]
      p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)%>%
        layout(title = mytitle,  margin = 0.04)
      
      print(p_seq)
      
    }
    
  }
  
}



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





CP_mean_by_k  <-function (list_CP_diff , Name){
  par(mfrow=c(1,2))
  c= 1 
  c_data_frame = as.data.frame(list_CP_diff[1])
  CP_diff_mean_df = data.frame("k" =unique(c_data_frame$K ))
  for (i in 1:length(list_CP_diff)){
    print(Name[i])
    c_data_frame = as.data.frame(list_CP_diff[i])
    print("here1")
    diff_CP2_CPN = abs(c_data_frame$CP2 - c_data_frame$CPN )
    c_data_frame = cbind(c_data_frame ,diff_CP2_CPN )
    print("here2")
    colnames(c_data_frame)[dim(c_data_frame)[2]] <- "Abs_diff"
    Mean_by_k =tapply(c_data_frame$Abs_diff ,c_data_frame$K, mean)
    CP_diff_mean_df <- cbind(CP_diff_mean_df ,Mean_by_k )
    colnames(CP_diff_mean_df)[dim(CP_diff_mean_df)[2]] <- Name[i]
    
  }
  print(dim(CP_diff_mean_df))
  Col = c()
  C=1
  for (i in 2:dim(CP_diff_mean_df)[2]){
    
    Col = c(Col,C)
    if (i ==2 ){
      plot(CP_diff_mean_df$k, CP_diff_mean_df[,i] , col=C , type="l"  , xlab = "k", ylab="mean(|CP2-CPN|)", main= 'centrality preservarion' )
      C =C+1
      Col = c(Col,C)
    }
    else{
      lines(CP_diff_mean_df$k, CP_diff_mean_df[,i] , col=C )
      C =C+1
      Col = c(Col,C)
    }
    
  }
  print(unique(Col))
  plot.new()
  legend("bottomleft",  legend = Name, col = c(unique(Col)), pch = 15)
  
  return
  
}




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
      plot(Set_diff_mean_df$k, Set_diff_mean_df[,i] , col=C  , xlab="k",ylab = "Mean Set Difference ", main="Means Set Difference view by level k ", type ='l' )
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
      plot(Seq_diff_mean_df$k, Seq_diff_mean_df[,i] , col=C  , xlab = "k" , ylab = "Mean Seq Diff", main = "Mean of Sequences difference view by level k", type ='l' )
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
  legend("bottomleft",  legend = Name,  col = c(unique(Col)), pch = 15)
}


