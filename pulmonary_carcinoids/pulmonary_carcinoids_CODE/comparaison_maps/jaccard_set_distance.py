#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
import pandas as pd
import numpy as np
import math 
import sys
import os

def creation_fichier(nom_fichier) :
    """Cette fonction permet de créer les fichiers textes contenant les résultats. Cette fonction prévient la redondance des fichiers,
    ainsi le nom de chaque fichier est unique. """

    nom_fichier=nom_fichier.replace("\n","")
    fichier_existe=True # Variable permettant de verifier que le fichier qu'on va creer n'en ecrase pas un preexistant.
    numero_fichier=0
    while fichier_existe: # Tant que le fichier "nom_fichier.png" existe le nom change.
        try:
            sortie=open(nom_fichier+"(%i).txt" % numero_fichier,'r') # Test si le fichier "nom_fichier.py" existe.
        except IOError:
        	fichier_existe=False
        else:
            sortie.close()
            numero_fichier+=1
            nom_fichier=nom_fichier.replace("(%i)" % (numero_fichier-1),"(%i)" % numero_fichier)  
    return(nom_fichier, numero_fichier)


def distance_matrix(data_coords) :
	# data_coords is a pandas data frame with sample , x , y
	n = data_coords.shape[0]
	d =  np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			c_x1 = data_coords.iloc[i , 1]
			c_y1 = data_coords.iloc[i , 2]
			c_x2 = data_coords.iloc[j , 1]
			c_y2 = data_coords.iloc[j , 2]
			d[i,j] = math.sqrt((c_x1-c_x2)**2  + (c_y1 - c_y2)**2)
	dist = pd.DataFrame(d, columns=data_coords.iloc[ : , 0], index=data_coords.iloc[  : , 0])
	dist = dist.reindex(sorted(dist.columns), axis=1)

	return dist


def set_difference(dist1 , dist2 , k):
	Jsim = pd.DataFrame(index=list(dist1.columns.values))
	Jsim["J_sim"] = range(len(list(dist1.columns.values)))
	if dist1.shape[0] == dist2.shape[0] :
		n = dist1.shape[0]
		for i in range(n):
			N1_df = pd.DataFrame(dist1.iloc[i,:], index=list(dist1.columns.values))
			N2_df = pd.DataFrame(dist2.iloc[i,:], index=list(dist2.columns.values))
			N1_df = N1_df.sort_values(by=list(N1_df.columns.values), ascending= True)
			N2_df = N2_df.sort_values(by=list(N2_df.columns.values), ascending= True)
			kNeighbors_N1 =list(N1_df.index)[:k]
			kNeighbors_N2 =list(N2_df.index)[:k]
			linter = len(set(kNeighbors_N1).intersection(set(kNeighbors_N2)))
			lunion = len(set(kNeighbors_N1)) + len(set(kNeighbors_N2)) - linter
			c_JS = 1- (linter / lunion)
			Jsim.iloc[i,0] = c_JS
		return Jsim
	else :
		print("Dim error")
	return 0


def centrality_preservation(dist1 , dist2 , K , filename):
	"""In this function the '2' refers to a the projection and 'N' to the original dimension  """

	name_centrality_preservation =  creation_fichier(filename)[0]      
	centrality_preservation_file = open(name_centrality_preservation  ,'a')
	header = "sample" + '\t' + "CP2"+ '\t' + "CPN" + '\t'  + 'K' + '\n'
	centrality_preservation_file.write(header)
	for k in K :
		CP2 = pd.DataFrame(index=list(dist1.columns.values))
		CP2["cp"] = range(len(list(dist1.columns.values)))
		CPN = pd.DataFrame(index=list(dist2.columns.values))
		CPN["cp"] = range(len(list(dist2.columns.values)))
		if dist1.shape[0] == dist2.shape[0] :
			n = dist1.shape[0]
			
			# Data frame of neighboors
			N2_df = pd.DataFrame(columns = range(k) ,index=list(dist1.columns.values))
			NN_df = pd.DataFrame(columns = range(k) ,index=list(dist1.columns.values))
			c =0
			for i in range(N2_df.shape[0]):
				N2_dist_df = pd.DataFrame(dist1.iloc[i,:], index=list(dist1.columns.values))
				NN_dist_df = pd.DataFrame(dist2.iloc[i,:], index=list(dist2.columns.values))
				N2_dist_df  = N2_dist_df.sort_values(by=list(N2_dist_df.columns.values), ascending= True)
				NN_dist_df = NN_dist_df.sort_values(by=list(NN_dist_df.columns.values), ascending= True)
				KNeighbors_N2 =list(N2_dist_df.index)[:k]
				KNeighbors_NN =list(NN_dist_df.index)[:k]
				N2_df.iloc[i,] = KNeighbors_N2
				NN_df.iloc[i,] = KNeighbors_NN
		
			for i in range(N2_df.shape[0]):
				# Current vertex
				c_vertex = N2_df.index.values[i]
				CP2_j = 0
				CPN_j = 0
				for j in range(N2_df.shape[0]):
					neighbors2_j = list( N2_df.iloc[j,])
					neighborsN_j = list(NN_df.iloc[j,])
					if c_vertex in neighbors2_j :
						CP2_j += k - neighbors2_j.index(c_vertex)	
					else :
						CP2_j += 0
					
					if c_vertex in neighborsN_j :
						CPN_j += k - neighborsN_j.index(c_vertex)	
					else :
						CPN_j += 0
				CP2.iloc[i,0]	= CP2_j
				CPN.iloc[i,0]	= CPN_j

			for l in range(CP2.shape[0]):
				#print('level k', k)
				line = str(CP2.index.values[l]) + '\t' + str(CP2.iloc[l,0]) + '\t' + str(CPN.iloc[l,0]) +'\t'+  str(k) + '\n'
				#print('Line ', line)
				centrality_preservation_file.write(line)
	
	return CP2 , CPN

	


def sequence_difference(dist1 , dist2 , k):
	Jsim = pd.DataFrame(index=list(dist1.columns.values))
	seq_diff = []
	if dist1.shape[0] == dist2.shape[0] :
		n = dist1.shape[0]
		for i in range(n):
			c_vertex =  dist1.index.values[i]
			N1_df = pd.DataFrame(dist1.iloc[i,:], index=list(dist1.columns.values))
			N2_df = pd.DataFrame(dist2.iloc[i,:], index=list(dist2.columns.values))
			N1_df = N1_df.sort_values(by=list(N1_df.columns.values), ascending= True)
			N1_df["rank_x"] = range(N1_df.shape[0])
			N1_df  = N1_df.iloc[:k,]
			N2_df = N2_df.sort_values(by=list(N2_df.columns.values), ascending= True)
			N2_df["rank_y"] = range(N2_df.shape[0])
			N2_df  = N2_df.iloc[:k,]
			N = pd.merge(N1_df, N2_df ,  how='inner', left_index=True, right_index=True)
			#print(N.head)
			#print(N.shape)
			s1 = 0
			s2 = 0
			for i in range(N.shape[0]):
				s1 += (k - N["rank_x"][i]) * abs(N["rank_x"][i] - N["rank_y"][i])   
				s2 += (k - N["rank_y"][i]) * abs(N["rank_x"][i] - N["rank_y"][i]) 
			S = 0.5 * s1 + 0.5 * s2
		
		 	seq_diff.append(S)
	else :
		"Dim error"

	Jsim['seq_diff'] = seq_diff

	return Jsim


def knn_high_dim(data_feature, filename ,k) :
	name_dist_file =  creation_fichier(filename)[0]      
	dist_file = open(name_dist_file  ,'a')
	header_dist = str(list(data_feature.index.values )) + '\n'
	dist_file.write(header_dist)
	# Euclidian distance in n dim 
	n = data_feature.shape[0] # sammple
	m = data_feature.shape[1] # gene
	d =  np.zeros((n, n))
	for i in range(n):
		for j in range(n):
			if i == j:
				d[i,j] = 0
			else :
				if d[i,j] != 0 :
					d[j,i] = d[i,j]
				else :
					s = 0
					for k in range(m):
						s+= (data_feature.iloc[i,k] - data_feature.iloc[j,k])**2
					d[i,j] = math.sqrt(s)
		dist_file.write(str(d[i,]) + '\n')
	dist = pd.DataFrame(d, columns=data_feature.index.values, index=data_feature.index.values )
	#dist = dist.reindex(sorted(dist.columns), axis=1)
	print(dist.head())
	return dist	

def main(df1, df2, k , filename_set_diff, filename_seq_diff):
	if k <= df1.shape[0] and df1.shape == df2.shape :
		dist1 = distance_matrix(df1)
		dist2 = distance_matrix(df2)
		name_set_diff_file =  creation_fichier(filename_set_diff)[0]      
		set_diff_file = open(name_set_diff_file  ,'a')
		name_sequence_diff_file = creation_fichier(filename_seq_diff)[0]  #
		seq_diff_file = open(name_sequence_diff_file  ,'a')
		header_set = "sample" + '\t' + "set_diff" + '\t'  + 'k' + '\n'
		set_diff_file.write(header_set)
		header_seq = "sample" + '\t' + "seq_diff" + '\t'  + 'k' + '\n'
		seq_diff_file.write(header_seq)

		print("A file named set_diff.txt have been created")
		for i in range(1, k+1):
			c_set_diff = set_difference(dist1 , dist2 , i)
			c_seq_diff = sequence_difference(dist1 , dist2 , i)
			for l in range(c_set_diff.shape[0]) : 
				line_c_set_diff =str(c_set_diff.index[l]) + '\t' + str(c_set_diff.iloc[l,0]) + '\t'  +str(i) + '\n'
				line_c_seq_diff =str(c_seq_diff.index[l]) + '\t' + str(c_seq_diff.iloc[l,0]) + '\t'  +str(i) + '\n'
				set_diff_file.write(line_c_set_diff)
				seq_diff_file.write(line_c_seq_diff)




if __name__ == '__main__': 
	#PCA_coords_df = pd.read_csv("Meso_pca_coords.tab", sep="\t")
	#TM_coords_df= pd.read_csv("Meso_tm_coords_v2.tab", sep="\t")
	#d1 =  distance_matrix(PCA_coords_df)
	#d2 = distance_matrix(TM_coords_df)
	#seq_diff = sequence_difference(d1,d2, 120)
	#centrality_preservation(d1,d2, [60 , 170,  260 ,280], 'CP_MesosomicsV2.txt')

	#main(PCA_coords_df, TM_coords_df, PCA_coords_df.shape[0] , "set_diff_meso" , "seq_diff_meso" ) #PCA_coords_df.shape[0]
	#D1_PCA = distance_matrix(PCA_coords_df)
	#D2_TM = distance_matrix(TM_coords_df)
	#set_diff_10 = set_difference(D1_PCA , D2_TM , 10)
	#print(set_diff_10)


	#  NN 
	Feature_data_df = pd.read_csv("feature_data_with_lv_2.tsv", sep="\t")
	Feature_data_df = Feature_data_df.transpose() 
	Feature_data_df.columns = Feature_data_df.iloc[0,]
	Feature_data_df= Feature_data_df.drop(Feature_data_df.index[0])
	Feature_data_df.to_csv('Feature_data_t.csv')



	#Feature_data_df= Feature_data_df.iloc[1:20,1:20]
	#knn_high_dim(Feature_data_df, 'dist_meso_data.txt',5)
	#print(Feature_data_df.head())
	#print(Feature_data_df.columns.values)
	#print(Feature_data_df.index.values)



