from __future__ import division
import pandas as pd
import numpy as np
import math 


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
			print(dim(N1_df))
			kNeighbors_N1 =list(N1_df.index)[:k]
			kNeighbors_N2 =list(N2_df.index)[:k]
			linter = len(set(kNeighbors_N1).intersection(set(kNeighbors_N2)))
			lunion = len(set(kNeighbors_N1)) + len(set(kNeighbors_N2)) - linter
			c_JS = 1- (linter / lunion)
			Jsim.iloc[i,0] = c_JS
		return Jsim
	else :
		"Dim error"
	return 0


def sequence_difference(dist1 , dist2 , k):
	Jsim = pd.DataFrame(index=list(dist1.columns.values))
	Jsim["J_sim"] = range(len(list(dist1.columns.values)))
	if dist1.shape[0] == dist2.shape[0] :
		n = dist1.shape[0]
		for i in range(n):
			N1_df = pd.DataFrame(dist1.iloc[i,:], index=list(dist1.columns.values))
			N2_df = pd.DataFrame(dist2.iloc[i,:], index=list(dist2.columns.values))
			N1_df = N1_df.sort_values(by=list(N1_df.columns.values), ascending= True)
			N2_df = N2_df.sort_values(by=list(N2_df.columns.values), ascending= True)
			#pd.merge(left, right, how='inner', on=None)
	else :
		"Dim error"
	return 0


if __name__ == '__main__': 
	PCA_coords_df = pd.read_csv("fig6A_PCA_coords.tab", sep="\t")
	print("PCA_coords_df", PCA_coords_df.shape)
	MOFA_coords_df = pd.read_csv("fig13A_MOFA_Coords.tab", sep="\t")
	print("MOFA_coords_df", MOFA_coords_df.shape)
	TM_coords_df = pd.read_csv("TM_coords_expr.tab", sep="\t")
	print("TM_coords_df", TM_coords_df.shape)
	sh = PCA_coords_df.shape
	D1 = distance_matrix(PCA_coords_df)
	D2 = distance_matrix(TM_coords_df)
	print(set_difference(D1 , D2 , 10))
	#print(distance_matrix(PCA_coords_df).iloc[1:10, 1:10])
	#print(PCA_coords_df.head())