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

	return dist


def jaccard_similarity(dist1 , dist2 , k):
	if dist1.shape[0] == dist2.shape[0] :
		n = dist1.shape[0]
		for i in range(n):
			N1_df = pd.DataFrame(dist1.iloc[i,:], columns=dist1.iloc[ : , 0])
			N2_df = pd.DataFrame(dist2.iloc[i,:], columns=dist2.iloc[ : , 0])
			N1_df.sort(axis = 0 , ascending=[1, 0])
			print(N1_df)
	else :
		"Dim error"


if __name__ == '__main__': 
	PCA_coords_df = pd.read_csv("fig6A_PCA_coords.tab", sep="\t")
	MOFA_coords_df = pd.read_csv("fig13A_MOFA_Coords.tab", sep="\t")
	TM_coords_df = pd.read_csv("TM_coords_expr.tab", sep="\t")
	sh = PCA_coords_df.shape
	print(distance_matrix(PCA_coords_df).iloc[1:10, 1:10])
	#print(PCA_coords_df.head())