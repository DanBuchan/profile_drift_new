import numpy as np

with open('/home/dbuchan/Data/distance_mat/rand_rep_distance_matrix.npy', 'rb') as f:
    dist_matrix = np.load(f)
    dom_list = np.load(f)

num_cells = dist_matrix.shape[0] * dist_matrix.shape[1]
print(num_cells)
vals_less_than = (dist_matrix <= 0.91).sum() 
print(vals_less_than)

percentage_of_close =(vals_less_than/num_cells)
print(percentage_of_close)