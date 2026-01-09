import numpy as np 
import sys
import sim
import os
from joblib import Parallel, delayed

#################### Invert ligand maps #######################
def inv_maps(lgnd_maps):

    lgnd_inv_maps_lst = []
    for i in range(lgnd_maps.shape[0]):
        matrix = lgnd_maps[i]
        inverted_matrix = matrix[::-1,:]
        lgnd_inv_maps_lst.append(inverted_matrix)

    lgnd_inv_maps = np.array(lgnd_inv_maps_lst)
    
    return lgnd_inv_maps                  
################################################################

### Find indices of maps next to selected surface dots #########
def find_map_indcs(atm_pair, coords_arr):

    map_indx_lst = []
    for atm_indx_i in range(len(atm_pair)):
        dot_i = atm_pair[atm_indx_i]
        dist_map_to_dot_i = np.sqrt((coords_arr[:, 0] - dot_i[0])**2 + (coords_arr[:, 1] - dot_i[1])**2 + (coords_arr[:, 2] - dot_i[2])**2)
        map_indx_i = np.argmin(dist_map_to_dot_i)
        map_indx_lst.append(map_indx_i)

    return map_indx_lst
################################################################

####################### Map-map comparison #####################
def srch_diff_v4(qry_maps_mat, sbjct_maps_mat):

        r_map_ij_lst = list()

        for q_map_i in range(qry_maps_mat.shape[0]):
                tmp_lgnd_map = qry_maps_mat[q_map_i, :]
    
                for r_map_j in range(sbjct_maps_mat.shape[0]):
                        tmp_rcptr_map = sbjct_maps_mat[r_map_j,:]
                        diff = sim.ssim(tmp_lgnd_map.astype(np.float32), \
                            tmp_rcptr_map.astype(np.float32), \
                            data_range=1, \
                            full=False)
                        r_map_ij_lst.append([-1*diff, q_map_i, r_map_j])

                r_map_ij_arr = np.array(r_map_ij_lst)
        
        return r_map_ij_arr
################################################################

################################################################

def map_analysis(sbjct_path, qry_path, rf_cutoff, atm_pair):

    in_lgnd = np.load(qry_path, allow_pickle=True)
    lgnd_coords = in_lgnd["coords"]
    lgnd_maps = in_lgnd["maps"]
    lgnd_nrmls = in_lgnd["nrmls"]

    map_indx_lst = find_map_indcs(atm_pair, lgnd_coords)
        
    lgnd_maps_sele = lgnd_maps[map_indx_lst]
    lgnd_coords_sele = lgnd_coords[map_indx_lst]
    lgnd_nrmls_sele = lgnd_nrmls[map_indx_lst]
    lgnd_inv_maps = inv_maps(lgnd_maps_sele)

    intr_qry_dist = np.linalg.norm(lgnd_coords_sele[0] - lgnd_coords_sele[1])
    dist_range = np.array([(intr_qry_dist-dist_tol),(intr_qry_dist+dist_tol)])

    npz_name = sbjct_path.split("/")[-1]
    in_rcpt = np.load(sbjct_path, allow_pickle=True)
    rcptr_coords = in_rcpt["coords"]
    rcptr_maps = in_rcpt["maps"]
    rcptr_nrmls = in_rcpt["nrmls"]

    r_ij = srch_diff_v4(lgnd_inv_maps, rcptr_maps)

    r_ij_0 = np.array([e for e in r_ij if int(e[1]) == 0])
    r_ij_1 = np.array([e for e in r_ij if int(e[1]) == 1])
    r_ij_0_sort = r_ij_0[r_ij_0[:,0].argsort()]
    r_ij_1_sort = r_ij_1[r_ij_1[:,0].argsort()]

    ###############################################################
    # select X % receptor maps with the lowest R-factor to each query
    ###############################################################
    n_top_hits = int(r_ij_0_sort.shape[0]*top_pct/100)

    top_qry_0 = r_ij_0_sort[:n_top_hits,:]
    top_qry_1 = r_ij_1_sort[:n_top_hits,:]

    top_qry_0_dict = {row[2]: row[0] for row in top_qry_0}
    top_qry_1_dict = {row[2]: row[0] for row in top_qry_1}

    ###############################################################
    # get pairwise dists for top hits from 2 queries
    ###############################################################
    top_qry_0_coords_lst = []
    for i in top_qry_0[:,2]:
        top_qry_0_coords_lst.append(rcptr_coords[int(i)])
    top_qry_0_coords = np.array(top_qry_0_coords_lst)

    top_qry_1_coords_lst = []
    for i in top_qry_1[:,2]:
        top_qry_1_coords_lst.append(rcptr_coords[int(i)])
    top_qry_1_coords = np.array(top_qry_1_coords_lst)

    dist_0_1 = np.sqrt(np.sum((top_qry_0_coords[:, np.newaxis] - top_qry_1_coords)**2, axis=2)) #pairwise distances

    dist_ij = np.empty((dist_0_1.size, 3))
    dist_ij[:, 0] = dist_0_1.flatten() 
    dist_ij[:, 1] = np.repeat(top_qry_0[:,2], len(top_qry_1))
    dist_ij[:, 2] = np.tile(top_qry_1[:,2], len(top_qry_0))

    dist_filt = (dist_ij[:, 0] >= dist_range[0]) & (dist_ij[:, 0] <= dist_range[-1])
    filt_map_pairs = dist_ij[dist_filt]

    if filt_map_pairs.size == 0:
        return np.array([])

    rf_avrg = []
    for a in range(len(filt_map_pairs)):
        rf_avrg.append((top_qry_0_dict[filt_map_pairs[a,1]]+top_qry_1_dict[filt_map_pairs[a,2]])/2)
        
    filt_map_pairs_rf = np.column_stack((filt_map_pairs, rf_avrg))
    # col0: distance between two subject maps;
    # col1: indx of subject map 1; 
    # col2: indx of subject map 2;
    # col3: average R-factor for two query patches

    filt_map_pairs_rf = filt_map_pairs_rf[filt_map_pairs_rf[:,3]<=rf_cutoff]

    if filt_map_pairs_rf.size == 0:
        return np.array([]) 
    rcptr_indcs_0 = filt_map_pairs_rf[:, 1].astype(int)
    rcptr_indcs_1 = filt_map_pairs_rf[:, 2].astype(int)

    coords_0 = rcptr_coords[rcptr_indcs_0] 
    coords_1 = rcptr_coords[rcptr_indcs_1]

    nrmls_0 = rcptr_nrmls[rcptr_indcs_0]
    nrmls_1 = rcptr_nrmls[rcptr_indcs_1]

    q_coords_0 = lgnd_coords[map_indx_lst[0]][np.newaxis, :] 
    q_coords_1 = lgnd_coords[map_indx_lst[1]][np.newaxis, :] 
    q_nrmls_0 = lgnd_nrmls[map_indx_lst[0]][np.newaxis, :]  
    q_nrmls_1 = lgnd_nrmls[map_indx_lst[1]][np.newaxis, :] 

    if coords_0.shape[0] > 1:
        q_coords_0 = np.tile(q_coords_0, (coords_0.shape[0], 1)) 
        q_coords_1 = np.tile(q_coords_1, (coords_1.shape[0], 1))
        q_nrmls_0 = np.tile(q_nrmls_0, (nrmls_0.shape[0], 1))
        q_nrmls_1 = np.tile(q_nrmls_1, (nrmls_1.shape[0], 1))

    filt_map_pairs_rf_coords = np.concatenate((filt_map_pairs_rf, \
                                               coords_0, \
                                               nrmls_0, \
                                               coords_1, \
                                               nrmls_1, \
                                               q_coords_0, \
                                               q_nrmls_0, \
                                               q_coords_1, \
                                               q_nrmls_1), axis=1)
    # col0: distance between two subject maps;
    # col1: indx of subject map 1; 
    # col2: indx of subject map 2;
    # col3: average R-factor for two query patches;
    # col4, col5, col6: coords of subject map 1;
    # col7, col8, col9: normal coords for subject map 1;
    # col10, col11, col12: coords of subject map 2;
    # col13, col14, col15: normal coords for subject map 2;
    # col16, col17, col18: coords of query map 1;
    # col19, col20, col21: normal coords for query map 1;
    # col22, col23, col24: coords of query map 2;
    # col25, col26, col27: normal coords for query map 2;

    filt_map_pairs_rf_coords_pdb = np.hstack((filt_map_pairs_rf_coords, np.full((filt_map_pairs_rf_coords.shape[0], 1), '{}'.format(npz_name), dtype=object)))

    print(npz_name)
    return filt_map_pairs_rf_coords_pdb

################################################################   

if __name__ == "__main__":

    try:
        npz_dir = sys.argv[1]

        sbjct_npzs = [file for file in os.listdir(npz_dir) if file.endswith(".npz")]

        qry_npz = sys.argv[2]

        rf_cutoff = float(sys.argv[3])
    except:
        raise Exception("""
    input:
        npz_dir := path where sbjct npzs are stored
        qry_npz := query npz file
        rf_cutoff := r-factor cutoff; min => -1, max => 1
            """)

    n_cores = 10

    atm0 = np.array([12.619, 18.25, -15.497])
    atm1 = np.array([9.304, 19.759, -4.274])
    atm2 = np.array([-2.767, 20.462, -7.706])

    atm_prs = [(atm0, atm1), \
               (atm1, atm2), \
               (atm2, atm0)]

    top_pct = 5 # fraction (in percents) of top hits for every query map
    dist_tol = .01 # distance tolerance between two subject maps

    rslts_lst = Parallel(n_jobs=n_cores)(
        delayed(map_analysis)(os.path.join(npz_dir, sbjct_npz), os.path.join(npz_dir, qry_npz), rf_cutoff, atm_pr)
        for sbjct_npz in sbjct_npzs for atm_pr in atm_prs)

    rslts_arr = np.vstack([rslt for rslt in rslts_lst if rslt.any()])
    np.save("/results/srch_rslts", rslts_arr)