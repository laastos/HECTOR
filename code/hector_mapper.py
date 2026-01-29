#Copyright (C) 2018 Mohammad ElGamacy / Max Planck Society - All Rights Reserved

import re
import sys
import glob
import time
import argparse
import numpy as np

###############################################################
def parse_ply(ply_fn, vrtx_cols=6, fc_cols=7):
  """A function to parse ply files
  This function is compatible with PLY files generated from EDTSurf or MeshLab

  parse_ply(ply_fn, vrtx_cols=6, fc_cols=7)

  Args:
  ply_fn (str): input PLY file name
  vrtx_cols (int): number of vertex information columns (optional)
  fc_cols (int): number of face information columns (optional)

  Returns:
  vrtx_arr (numpy.ndarray): an array of floats with shape
  N_vertices*vrtx_cols)
  fc_arr (numpy.ndarray): an array of ints with shape (N_faces*fc_cols) shape
  """
  ply_fh = open(ply_fn, 'r')
  ply_str = ply_fh.read()
  ply_fh.close()

  splt_str = ply_str.split("end_header\n")
  vrtx_cnt = int(re.search("element\s+vertex\s+(\d+)", splt_str[0]).group(1))
  fc_cnt = int(re.search("element\s+face\s+(\d+)", splt_str[0]).group(1))
  # split raw text and extract vertex and face information as lists of strings
  ntry_strlst = splt_str[1].splitlines()
  vrtx_str = "\n".join(ntry_strlst[:vrtx_cnt])
  fc_str = "\n".join(ntry_strlst[vrtx_cnt:])

  # extract all vertex columns
  vrtx_arr = np.fromstring(vrtx_str, dtype=np.float64, sep=" ").reshape(-1, vrtx_cols)
  # extract all face columns
  fc_arr = np.array([[int(float(idx)) for idx in re.split("\s+", idcs.strip())] for idcs in ntry_strlst[vrtx_cnt:]])

  return vrtx_arr, fc_arr
###############################################################

###############################################################
def normalise(in_arr):
  """ Vector-normaliser 
  
  normalise(in_arr)

  Args: 
    in_arr (numpy.ndarray): a matrix of surface triangles surface normal vectors
    
  Returns
    norm_arr (numpy.ndarray): a matrix of surface triangles unit surface normal vectors
  """
  
  norm_arr = np.copy(in_arr)
  norms = np.linalg.norm(in_arr, axis=1)

  norm_arr[:, 0] /= norms
  norm_arr[:, 1] /= norms
  norm_arr[:, 2] /= norms

  return norm_arr 
###############################################################
          
###############################################################
# Spin map generator 
# VECTORISED IMPLEMENTATION
###############################################################
def spn_map_srfc(srfc_vrtx_arr, srfc_nrml_arr, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq):

  """The core fingerprinting function
  This function maps given vertices depending on the direction of the input surface normals. 

  Args: 
    srfc_vrtx_arr (numpy.ndarray): a matrix of input vertices position vectors
    srfc_nrml_arr (numpy.ndarray): a matrix of input vertices surface normal vectors
    spprt_dstnc_x (float): support distance alpha; radius of the cylinderical basis in Angstroms
    spprt_dstnc_y (float): support distance beta; height of the cylinderical basis in Angstroms 
    bin_sz (float): size of the bin; i.e. resolution of the map in Angstroms
    dot_skp_frq (int): dot-skipping frequency; the sparsity of sampling for references vertices 
    dot_map_frq (int): dot-skipping frequency; the sparsity of sampling for mapped vertices 
    
  Returns
    vrtx_spn_maps (list(numpy.ndarray)): A list of the final fingerprints 
  """

  # S_O: R^3 -> R^2 
  # S_O(x) -> (alpha, beta) = ( (||x-p||^2 - (n.(x-p))^2)^0.5 , (n.(x-p)) )
  #
  # x := ovrtx; accumulated vertex coords
  # p := vrtx; point-of-view (i.e. reference) vertex coords
  # n := nrml; point-of-view (reference) unit surface normal
  # oriented point O := surface mesh vertex with position vector p and unit surface normal n 

  print("generating local maps ...")
  # total number of bins along one dimension
  # int() floors and does not round; hence the "+ 2" to pad
  # tot_bins:= matrix number of rows and columns
  tot_bins_x = int(spprt_dstnc_x/bin_sz) + 2
  tot_bins_y = int(spprt_dstnc_y/bin_sz) + 2 
  half_spprt_dstnc = spprt_dstnc_y/2
  vrtx_spn_maps = list()

  c=0
  for vrtx, n in zip(srfc_vrtx_arr[::dot_skp_frq,:], srfc_nrml_arr[::dot_skp_frq,:]): 
  # loop over query vertices and their normals
    c+=1; 
    if c % 1000 == 0: print(c) 
    spn_img = np.zeros((tot_bins_y,tot_bins_x), dtype=np.float64) 
    # create an empty spin image as a tot_bins_y*tot_bins_x matrix

    # xp = ovrtx - vrtx # x - p 
    # beta = np.dot(n, xp) 
    # beta can be positive or negative
    # dst_sq = (np.linalg.norm(xp))**2 - (beta**2)
    # if abs(beta) <= half_spprt_dstnc: 
    # limit sampling sphere radius by half of the support distance
    # alpha =  np.sqrt( (np.linalg.norm(xp))**2 - (beta**2) )
    # alpha is always positive
    # if alpha <= spprt_dstnc: 

    # Array of "x - p" 
    xp_arr = srfc_vrtx_arr[::dot_map_frq,:] - vrtx
    # Array of all beta's with respect to reference vertex
    beta_arr = np.dot(n, xp_arr.T)

    lmt_idx_arr = np.where(np.abs(beta_arr) <= half_spprt_dstnc)[0]
    lmt_xp_arr = xp_arr[lmt_idx_arr]
    lmt_beta_arr = beta_arr[lmt_idx_arr]
    lmt_xp_norm_arr = np.linalg.norm(lmt_xp_arr, axis=1)
    alpha_arr = np.sqrt((lmt_xp_norm_arr*lmt_xp_norm_arr) - (lmt_beta_arr*lmt_beta_arr))
    lmt_lmt_idx_arr = np.where(alpha_arr <= spprt_dstnc_x)
    lmt_lmt_alpha_arr = alpha_arr[lmt_lmt_idx_arr]
    lmt_lmt_beta_arr = lmt_beta_arr[lmt_lmt_idx_arr]

    #limit sampling sphere radius by the support distance
    # x = (half_spprt_dstnc - beta)/bin_sz
    # y = alpha/bin_sz
    # i = int(x) #flooring down to the i coordinate 
    # j = int(y) #flooring down to the j coordinate

    x_arr = (half_spprt_dstnc - lmt_lmt_beta_arr)/bin_sz
    y_arr = lmt_lmt_alpha_arr/bin_sz
    i_arr = x_arr.astype(int)
    j_arr = y_arr.astype(int)  

    for i, j in zip(i_arr, j_arr): 

      spn_img[i][j] += 1

    vrtx_spn_maps.append(spn_img[1:-1,:-2])
    # print("spn_img shape: ", vrtx_spn_maps[-1].shape)
##############################################################

  return vrtx_spn_maps

###############################################################

################################################################
############ Sigmoid radial fading of the maps #################

def generate_circle_matrix_3d(input_3d_matrix, radius_fade, slope_fade):
    sizes = input_3d_matrix.shape[1:]
    filtermap = np.zeros(sizes, dtype=float)
    bestypoint = sizes[0] // 2
    for x in range(sizes[0]):
        for y in range(sizes[1]):
            distance = np.sqrt((bestypoint - x) ** 2 + (y) ** 2)
            filtermap[x, y] = max(0, 1 - logistic(distance - radius_fade, k=slope_fade))
    filtermap = 1 - filtermap
    
    return input_3d_matrix * filtermap
################################################################

################################################################
def logistic(x, k=10):
    return 1 - 1 / (1 + np.exp(-k * x))
################################################################

################# min-max map normalization ####################
def minmax_norm(in_arr_maps):
    maps_norm_lst = []
    for map_i in range(in_arr_maps.shape[0]):
        map_i_norm = (in_arr_maps[map_i] - np.min(in_arr_maps[map_i]))/(np.max(in_arr_maps[map_i]) - np.min(in_arr_maps[map_i]))
        maps_norm_lst.append(map_i_norm)

    maps_norm_arr = np.array(maps_norm_lst)

    return maps_norm_arr
################################################################

###############################################################
def hector_mapper_vctrsd(ply_fn, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq, sign_flg, infltn, radius_fade, slope_fade):

  """The mapping and PLY-preprocessing wrapper function

  Args: 
    ply_fn (str): input PLY file name 
    spprt_dstnc_x (float): support distance alpha; radius of the cylinderical basis in Angstroms
    spprt_dstnc_y (float): support distance beta; height of the cylinderical basis in Angstroms 
    bin_sz (float): size of the bin; i.e. resolution of the map in Angstroms
    dot_skp_frq (int): dot-skipping frequency; the sparsity of sampling for references vertices 
    dot_map_frq (int): dot-skipping frequency; the sparsity of sampling for mapped vertices 
    sign_flg (str): the mapping polarity; i.e. ligand or receptor mapping; only ("lgnd","rcpt") flags are valid
    infltn (float): surface inflation in Angstroms
    radius_fade (float): fading radius in bins
    slope_fade (float): fading slope 
  """

  t_strt = time.time()

  sign = {"lgnd": -1, "rcpt": 1}

  vrtx_arr, fc_arr = parse_ply(ply_fn, vrtx_cols=6, fc_cols=7)
  # extract only the coordinates
  ply_vrtx_arr = vrtx_arr[:, :3]
  # extract only the three-point cliques
  ply_fc_arr = fc_arr[:, 1:4]

  #Create an indexed view into the vertex array using the array of three indices for triangles
  tris = ply_vrtx_arr[ply_fc_arr]
  
  #Calculate the normal for all the triangles, by taking the cross product of the vectors v1-v0, and v2-v0 in each triangle             
  tris_nrmls = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
  
  # n is now an array of normals per triangle. The length of each normal is dependent the vertices, 
  # these must be normalised, so weigh each normal equally.
  tris_unt_nrmls = normalise(tris_nrmls)

  # triangle centres will be the representative vertices now
  # this reduces the overall dot density by 1/3; since each clique is 3-membered
  tris_cntrs = np.average(ply_vrtx_arr[ply_fc_arr], axis=1)

  #inflate surface by infltn_fctr in Å
  infltn_fctr = infltn
  tris_cntrs = tris_cntrs + infltn_fctr*tris_unt_nrmls

  # perform the actual mapping with reference and map vertex-skipping
  vrtx_spn_maps = spn_map_srfc(tris_cntrs, tris_unt_nrmls*sign[sign_flg], spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq)

  g_spn_maps = np.array(vrtx_spn_maps)

  g_spn_maps_fade = generate_circle_matrix_3d(g_spn_maps, radius_fade, slope_fade)
  g_spn_maps_fade = g_spn_maps_fade.astype(np.float16)

  g_spn_maps_fade_norm = minmax_norm(g_spn_maps_fade)
  # g_flt_spn_maps = g_spn_maps.reshape(g_spn_maps.shape[0], g_spn_maps.shape[-1]**2)
  cmnt_f = """input surface: %s
  point-of-view: outwards
  surface description density: projections were generated with a dot-skipping frequency of %d dot^-1 and dot map frequency of %d dot^-1
  dimensions: every row represents a flattened matrix from a single dot, concatening rows back-to-back; i.e. the 3D matrix is rows x cols^0.2 x cols^0.2
  support distance: %f x %f Angstroms
  resolution: %f Angstrom
  fading radius: %f bins
  fading slope: %f
  """ % (ply_fn, dot_skp_frq, dot_map_frq, spprt_dstnc_x, spprt_dstnc_y, bin_sz, radius_fade, slope_fade)
  np.savez_compressed("results/" + ply_fn[:-4].split("/")[-1] + "_" + sign_flg, \
           coords=tris_cntrs[::dot_skp_frq,:], \
           nrmls=tris_unt_nrmls[::dot_skp_frq,:], \
           maps=g_spn_maps_fade_norm, \
           comments=cmnt_f)
  tst_full_coords_nrmls = [(i,j) for i, j in zip(tris_cntrs[-100:], tris_unt_nrmls[-100:])]
  tst_skpd_coords_nrmls = [(i,j) for i, j in zip(tris_cntrs[::dot_skp_frq,:][-100:], tris_unt_nrmls[::dot_skp_frq,:][-100:])]

  t_end = time.time()
  t_elpsd = t_end - t_strt

  print("mapped %d patches of " % len(vrtx_spn_maps) + ply_fn[:-4] + " in " + str(t_elpsd) + " sec")


################################################################
################################################################
def usage():
  hlp_str = """HECTOR Mapper (CPU) version a0.706  

mandatory input: 
ply:= \t\tPLY file name <str>
dist_x:=\t\tsupport distance alpha <float>
dist_y:=\t\tsupport distance beta <float>
bin:=\t\tbin size <float>
dot_skp:= \t\t\treference vertix-skipping frequency <int>
map_skp:=\t\t\tmapped vertix-skipping frequency <int>
infltn:=\t\tsurface inflation <float>
radius_fade:=\t\tfading radius <float>
slope_fade:=\t\tfading slope <float>
sign:=\t\t"lgnd|rcpt" <str>

# example usage:

python hector_mapper_v10.py in_qry_surf.ply 10.0 20.0 0.2 5 1 0.5 40 0.3 lgnd # for query surface

# or 

python hector_mapper_v10.py in_qry_surf.ply 10.0 20.0 0.2 5 1 0.5 40 0.3 rcpt # for subject surface

# For generating PLY files, use EDTSurf as follows 
# EDTSurf -i input_structure.pdb -o surface_file_name.ply -s 3
# EDTSurf code is available on https://zhanggroup.org/EDTSurf/EDTSurf.zip 
"""
  print(hlp_str)
###################################################################

if __name__ == "__main__":
  
  try:
    parser = argparse.ArgumentParser(description="HECTOR Mapper (CPU) v13")
    parser.add_argument("ply", type=str, help="input PLY file name <str>")
    parser.add_argument("dist_x", type=float, help="support distance alpha <float>", default=12.0)
    parser.add_argument("dist_y", type=float, help="support distance beta <float>", default=6.0)
    parser.add_argument("bin", type=float, help="bin size <float>", default=0.4)
    parser.add_argument("dot_skp", type=int, help="reference vertix-skipping frequency <int>", default=40)
    parser.add_argument("map_skp", type=int, help="mapped vertix-skipping frequency <int>", default=8)
    parser.add_argument("infltn", type=float, help="surface inflation factor in Angstrom <float>", default=0.0)
    parser.add_argument("radius_fade", type=float, help="radius for sigmoid radial fading in bins <float>", default=40.0)
    parser.add_argument("slope_fade", type=float, help="slope of sigmoid radial fading <float>", default=0.3)
    parser.add_argument("sign", type=str, help="\"lgnd\" or \"rcpt\" <str>")


    args = parser.parse_args()
    print(args)
    ply_fn=args.ply
    spprt_dstnc_x=args.dist_x
    spprt_dstnc_y=args.dist_y
    bin_sz=args.bin
    dot_skp_frq=args.dot_skp
    dot_map_frq=args.map_skp
    sign_flg=args.sign
    infltn=args.infltn
    radius_fade=args.radius_fade
    slope_fade=args.slope_fade

    if (sign_flg != "lgnd") and (sign_flg != "rcpt"): 
      usage()
      sys.exit(2)

  except getopt.ArgumentError as e:
    print(str(e))
    usage()
    sys.exit(2)

  try:
    print("HECTOR Mapper (CPU) v13  - Copyright (C) 2018 Mohammad ElGamacy / Max Planck Society - All Rights Reserved")
    srl_rslts = hector_mapper_vctrsd(ply_fn, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq, sign_flg, infltn, radius_fade, slope_fade) 
  
  except Exception as e:
    print(str(e))
    print("--ERROR ENCOUNTERED - REVISE USAGE")
    usage()
    sys.exit(2)
