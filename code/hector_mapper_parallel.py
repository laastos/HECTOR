#!/usr/bin/env python
#Copyright (C) 2018 Mohammad ElGamacy / Max Planck Society - All Rights Reserved
#Parallelized version using joblib

import sys
import time
import argparse
import numpy as np
from joblib import Parallel, delayed
import multiprocessing

# Import functions from original hector_mapper
from hector_mapper import parse_ply, normalise, generate_circle_matrix_3d, minmax_norm, logistic

###############################################################
# Parallelized spin map generator
###############################################################

def process_single_vertex(vrtx, n, srfc_vrtx_arr, dot_map_frq, spprt_dstnc_x, spprt_dstnc_y, bin_sz, half_spprt_dstnc, tot_bins_x, tot_bins_y):
    """Process a single reference vertex to generate its spin map

    Args:
        vrtx: reference vertex coordinates
        n: reference vertex normal
        srfc_vrtx_arr: all surface vertices
        dot_map_frq: mapped vertex-skipping frequency
        spprt_dstnc_x: support distance alpha
        spprt_dstnc_y: support distance beta
        bin_sz: bin size
        half_spprt_dstnc: half of support distance y
        tot_bins_x: total bins in x direction
        tot_bins_y: total bins in y direction

    Returns:
        spn_img: the generated spin image for this vertex
    """
    spn_img = np.zeros((tot_bins_y, tot_bins_x), dtype=np.float64)

    # Array of "x - p"
    xp_arr = srfc_vrtx_arr[::dot_map_frq, :] - vrtx
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

    x_arr = (half_spprt_dstnc - lmt_lmt_beta_arr)/bin_sz
    y_arr = lmt_lmt_alpha_arr/bin_sz
    i_arr = x_arr.astype(int)
    j_arr = y_arr.astype(int)

    for i, j in zip(i_arr, j_arr):
        spn_img[i][j] += 1

    return spn_img[1:-1, :-2]


def spn_map_srfc_parallel(srfc_vrtx_arr, srfc_nrml_arr, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq, n_jobs=-1):
    """Parallelized version of the core fingerprinting function

    Args:
        srfc_vrtx_arr (numpy.ndarray): a matrix of input vertices position vectors
        srfc_nrml_arr (numpy.ndarray): a matrix of input vertices surface normal vectors
        spprt_dstnc_x (float): support distance alpha; radius of the cylinderical basis in Angstroms
        spprt_dstnc_y (float): support distance beta; height of the cylinderical basis in Angstroms
        bin_sz (float): size of the bin; i.e. resolution of the map in Angstroms
        dot_skp_frq (int): dot-skipping frequency; the sparsity of sampling for references vertices
        dot_map_frq (int): dot-skipping frequency; the sparsity of sampling for mapped vertices
        n_jobs (int): number of parallel jobs (-1 = all CPUs)

    Returns
        vrtx_spn_maps (list(numpy.ndarray)): A list of the final fingerprints
    """

    tot_bins_x = int(spprt_dstnc_x/bin_sz) + 2
    tot_bins_y = int(spprt_dstnc_y/bin_sz) + 2
    half_spprt_dstnc = spprt_dstnc_y/2

    # Get reference vertices and normals
    ref_vertices = srfc_vrtx_arr[::dot_skp_frq, :]
    ref_normals = srfc_nrml_arr[::dot_skp_frq, :]

    total_vertices = len(ref_vertices)
    print(f"Generating local maps for {total_vertices} reference vertices using {n_jobs if n_jobs > 0 else multiprocessing.cpu_count()} cores...")

    # Process vertices in parallel
    vrtx_spn_maps = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_single_vertex)(
            vrtx, n, srfc_vrtx_arr, dot_map_frq,
            spprt_dstnc_x, spprt_dstnc_y, bin_sz,
            half_spprt_dstnc, tot_bins_x, tot_bins_y
        )
        for vrtx, n in zip(ref_vertices, ref_normals)
    )

    return vrtx_spn_maps


###############################################################
def hector_mapper_vctrsd_parallel(ply_fn, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq, sign_flg, infltn, radius_fade, slope_fade, n_jobs=-1):
    """The parallelized mapping and PLY-preprocessing wrapper function

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
        n_jobs (int): number of parallel jobs (-1 = all CPUs)
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

    # perform the actual mapping with reference and map vertex-skipping (PARALLELIZED)
    vrtx_spn_maps = spn_map_srfc_parallel(tris_cntrs, tris_unt_nrmls*sign[sign_flg],
                                          spprt_dstnc_x, spprt_dstnc_y, bin_sz,
                                          dot_skp_frq, dot_map_frq, n_jobs)

    g_spn_maps = np.array(vrtx_spn_maps)

    g_spn_maps_fade = generate_circle_matrix_3d(g_spn_maps, radius_fade, slope_fade)
    g_spn_maps_fade = g_spn_maps_fade.astype(np.float16)

    g_spn_maps_fade_norm = minmax_norm(g_spn_maps_fade)

    cmnt_f = """input surface: %s
    point-of-view: outwards
    surface description density: projections were generated with a dot-skipping frequency of %d dot^-1 and dot map frequency of %d dot^-1
    dimensions: every row represents a flattened matrix from a single dot, concatening rows back-to-back; i.e. the 3D matrix is rows x cols^0.2 x cols^0.2
    support distance: %f x %f Angstroms
    resolution: %f Angstrom
    fading radius: %f bins
    fading slope: %f
    """ % (ply_fn, dot_skp_frq, dot_map_frq, spprt_dstnc_x, spprt_dstnc_y, bin_sz, radius_fade, slope_fade)

    np.savez_compressed("/results/" + ply_fn[:-4].split("/")[-1] + "_" + sign_flg, \
                        coords=tris_cntrs[::dot_skp_frq,:], \
                        nrmls=tris_unt_nrmls[::dot_skp_frq,:], \
                        maps=g_spn_maps_fade_norm, \
                        comments=cmnt_f)

    t_end = time.time()
    t_elpsd = t_end - t_strt

    print("mapped %d patches of " % len(vrtx_spn_maps) + ply_fn[:-4] + " in " + str(t_elpsd) + " sec")


################################################################
################################################################
def usage():
    hlp_str = """HECTOR Mapper (CPU-Parallel) version a0.706-parallel

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

optional:
n_jobs:=\t\tnumber of parallel jobs (default: -1 = all CPUs) <int>

# example usage:

python hector_mapper_parallel.py in_qry_surf.ply 10.0 20.0 0.2 5 1 0.5 40 0.3 lgnd --n_jobs 8

# For generating PLY files, use EDTSurf as follows
# EDTSurf -i input_structure.pdb -o surface_file_name.ply -s 3
# EDTSurf code is available on https://zhanggroup.org/EDTSurf/EDTSurf.zip
"""
    print(hlp_str)
###################################################################

if __name__ == "__main__":

    try:
        parser = argparse.ArgumentParser(description="HECTOR Mapper (CPU-Parallel) v13")
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
        parser.add_argument("--n_jobs", type=int, help="number of parallel jobs (-1 = all CPUs)", default=-1)

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
        n_jobs=args.n_jobs

        if (sign_flg != "lgnd") and (sign_flg != "rcpt"):
            usage()
            sys.exit(2)

    except Exception as e:
        print(str(e))
        usage()
        sys.exit(2)

    try:
        print("HECTOR Mapper (CPU-Parallel) v13  - Copyright (C) 2018 Mohammad ElGamacy / Max Planck Society - All Rights Reserved")
        print(f"Using {n_jobs if n_jobs > 0 else multiprocessing.cpu_count()} CPU cores")
        srl_rslts = hector_mapper_vctrsd_parallel(ply_fn, spprt_dstnc_x, spprt_dstnc_y, bin_sz,
                                                   dot_skp_frq, dot_map_frq, sign_flg, infltn,
                                                   radius_fade, slope_fade, n_jobs)

    except Exception as e:
        print(str(e))
        print("--ERROR ENCOUNTERED - REVISE USAGE")
        usage()
        sys.exit(2)
