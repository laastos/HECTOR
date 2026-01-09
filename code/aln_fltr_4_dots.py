import sys
import copy
import numpy as np
from tensors import *
from joblib import Parallel, delayed
from scipy.spatial import KDTree

################### ATOM CLASS - from PDB str ######################
class atom():

    def __init__(self, pdbatom_str):

        self.idx = int(pdbatom_str[6:11])
        self.name = pdbatom_str[12:16].strip()
        self.resn = pdbatom_str[17:20].strip()
        self.chain_id = pdbatom_str[21].strip()
        self.resid = int(pdbatom_str[22:26])
        x, y, z = float(pdbatom_str[30:38]), float(pdbatom_str[38:46]), \
            float(pdbatom_str[46:54])
        self.coords = np.array([x, y, z])

################### ATOM CLASS - from PDB str ######################

class pdb():

	def __init__(self, atms_lst):
		self.atms_lst = atms_lst
		self.coords = np.array([atm.coords for atm in atms_lst])

	def update(self):
		for i in range(len(self.atms_lst)):
			self.atms_lst[i].coords = self.coords[i]

	def transform(self, \
				  rot_mat = np.array([[1.0, 0.0, 0.0],  \
							    	  [0.0, 1.0, 0.0],  \
							    	  [0.0, 0.0, 1.0]]),\
				  trans_vec = np.array([0.0, 0.0, 0.0]),\
				  centre = True):
		if centre:
			centroid = np.mean(self.coords , axis=0)
			self.coords = self.coords - centroid
		self.coords = np.dot(self.coords , rot_mat) + trans_vec
		self.update()


############################ READ PDB ############################
def read_pdb(fn):
    
    fh_in = open(fn, 'r')
    atms_lst = [atom(line) for line in fh_in.readlines() if "ATOM" in line[:5]]

    return atms_lst
############################ READ PDB ############################

############################ GEN PDB ATOM ############################
def gen_pdb_rcrd(atm, idx=None):

    c1 = "ATOM"
    if idx: c2 = idx 
    else: c2 = atm.idx
    c3 = atm.name if len(atm.name) == 4 else " "+atm.name
    c4 = ' ' #alt-loc
    c5 = atm.resn
    c6 = atm.chain_id
    c7 = atm.resid
    c8 = ' ' #residue insertion code
    c9 = atm.coords[0]
    c10 = atm.coords[1]
    c11 = atm.coords[2]
    c12 = 1.0 # occupancy
    c13 = 0.0 # B factor
    c14 = ' ' # element
    c15 = ' ' # charge
    atm_str = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

    return atm_str
############################ GEN PDB ATOM ############################

############################ WRITE PDB ############################
def write_pdb(o_fn, atms_lst, rmrk=None):

    o_fh = open(("/results/"+o_fn), "w")
    if rmrk:
        o_str = "REMARK " + "\nREMARK ".join(rmrk.split("\n")) + "\n"
    else:
        o_str = ""
    o_str += '\n'.join([gen_pdb_rcrd(atm) for atm in atms_lst])
    o_fh.write(o_str)
    o_fh.close()

    return 0
############################ WRITE PDB ############################

############################ ALIGN COORDS ############################
def kabsch_v2(P, Q):
	"""
	The Kabsch Algorithm aligns two sets of points P and Q in 3D space.
	It returns the optimal rotation matrix U that minimizes the RMSD between P and Q.

	Parameters:
	P -- a numpy array of shape (N,3) representing the first set of points
	Q -- a numpy array of shape (N,3) representing the second set of points

	Returns:
	U -- the optimal rotation matrix
	RMSD -- the root mean square deviation after alignment
	"""

	# Ensure that the two sets of points have the same centroid
	centroid_P = np.mean(P, axis=0)
	centroid_Q = np.mean(Q, axis=0)
	
	P_centered = P - centroid_P
	Q_centered = Q - centroid_Q

	trnsltn_vec = centroid_Q - centroid_P
	tvec_ref = centroid_Q - np.mean(Q_centered, axis=0)
	tvec_mbl = centroid_P - np.mean(P_centered, axis=0)

	# Compute the covariance matrix
	H = np.dot(P_centered.T, Q_centered)

	# Perform Singular Value Decomposition (SVD)
	U, S, Vt = np.linalg.svd(H)

	# Compute the rotation matrix
	d = (np.linalg.det(U) * np.linalg.det(Vt)) < 0.0
	if d:
		S[-1] = -S[-1]
		U[:, -1] = -U[:, -1]

	# Create the rotation matrix
	rotation_matrix = np.dot(U, Vt)

	# Calculate the aligned points
	P_aligned = np.dot(P_centered, rotation_matrix)

	# Calculate the RMSD
	P_alnd_trnsltd = P_aligned + tvec_ref

	pre_rmsd = np.sqrt(np.mean(np.linalg.norm(Q_centered - P_centered, axis=1) ** 2))
	rmsd = np.sqrt(np.mean(np.linalg.norm(Q_centered - P_aligned, axis=1) ** 2))
	rmsd_mbl = np.sqrt(np.mean(np.linalg.norm(Q - P_alnd_trnsltd, axis=1) ** 2))

	return P_alnd_trnsltd, rotation_matrix, tvec_ref, tvec_mbl, rmsd
############################ ALIGN COORDS ############################

############################ ALIGN STRUCT ############################
def anchor_4_dots(s_mbl, q_ref, in_mbl_fn, dmp_alnd_pdb=False):

	mbl_alnd_trnsltd, rotation_matrix, tvec_ref, tvec_mbl, rmsd = kabsch_v2(s_mbl, q_ref)

	try:
		# hit_name = hit[-1][:6]
		# in_mbl_fn = pdbs_dir+hit_name+".pdb"
		in_mbl_pdb = pdb(read_pdb(in_mbl_fn))

		in_mbl_pdb.transform(trans_vec=-tvec_mbl, centre=False)
		in_mbl_pdb.transform(rot_mat=rotation_matrix, centre=False)
		in_mbl_pdb.transform(trans_vec=tvec_ref, centre=False)

		return (in_mbl_pdb, rmsd)
		if dmp_alnd_pdb:
			out_mbl_fn = hit[-1][:6]+"_alnd.pdb"
			write_pdb(out_mbl_fn, in_mbl_pdb.atms_lst)
	except Exception as e:
		return None
############################ ALIGN STRUCT ############################

############################ QUANTIFY OVERLAPS ############################
def calc_ovrlp(qry_atms_lst, \
			   sbj_atms_lst, \
			   slv_rsltn, \
			   slv_n_vxls, \
			   slv_dict, \
			   spprt_dstnc):

	for i in range(len(qry_atms_lst)):
		qry_atms_lst[i].chain_id = "A"

	for i in range(len(sbj_atms_lst)):
		sbj_atms_lst[i].chain_id = "B"

	atms_lst = list()
	atms_lst.extend(qry_atms_lst)
	atms_lst.extend(sbj_atms_lst)

	cmb_cst = cast_bsc(atms_lst)

	sbj_ca_coords = np.array([tmp_res.atms_nms_dict["CA"].coords \
					 for tmp_res in cmb_cst.res_lst \
					 if tmp_res.chain_id == "B" and "CA" in tmp_res.atms_nms_dict])

	# identify proximal interface residues in chain A
	ref_resids = list()
	for ref_res in cmb_cst.res_lst:
		if ref_res.chain_id == "A":
			tmp_ca_coords = ref_res.atms_nms_dict["CA"].coords
			tmp_dists = np.linalg.norm(tmp_ca_coords - sbj_ca_coords, axis=1)
			if np.min(tmp_dists) < 11.0:
				ref_resids.append(ref_res.resid)

	ovrlp = 0
	for ref_resid in ref_resids:

		alnd_atms_lst, trgt_res = init_ref_cst_nodel(cmb_cst, \
												 	 ref_resid, \
												 	 "A")
		qry_atms = [atm for atm in alnd_atms_lst if atm.chain_id == "A"]
		sbj_atms = [atm for atm in alnd_atms_lst if atm.chain_id == "B"]
		qry_msk_mat = prjct_msk_solv( \
									qry_atms, \
									slv_rsltn, \
									slv_n_vxls, \
									slv_dict, \
									spprt_dstnc)
		sbj_msk_mat = prjct_msk_solv( \
									sbj_atms, \
									slv_rsltn, \
									slv_n_vxls, \
									slv_dict, \
									spprt_dstnc)
		ovrlp += np.sum((qry_msk_mat-1)*(sbj_msk_mat-1))

	return ovrlp, len(ref_resids)
############################ QUANTIFY OVERLAPS ############################

def dock_hits(hit, \
			  pdbs_dir, \
			  in_qry_pdb_fn, \
			  ovrlp_co, \
			  intrfc_resids_co, \
			  rmsd_co, \
			  slv_rsltn, \
			  slv_n_vxls, \
			  slv_dict, \
			  spprt_dstnc):

	q_ref = np.array([hit[16:19], \
					  hit[22:25], \
					  hit[16:19] + hit[19:22], \
					  hit[22:25] + hit[25:28]], \
					  dtype=float)

	s_mbl = np.array([hit[4:7], \
					  hit[10:13], \
					  hit[4:7] - hit[7:10], \
					  hit[10:13] - hit[13:16]], \
					  dtype=float)
	in_mbl_fn = pdbs_dir+hit[-1][:6]+".pdb"
	mbl_rms_tup = 		anchor_4_dots( \
								s_mbl, \
								q_ref, \
							in_mbl_fn, \
					 dmp_alnd_pdb=False)
	if mbl_rms_tup == None:
		return None
	elif mbl_rms_tup[1] > rmsd_co:
		return None
	else:
		in_mbl_pdb = mbl_rms_tup[0]
		rmsd = mbl_rms_tup[1]

	qry_atms_lst = read_pdb(in_qry_pdb_fn)
	sbj_atms_lst = in_mbl_pdb.atms_lst
	try:
		ovrlp, intrfc_resids = calc_ovrlp(copy.deepcopy(qry_atms_lst), \
						   copy.deepcopy(sbj_atms_lst), \
						   slv_rsltn, \
						   slv_n_vxls, \
						   slv_dict, \
						   spprt_dstnc)
	except:
		return None

	if ovrlp <= ovrlp_co and \
	   intrfc_resids >= intrfc_resids_co and \
	   rmsd <= rmsd_co:

		out_mbl_fn = hit[-1][:6]+\
					 "_%.3d" % np.random.randint(0,999)+\
					 "_alnd.pdb"
		write_pdb(out_mbl_fn, in_mbl_pdb.atms_lst)
		print("hit:", hit[-1][:6], \
				  " rf:", hit[3], \
				  " rmsd:", rmsd, \
				  " ovrlp:",ovrlp, \
				  " intrfc:", intrfc_resids)
		return [hit, rmsd, ovrlp, intrfc_resids, in_qry_pdb_fn]

###### hits_arr organization ##############################
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
# col28: npz file name
###########################################################

hits_fn = sys.argv[1]
pdbs_dir = sys.argv[2]+"/"
in_qry_pdb_fn = "/data/scaffolds_db/il7ra.pdb"
spprt_dstnc = 22
solv_krnls_fn = r"./solv_krnls_0.50A_24vxl.npz"
slv_rsltn, slv_n_vxls, slv_dict = load_solv_krnls(solv_krnls_fn)
ovrlp_co = 125000
intrfc_resids_co = 25
rmsd_co = 0.5

hits_arr = np.load(hits_fn, allow_pickle=True)

n_cores = 10
rslts_lst = Parallel(n_jobs=n_cores)(delayed(dock_hits)(hit, \
			  pdbs_dir, \
			  in_qry_pdb_fn, \
			  ovrlp_co, \
  			  intrfc_resids_co, \
			  rmsd_co, \
			  slv_rsltn, \
			  slv_n_vxls, \
			  slv_dict, \
			  spprt_dstnc) \
			  for hit in hits_arr)

rslts_arr = np.array([rslt for rslt in rslts_lst if rslt], dtype=object)
np.save("/results/dock_rslts", rslts_arr)
