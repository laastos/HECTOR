# An example presented here shows the selection of design scaffolds that are complementary to the target epitope on the surface of the interleukin-7 receptor-α (IL-7Rα)

###################################################################################################
# Step 1: Dot surface of the query protein (IL-7Rα) is forward-mapped. Later (in step 2) for comparison with subject maps, the query maps will be inverted.  

# For generating PLY files, use EDTSurf as follows 
# EDTSurf -i input_structure.pdb -o surface_file_name.ply -s 3
# EDTSurf code is available on https://zhanggroup.org/EDTSurf/EDTSurf.zip

python hector_mapper.py /data/il7ra.ply 10 20 0.2 5 1 0.5 40 0.3 rcpt

# Output file: il7ra_rcpt.npz
###################################################################################################

###################################################################################################
# Step 2: Coordinates of three points at the target epitope are used as input for the HECTOR search. Pairwise combinations of query maps corresponding to the selected coordinates are compared against all subject maps from the scaffold database (two-vs-all search). Pairs of subject maps that satisfy the inter-patch distance requirement and have an average R-factor lower than -0.82 are identified as HECTOR hits.

python maps_analysis_pairs_vs_all.py /data/scaffolds_db/ il7ra_rcpt.npz -0.82

# Output file: srch_rslts.npy
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
###################################################################################################

# Step 3: Identified scaffolds from step 2 are docked against target epitope. Hits are filtered by RMSD, surface overlap and number of interface residues.

python aln_fltr_4_dots.py /results/srch_rslts.npy /data/scaffolds_db

# Output file: dock_rslts.npy containing information on hit, rmsd, surface overlap, interface residues, and query PDB
###################################################################################################