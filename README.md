# Data_Mining_for_Hydrogen_Bond_Distance-Python
This is an example to extract hydrogen bond distance of protein structures.

First, the protein structures having the following critera were selected.
1) resolution of <=1.5 Angstrom and
2) secondary structure has: 10 or more percent of elements are Alpha Helical 
The PDBid of the selected protein structures are saved in "all_helix_pdb_0.1_r0.15"

Then, the strict_calculate_helix.py will extract the N-O of helix structures, and writes 
the resulting information in "helix_hb_0.15r.dat"

Finally, the resulting distribution is shown in hydrogen_bond_distance.png.

