import sys,os
from Bio.PDB.PDBParser import PDBParser 
from numpy import loadtxt
from Bio.PDB import *
import pross as p 

def trimfile(in_fp, out_fp):
    if not os.path.exists(in_fp):
        return


#rmovef='rm '+'helix_hb_0.15r.dat' + ' hb.log'
#os.system(rmovef)

sys.stdout = open('hb.log', 'w')
## open data
f = open('helix_hb_0.15r.dat', 'w')
pf = open('success.list', 'w')
nf = open('fail.list', 'w')

negative=0
positive=0

## read the PDBID from Protein Data Bank 
## In this example, the protein structures matches: 1) resolution of <=1.5 Angstrom and
## 2) secondary structure has: 10 or more percent of elements are Alpha Helical 
pdbfile = open('all_helix_pdb_0.1_r0.15', 'r')     
temp = pdbfile.readline()
list1 = temp.split()

for name in list1:
	pdbl = PDBList()
	name2=name.lower()
	pdb1 ='pdb'+name2+".ent"
	try:
		pdbl.retrieve_pdb_file(name, pdir=".")
		pdb = p.PDBFile(pdb1)
		chains = pdb.read(as_protein=1)
		for chain in chains.elements:
			(a, b, c, sst_list) = p.rc_ss(chain)
			idx = 0
			for res_class in sst_list:
				if res_class not in ["H"]:     
					chain.elements[idx] = []
				idx += 1
		mol = p.molMol(chains=chains)

		p.write_pdb(mol, pdb1) # write modified PDB to disk
		parser=PDBParser(PERMISSIVE=1)
		structure = parser.get_structure(name,pdb1)

		for model in structure:
			for chain in model:
				for residue in chain:
					try:
						rlist=residue.get_id()
						i = int(rlist[1])
						j=i+4
						residue_1 = chain[i]
						residue_2 = chain[j]
						atom_1 = residue_1['O'] 
						atom_2 = residue_2['N'] 
						distance = atom_1-atom_2 
						f.write(str(distance)+"\t"+str(rlist[1])+"\t"+name+"\n")
					except:
						pass
		os.system('rm %(pdb1)s' % locals())
		pf.write(str(name)+ "\n")
		positive=positive+1
	except:
		negative=negative+1
		nf.write(str(name)+ "\n")
		os.system('rm %(pdb1)s' % locals())
		pass
f.write("#positive: "+ str(positive) +"\n")
f.write("#negative: "+ str(negative) +"\n")
f.close()
nf.close()
sys.exit(0)



