from pymol import cmd, stored
from glob import glob
import sys, os

def average_b(selection):
	stored.tempfactor = 0
	stored.atomnumber = 0
	cmd.iterate(selection, "stored.tempfactor = stored.tempfactor + b")
	cmd.iterate(selection, "stored.atomnumber = stored.atomnumber + 1")
	averagetempfactor = stored.tempfactor / stored.atomnumber

	return averagetempfactor
	# print("Your selection: %s" % selection)
	# print("sum of B factors: %s" % stored.tempfactor)
	# print("number of atoms: %s" % stored.atomnumber)
	# print("average B of '%s': %s" % (selection, averagetempfactor))

# Edit these two variables with the paths to your files
file_glob = "pdbs/*.pdb"

# load the comparison pdbs
cmd.load("1LU4_orig.pdb","wt_1")
cmd.load("1LU4_alt.pdb","wt_2")

print('filename\tRMSD_orig\tRMSD_alt\tpLDDT')

# loop through the files
for file in glob(file_glob):
	cmd.load(file,"decoy")
	rms_1 = cmd.align("wt_1","decoy")[0]
	rms_2 = cmd.align("wt_2","decoy")[0]
	b = average_b("decoy")
	print("%s\t%.3f\t%.3f\t%.3f" % (os.path.basename(file), rms_1, rms_2, b ))
	cmd.delete("decoy")
