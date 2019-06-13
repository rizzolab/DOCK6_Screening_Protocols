import sys
from mol2 import *
import os.path
def main():
	file1 = sys.argv[1]
	file2 = sys.argv[2]
	if (not os.path.exists(file1)):
		print "Missing "+file1
	elif (not os.path.exists(file2)):	
		print "Missing "+file2
	else:
        	print "%.10f" % (heavy_atom_RMSD(read_Mol2_file(sys.argv[1])[0], read_Mol2_file(sys.argv[2])[0]))
        return
####################################################################################################

main()
