# This script reads in a mol2 file, and writes another one out to make sure 
# the grid program reads it in correctly

import sys
from mol2 import *
import os.path
def main():
	inmol2 = sys.argv[1]
	outmol2 = sys.argv[2]
	in_molecule = read_Mol2_file(inmol2)[0]
	write_mol2(in_molecule,outmol2)

        return
####################################################################################################

main()
