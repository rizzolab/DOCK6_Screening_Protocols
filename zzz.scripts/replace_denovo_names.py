### This script replaces the arbitrary zinc names of de novo molecules with unique identifiers
#!/usr/bin/python

import math, sys
import os.path
import subprocess
import fileinput
import string
import shutil

counter=1
with open(sys.argv[1],"r") as fin:
    with open(sys.argv[2],"w")as fout:
        for line in fin:
	    linesplit = line.split()
            #if (len(linesplit) == 1 and linesplit[0] == "@<TRIPOS>MOLECULE"):
            if (len(linesplit) == 1 and line.startswith("ZINC")):
		filedata = line.replace( str(line) , sys.argv[3] + "_"+ str(counter) + " \n ")
		fout.write(filedata)
		counter+=1
		#print "counter" + str(counter)
	    
            else:        
                fout.write(line)
fin.close()
fout.close()
