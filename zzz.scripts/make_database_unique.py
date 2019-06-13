#! usr/bin/python
import sys

inputfilename = sys.argv[1]
inputfile = open(inputfilename, 'r')
count = 0
for line in inputfile:
    if ("MOLECULE" in line):
	skip = 1
    elif (count == 0) and ("ZINC" in line):
        count = 1
        name = line.split()[0] + '.mol2'
        outfile = file(str(name), 'w')
        a = line
	outfile.write('@<TRIPOS>MOLECULE\n')
        outfile.write(a)
    else:
         if (count > 0) and ("ZINC" not in line):
            b = line
            outfile.write(b)
         elif (count > 0) and ("ZINC" in line):
            outfile.close()
            name = line.split()[0] + '.mol2'
            outfile = file(str(name), 'w')
            a = line
	    outfile.write('@<TRIPOS>MOLECULE\n')
            outfile.write(a)
outfile.close()
