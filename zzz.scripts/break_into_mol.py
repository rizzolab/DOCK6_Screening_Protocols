#!/usr/bin/env python
import sys,os

def checkexist(filename):
  if not os.path.isfile(filename):
    print " The file [ %s ] does not exist" % filename
    quit()


inputfilename = sys.argv[1]
checkexist(inputfilename)
inputfile = open(inputfilename, 'r')
count = 0
for line in inputfile:
    if (count == 0) and ("##########                                Name:" in line):
        count = 1
        name = line.split()[2] + '.mol2'
        outfile = file(str(name), 'w')
        a = line
        outfile.write(a)
    else:
         if (count > 0) and ("##########                                Name:" not in line):
            b = line
            outfile.write(b)
         elif (count > 0) and ("##########                                Name:" in line):
            outfile.close()
            name = line.split()[2] + '.mol2'
            outfile = file(str(name), 'w')
            a = line
            outfile.write(a)
outfile.close()
