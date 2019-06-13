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
    if (count == 0) and ("###  Molecule:" in line):
        count = 1
        name = line.split()[2] + '.txt'
        outfile = file(str(name), 'w')
        a = line
        outfile.write(a)
    else:
         if (count > 0) and ("###  Molecule:" not in line):
            b = line
            outfile.write(b)
         elif (count > 0) and ("###  Molecule:" in line):
            outfile.close()
            name = line.split()[2] + '.txt'
            outfile = file(str(name), 'w')
            a = line
            outfile.write(a)
outfile.close()

#fileName = sys.argv[2]
#ZINClist = open(fileName,'r')
#ZINCids = ZINClist.readlines()
#ZINClist.close()
#
#for i in range(len(ZINCids)):
#	ZINCids[i] = ZINCids[i].strip()
#
#ZINCids.sort()
#print ZINCids
