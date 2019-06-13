#!/usr/bin/env python


#This programs reads in the file ${system}.rec.clust.close.pdb
#and adds a TER card after each line. This is only a cosmetic
#fix so that chimera does not try to place bonds in between the
#spheres.

import os,sys

def usage():
  if len(sys.argv)==1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print "add_TER.py [INPUT pdb] [OUTPUT pdb]"
    sys.exit()

def main():
  usage()
  if not os.path.isfile(sys.argv[1]):
    print "The pdb file %s does not exist" % sys.argv[1]
    sys.exit()

  clustpdb=open(sys.argv[1], "r")
  outpdb=open(sys.argv[2], "w")

  for line in clustpdb:
    outpdb.write(line)
    outpdb.write("TER\n") 

  clustpdb.close()
  outpdb.close()

main()
