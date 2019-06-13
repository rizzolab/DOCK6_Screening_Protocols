### This script will renumber the residue numbers from de novo growth (column 7) so that
### all molecules are one residue for the purposes of moe clustering 

import shutil
import os,sys

with open(sys.argv[1], "r") as in_file:
    buf = in_file.readlines()

with open(sys.argv[2], "w") as out_file:
    bool = 0
    for line in buf:
       if line.startswith('@<TRIPOS>ATOM'):
            out_file.write(str(line))
            bool = 1
       if line.startswith('@<TRIPOS>BOND'):
            bool = 0

       if (bool == 1 ) and line.startswith(' ') :
           line = line[0:55] + "1" + line[57:]
           out_file.write(str(line))
       if (bool == 0):
            out_file.write(str(line))    



in_file.close()
out_file.close()
