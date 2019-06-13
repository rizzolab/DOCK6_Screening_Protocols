#! /usr/bin/python
import math, sys
import os.path
from math import sqrt

## data structre that contains the following:
## (1) Name (2) vdw (3) es (4) hbond (5) Int_energy (6) vdw fp (7) es fp (8) hbond fp
## (9) Total score (vdw+es) (10) vdw_fp num res (11) es fp num res (12) hb fp num res

# data structure to store information about each ligand.
class LIG_DATA:
    def __init__(self, name_dock, vdw,es, hbond, int_energy, vdw_fp, es_fp, hb_fp, totalScore, vdw_fp_numres, es_fp_numres, hb_fp_numres):
        self.name_dock     = str(name_dock)
        self.vdw           = float(vdw)
        self.es            = float(es)
        self.int_energy    = float(int_energy)
        self.hbond         = int(hbond)
        self.vdw_fp        = float(vdw_fp)
        self.es_fp         = float(es_fp)
        self.hb_fp         = float(hb_fp)
        self.totalScore    = float(totalScore) 
        self.vdw_fp_numres = int(vdw_fp_numres) 
        self.es_fp_numres  = int(es_fp_numres) 
        self.hb_fp_numres  = int(hb_fp_numres) 
         
    def __cmp__(self, other):
        return cmp(self.totalScore, other.totalScore)

def write_LIG_DATA(data, file):
    file.write( str(data.name_dock) + "," + str(data.vdw) + "," + str(data.es) + "," + str(data.int_energy) + "," + str(data.hbond) + "," + 
                str(data.vdw_fp) + "," + str(data.es_fp) + "," + str(data.hb_fp) + "," + str(data.totalScore) + "," + str(data.vdw_fp_numres) + "," +
                str(data.es_fp_numres) + "," + str(data.hb_fp_numres)+ "\n")

def read_dock_fp_out(file, max_num):
    file1 = open(file, 'r')
    list = []

    for line in file1:
         linesplit = line.split() #split on white space
         if (len(linesplit) > 2):
             if (linesplit[1] == "Name:"):
                 temp_name_dock = str(linesplit[2])
             if (linesplit[1] == "vdw:"):
                 temp_vdw = float(linesplit[2])
             if (linesplit[1] == "es:"):
                 temp_es = float(linesplit[2])
             if (linesplit[1] == "vdw+es:"):
                 temp_totalScore = float(linesplit[2]) 
             if (linesplit[1] == "hbond:"):
                 temp_hbond = int(linesplit[2])
             if (linesplit[1] == "Int_energy:"):
                 temp_internal_energy = float(linesplit[2])
             if (linesplit[1] == "vdw_fp:"):
                 temp_vdw_fp = float(linesplit[2])
             if (linesplit[1] == "es_fp:"):
                 temp_es_fp = float(linesplit[2])
             if (linesplit[1] == "hb_fp:"):
                 temp_hb_fp = float(linesplit[2])
             if (linesplit[1] == "vdw_fp_numres:"):
                 temp_vdw_fp_numres = float(linesplit[2])
             if (linesplit[1] == "es_fp_numres:"):
                 temp_es_fp_numres = float(linesplit[2])
             if (linesplit[1] == "hb_fp_numres:"):
                 temp_hb_fp_numres = float(linesplit[2])
                 data = LIG_DATA(temp_name_dock, temp_vdw, temp_es, temp_hbond, temp_internal_energy, temp_vdw_fp, temp_es_fp,
                                 temp_hb_fp, temp_totalScore, temp_vdw_fp_numres, temp_es_fp_numres, temp_hb_fp_numres)
                 list.append(data)

    file1.close()
    del(file1)
    del(line)

    # sort by the total score (vdw+es)
    list.sort()

#    shortlist = []
#    for i in range(max_num):
#        if i in range(len(list)):
#            shortlist.append(list[i])

    shortlist = []
    zinc_names = []
    i = 0
    j = 0

    while (i < max_num):
        if j in range(len(list)):
            if not (list[j].name_dock in zinc_names):
                zinc_names.append(list[j].name_dock)
                shortlist.append(list[j])
                j += 1
                i += 1
            else:
                j += 1
        else:
            break

    return list, shortlist;

def write_csv(data, filename):
    filehandle = open(filename, 'w')
    filehandle.write("name_dock,vdw,es,int_energy,hbond,vdw_fp,es_fp,hb_fp,totalScore,vdw_fp_numres,es_fp_numres,hb_fp_numres\n")
    for i in range(len(data)):
        write_LIG_DATA(data[i],filehandle)
    filehandle.close()
    return;

def main():
    filename         = sys.argv[1];
    outputfileprefix = sys.argv[2];
    max_number       = sys.argv[3]
    max_num          = int(max_number)

    list, shortlist           = read_dock_fp_out(filename, max_num);
    outputfilename_list       = outputfileprefix + "_sorted_totScore_all.csv";
    outputfilename_shortlist  = outputfileprefix + "_sorted_totScore_" + str(max_num) + ".csv";

    write_csv(list, outputfilename_list)
    write_csv(shortlist, outputfilename_shortlist)
    return; 

main()
