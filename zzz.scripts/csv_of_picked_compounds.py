#!/usr/bin/env python
import math, sys
import os.path
import subprocess
from math import sqrt

## This script will pull all these descriptors from the header 

# data structure to store information about each ligand.
class LIG_DATA:
    def __init__(self, name_dock, which_score, which_rank, cont_score, cont_vdw, cont_es, hbond, int_energy, fps, vdw_fp, es_fp, hb_fp, totalScore, vdw_fp_numres, es_fp_numres, hb_fp_numres, rot_bond, ph4, vo, tan, hun, dscrpt, slogp, logs, chiralcenters, lipdonors, lipaccept, lipdrug, lipviol, formcharge, ligeff, clustersize, smiles):
        self.name_dock     = str(name_dock)
	self.which_score   = str(which_score)
        self.which_rank    = int(which_rank)
	self.cont_score    = float(cont_score)
        self.cont_vdw      = float(cont_vdw)
        self.cont_es       = float(cont_es)
        self.int_energy    = float(int_energy)
        self.hbond         = int(hbond)
	self.fps	   = float(fps)
        self.vdw_fp        = float(vdw_fp)
        self.es_fp         = float(es_fp)
        self.hb_fp         = float(hb_fp)
        self.totalScore    = float(totalScore) 
        self.vdw_fp_numres = int(vdw_fp_numres) 
        self.es_fp_numres  = int(es_fp_numres) 
        self.hb_fp_numres  = int(hb_fp_numres) 
        self.rot_bond      = int(rot_bond)
	self.ph4           = float(ph4)
	self.vo            = float(vo)
	self.tan           = float(tan)
        self.hun           = float(hun)
        self.dscrpt        = float(dscrpt)
	self.slogp         = float(slogp)
	self.logs          = str(logs)
	self.chiralcenters = int(chiralcenters)
	self.lipdonors     = int(lipdonors)
	self.lipaccept     = int(lipaccept)
	self.lipdrug       = int(lipdrug)
	self.lipviol       = int(lipviol)
	self.formcharge    = int(formcharge)
	self.ligeff        = float(ligeff)
	self.clustersize   = int(clustersize)
        self.smiles        = str(smiles) 

def write_LIG_DATA(data, file):
    file.write( str(data.name_dock) + ","  + str(data.which_score) + ","  + str(data.which_rank) + ","  + str(data.cont_score) + "," + str(data.cont_vdw) + "," + str(data.cont_es) + "," + str(data.int_energy) + "," + str(data.hbond) + "," + str(data.fps) + "," + str(data.vdw_fp) + "," + str(data.es_fp) + "," + str(data.hb_fp) + "," + str(data.totalScore) + "," + str(data.vdw_fp_numres) + "," + str(data.es_fp_numres) + "," + str(data.hb_fp_numres) + "," + str(data.rot_bond) + "," + str(data.ph4) + "," + str(data.vo) + "," + str(data.tan) + "," + str(data.hun) + "," + str(data.dscrpt) +  "," + str(data.slogp) +  "," + str(data.logs) +  "," + str(data.chiralcenters) +  "," + str(data.lipdonors) +  "," + str(data.lipaccept) +  "," + str(data.lipdrug) +  "," + str(data.lipviol) +  "," + str(data.formcharge) +  "," + str(data.ligeff) +  "," + str(data.clustersize) +  "," + str(data.smiles) + "\n")

def read_dock_fp_out(file, max_num):
    file1 = open(file, 'r')
    list = []

    for line in file1:
         linesplit = line.split() #split on white space
         if (len(linesplit) > 2):
             if (linesplit[1] == "Name_DOCK:"):
                 temp_name_dock = str(linesplit[2])
             if (linesplit[1] == "From_List:"):
                 temp_which_score = str(linesplit[2])
             if (linesplit[1] == "List_Rank:"):
                 temp_which_rank = str(linesplit[2])
             if (linesplit[1] == "DOCK_rot_bonds:"):
                 temp_rot_bond = int(linesplit[2])
             if (linesplit[1] == "Continuous_Score:"):
                 temp_cont_score = float(linesplit[2])
             if (linesplit[1] == "Continuous_es_energy:"):
                 temp_cont_es = float(linesplit[2])
             if (linesplit[1] == "Continuous_vdw_energy:"):
                 temp_cont_vdw = float(linesplit[2])
             if (linesplit[1] == "TotalScore_(FPS+DCE):"):
                 temp_totalScore = float(linesplit[2]) 
             if (linesplit[1] == "Num_H-bonds:"):
                 temp_hbond = int(linesplit[2])
             if (linesplit[1] == "Footprint_Similarity_Score:"):
                 temp_fps = float(linesplit[2])
             if (linesplit[1] == "FPS_vdw_fps:"):
                 temp_vdw_fp = float(linesplit[2])
             if (linesplit[1] == "FPS_es_fps:"):
                 temp_es_fp = float(linesplit[2])
             if (linesplit[1] == "FPS_hb_fps:"):
                 temp_hb_fp = float(linesplit[2])
             if (linesplit[1] == "FPS_vdw_fp_numres:"):
                 temp_vdw_fp_numres = float(linesplit[2])
             if (linesplit[1] == "FPS_es_fp_numres:"):
                 temp_es_fp_numres = float(linesplit[2])
             if (linesplit[1] == "FPS_hb_fp_numres:"):
                 temp_hb_fp_numres = float(linesplit[2])
             if (linesplit[1] == "Pharmacophore_Score:"):
	         temp_ph4 = float(linesplit[2])
             if (linesplit[1] == "Property_Volume_Score:"):
                 temp_vo = float(linesplit[2])
             if (linesplit[1] == "Tanimoto_Score:"):
                 temp_tan = float(linesplit[2])
             if (linesplit[1] == "Hungarian_Matching_Similarity_Score:"):
                 temp_hun = float(linesplit[2])
             if (linesplit[1] == "Descriptor_Score:"):
                 temp_dscrpt = float(linesplit[2])
             if (linesplit[1] == "Internal_energy_repulsive:"):
                 temp_internal_energy = float(linesplit[2])
	     if (linesplit[1] == "SlogP:"):
                 temp_slogp = float(linesplit[2])
             if (linesplit[1] == "logS:"):
                 temp_logs = float(linesplit[2])
             if (linesplit[1] == "Num_chiral_centers:"):
                 temp_chiralcenters = float(linesplit[2])
             if (linesplit[1] == "Lipinski_donors:"):
                 temp_lipdonors = float(linesplit[2])
             if (linesplit[1] == "Lipinski_acceptors:"):
                 temp_lipaccept = float(linesplit[2])
             if (linesplit[1] == "Lipinski_druglike:"):
                 temp_lipdrug = float(linesplit[2])
             if (linesplit[1] == "Lipinski_violations:"):
                 temp_lipviol = float(linesplit[2])
             if (linesplit[1] == "Formal_charge:"):
                 temp_formcharge = int(linesplit[2])
             if (linesplit[1] == "Ligand_efficiency:"):
                 temp_ligeff = float(linesplit[2])
             if (linesplit[1] == "Cluster_size:"):
                 temp_clustersize = float(linesplit[2])
             if (linesplit[1] == "SMILES:"):
                 temp_smiles = str(linesplit[2])

                 data = LIG_DATA(temp_name_dock, temp_which_score, temp_which_rank, temp_cont_score, temp_cont_vdw, temp_cont_es, temp_hbond, temp_internal_energy, temp_fps, temp_vdw_fp, temp_es_fp, temp_hb_fp, temp_totalScore, temp_vdw_fp_numres, temp_es_fp_numres, temp_hb_fp_numres, temp_rot_bond, temp_ph4, temp_vo, temp_tan, temp_hun, temp_dscrpt, temp_slogp, temp_logs, temp_chiralcenters, temp_lipdonors, temp_lipaccept, temp_lipdrug, temp_lipviol, temp_formcharge, temp_ligeff, temp_clustersize, temp_smiles)
                 list.append(data)
    
    file1.close()
    del(file1)
    del(line)


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
   

    return list;


def write_csv(data, filename):
    filehandle = open(filename, 'w')
    filehandle.write("name_dock,from_list,rank_in_list,continuous_score,continuous_vdw,continuous_es,internal_energy_repulsive,hbond,fps,vdw_fps,es_fps,hb_fps,totalScore,vdw_fps_numres,es_fps_numres,hb_fps_numres,rot_bond,fms_score,vo_score,tan_score,hms_score,descriptor_score,SlogP,logS,Num_chiral_centers,Lipinski_donors,Lipinski_acceptors,Lipinski_druglike,Lipinski_violations,Formal_charge,Ligand_efficiency,Cluster_size,SMILES\n")
    for i in range(len(data)):
        write_LIG_DATA(data[i],filehandle)
    filehandle.close()
    return;

def checkexist(filename):
  if not os.path.isfile(filename):
    print " The file [ %s ] does not exist" % filename
    quit()



def main():
    filename         = sys.argv[1]
    outputfileprefix = sys.argv[2]
    max_number       = sys.argv[3]
    max_num          = int(max_number)
    checkexist(filename)

    list                      = read_dock_fp_out(filename, max_num);
    outputfilename_list       =  outputfileprefix + "_picked_compounds_" + str(max_num) + ".csv";

    write_csv(list, outputfilename_list)
    return; 

main()
