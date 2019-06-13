#! /usr/bin/python
import math, sys, string, csv
import os.path
from math import sqrt


### Combine the two csvs from MOE and DOCK, then sort them based on total score, keeping molecules
### within same cluster together.

### DOCK csv contains the following:
###  name_dock,vdw,es,int_energy,hbond,vdw_fp,es_fp,hb_fp,totalScore,vdw_fp_numres,es_fp_numres,hb_fp_numres
###  ZINC05935212,-32.205696,-0.606634,126.884537,0,4.105981,5.769107,2.44949,-32.812328,118,118,118

### MOE csv contains the following:
###  name,SMILES,CLUSTER_NO,Weight,b_rotN,lip_don,lip_acc,lip_druglike,lip_violation,SlogP,FCharge,logS,chiral
###  ZINC05935212,Nc1c2c3c(ccc2cc2c1ccc1c2cccc1)cccc3,2,293.36899,0,2,1,1,1,5.8815999,0,-8.61728,0


### Data structure to store information about each ligand
class LIG_DATA:
    def __init__(self, name_dock, name_moe, cluster, clustersize, totalScore, dce_sum, dce_vdw, dce_es, int_energy, fps_sum, fps_vdw, fps_es, fps_hb, vdw_fp_numres, es_fp_numres, hb_fp_numres, hbond, weight, b_rotN, chiral, lip_don, lip_acc, lip_druglike, lip_violation, SlogP, FCharge, logS, lig_eff, smiles):
        self.name_dock     = str(name_dock)
        self.name_moe      = str(name_moe)
        self.cluster       = int(cluster)
        self.clustersize   = int(clustersize)
        self.totalScore    = float(totalScore)
        self.dce_sum       = float(dce_sum)
        self.dce_vdw       = float(dce_vdw)
        self.dce_es        = float(dce_es)
        self.int_energy    = float(int_energy)
        self.fps_sum       = float(fps_sum)
        self.fps_vdw       = float(fps_vdw)
        self.fps_es        = float(fps_es)
        self.fps_hb        = float(fps_hb)
        self.vdw_fp_numres = int(vdw_fp_numres)
        self.es_fp_numres  = int(es_fp_numres)
        self.hb_fp_numres  = int(hb_fp_numres)
        self.hbond         = int(hbond)
        self.weight        = float(weight)
        self.b_rotN        = int(b_rotN)
        self.chiral        = int(chiral)
        self.lip_don       = int(lip_don)
        self.lip_acc       = int(lip_acc)
        self.lip_druglike  = int(lip_druglike)
        self.lip_violation = int(lip_violation)
        self.SlogP         = float(SlogP)
        self.FCharge       = int(FCharge)
        self.logS          = float(logS)
        self.lig_eff       = float(lig_eff)
        self.smiles        = str(smiles)
    def __cmp__(self, other):
        return cmp(self.totalScore, other.totalScore)


### Some sort definitions
def by_fps_vdw(x, y):
    return cmp(x.fps_vdw, y.fps_vdw)
def by_fps_es(x, y):
    return cmp(x.fps_es, y.fps_es)
def by_fps_sum(x, y):
    return cmp(x.fps_sum, y.fps_sum)
def by_dce_sum(x, y):
    return cmp(x.dce_sum, y.dce_sum)
def by_totalScore(x, y):
    return cmp(x.totalScore, y.totalScore)

### Functions to print and write the LIG_DATA object
def print_LIG_DATA(data):
    print str(data.name_dock) + "," + str(data.name_moe) + "," + str(data.cluster) + "," + str(data.clustersize) + "," + str(data.totalScore) + "," + str(data.dce_sum) + "," + str(data.dce_vdw) + "," + str(data.dce_es) + "," + str(data.int_energy) + "," + str(data.fps_sum) + "," + str(data.fps_vdw) + "," + str(data.fps_es) + "," + str(data.fps_hb) + "," + str(data.vdw_fp_numres) + "," + str(data.es_fp_numres) + "," + str(data.hb_fp_numres) + "," + str(data.hbond) + "," + str(data.weight) + "," + str(data.b_rotN) + "," + str(data.chiral) + "," + str(data.lip_don) + "," + str(data.lip_acc) + "," + str(data.lip_druglike) + "," + str(data.lip_violation) + "," + str(data.SlogP) + "," + str(data.FCharge) + "," + str(data.logS) + "," + str(data.lig_eff) + "," + str(data.smiles)

def write_LIG_DATA(data,file):
    file.write(str(data.name_dock) + "," + str(data.name_moe) + "," + str(data.cluster) + "," + str(data.clustersize) + "," + str(data.totalScore) + "," + str(data.dce_sum) + "," + str(data.dce_vdw) + "," + str(data.dce_es) + "," + str(data.int_energy) + "," + str(data.fps_sum) + "," + str(data.fps_vdw) + "," + str(data.fps_es) + "," + str(data.fps_hb) + "," + str(data.vdw_fp_numres) + "," + str(data.es_fp_numres) + "," + str(data.hb_fp_numres) + "," + str(data.hbond) + "," + str(data.weight) + "," + str(data.b_rotN) + "," + str(data.chiral) + "," + str(data.lip_don) + "," + str(data.lip_acc) + "," + str(data.lip_druglike) + "," + str(data.lip_violation) + "," + str(data.SlogP) + "," + str(data.FCharge) + "," + str(data.logS) + "," + str(data.lig_eff) + "," + str(data.smiles) + "\n")

####################################################################################################

def read_csv_file(inputfilename):

    #### Open the file and read all lines
    file1 = open(inputfilename,'r')
    lines = file1.readlines()

    ### Declare an empty matrix
    matrix = []
    row = 1
    heading = 0
    smiles = 0

    ### For every line in the file
    for line in lines:

        ### Split the file by commas
        text = line.split(",")
        row_of_m = []

        ### Check to see if this is the header line, and whether it is the MOE or DOCK file
        if (text[0] == "name") or (text[0] == "name_dock"):
            heading = 1
            if (text[1] == 'SMILES'):
                smiles = 1

        ### If it is not the header
        else:

            ### For every column in the line
            for column in range(len(text)):

                ### If the field is empty, assign it the value '0'
                if (len(text[column]) == 0):
                   text[column] = '0'
                ### The first column in each file is the name
                if (column == 0):
                   entry = str(text[column])
                ### Only the MOE file has another string
                elif ((smiles == 1) and (column == 1)):
                    entry = str(text[column])
                ### All other entries can be interpreted as floats
                else:
                   entry = float(text[column])
                row_of_m.append(entry)
            
            ### Append the finished row to the matrix
            matrix.append(row_of_m)
            row = row + 1

    file1.close()

    return matrix

####################################################################################################

def write_csv(data,filename):

    ### Open up the filehandle
    filehandle = open(filename,'w')

    ### Write the header information
    filehandle.write("name_dock,name_moe,cluster,clustersize,totalScore(fps+dce),dce_sum(vdw+es),dce_vdw,dce_es,int_energy,fps_sum(vdw+es),fps_vdw,fps_es,fps_hb,vdw_fp_numres,es_fp_numres,hb_fp_numres,hbond,weight,b_rotN,chiral,lip_don,lip_acc,lip_druglike,lip_violation,SlogP,FCharge,logS,lig_eff,smiles\n")

    ### Write every LIG_DATA object in the given list
    for i in range(len(data)):
        write_LIG_DATA(data[i],filehandle)
    filehandle.close()

    return

####################################################################################################

def combine_and_sort(dock_inputfilename,moe_inputfilename):

    ### Read the two csv files and put them in matrix format
    dock_matrix = read_csv_file(dock_inputfilename)
    moe_matrix  = read_csv_file(moe_inputfilename)

### DOCK csv contains the following:
###      0      1   2     3        4     5      6     7       8           19           10           11
###  name_dock,vdw,es,int_energy,hbond,vdw_fp,es_fp,hb_fp,totalScore,vdw_fp_numres,es_fp_numres,hb_fp_numres
###  ZINC05935212,-32.205696,-0.606634,126.884537,0,4.105981,5.769107,2.44949,-32.812328,118,118,118
### MOE csv contains the following:
###      0       1       2        3       4      5       6       7              8          9      10    11    12
###  name_moe,SMILES,CLUSTER_NO,Weight,b_rotN,lip_don,lip_acc,lip_druglike,lip_violation,SlogP,FCharge,logS,chiral
###  ZINC05935212,Nc1c2c3c(ccc2cc2c1ccc1c2cccc1)cccc3,2,293.36899,0,2,1,1,1,5.8815999,0,-8.61728,0

    ### Create an empty list that will be populated by the combined csv elements
    list = []
    i = 0 # dockmatrix count
    j = 0 # moematrix count

    if range(len(dock_matrix)) != range(len(moe_matrix)) :
       print "UNEQUAL ENTRIES IN DOCK AND MOE CSVS"

    while i < len(dock_matrix):
        temp_name_dock         = str(dock_matrix[i][0])
        temp_dce_vdw           = '%.3f' % (float(dock_matrix[i][1]))
        temp_dce_es            = '%.3f' % (float(dock_matrix[i][2]))
        temp_dce_sum           = '%.3f' % (float(temp_dce_vdw)+float(temp_dce_es))
        temp_int_energy        = '%.3f' % (float(dock_matrix[i][3]))
        temp_hbond             = int(dock_matrix[i][4])
        temp_fps_vdw           = '%.3f' % (float(dock_matrix[i][5]))
        temp_fps_es            = '%.3f' % (float(dock_matrix[i][6]))
        temp_fps_hb            = '%.3f' % (float(dock_matrix[i][7]))
        temp_fps_sum           = '%.3f' % (float(temp_fps_vdw)+float(temp_fps_es))
        temp_totalScore        = '%.3f' % (float(temp_dce_sum)+float(temp_fps_sum))
        temp_vdw_fp_numres     = int(dock_matrix[i][9])
        temp_es_fp_numres      = int(dock_matrix[i][10])
        temp_hb_fp_numres      = int(dock_matrix[i][11])
        temp_name_moe          = str(moe_matrix[j][0])
        temp_smiles            = str(moe_matrix[j][1])
        temp_cluster           = int(moe_matrix[j][2])
        temp_clustersize       = 0
        temp_weight            = '%.3f' % (float(moe_matrix[j][3]))
        temp_b_rotN            = int(moe_matrix[j][4])
        temp_lip_don           = int(moe_matrix[j][5])
        temp_lip_acc           = int(moe_matrix[j][6])
        temp_lip_druglike      = int(moe_matrix[j][7])
        temp_lip_violation     = int(moe_matrix[j][8])
        temp_SlogP             = '%.3f' % (float(moe_matrix[j][9]))
        temp_FCharge           = int(moe_matrix[j][10])
        temp_logS              = float(moe_matrix[j][11])
        temp_chiral            = int(moe_matrix[j][12])
        temp_lig_eff           = '%.3f' % (float(dock_matrix[i][8]/moe_matrix[j][3])* (-1))

        ### If the DOCK name and MOE name match for this entry, add all data to a LIG_DATA object
        if (temp_name_dock == temp_name_moe) :
            data = LIG_DATA(temp_name_dock, temp_name_moe, temp_cluster, temp_clustersize, temp_totalScore, temp_dce_sum, temp_dce_vdw, temp_dce_es, temp_int_energy, temp_fps_sum, temp_fps_vdw, temp_fps_es, temp_fps_hb, temp_vdw_fp_numres, temp_es_fp_numres, temp_hb_fp_numres, temp_hbond, temp_weight, temp_b_rotN, temp_chiral, temp_lip_don, temp_lip_acc, temp_lip_druglike, temp_lip_violation, temp_SlogP, temp_FCharge, temp_logS, temp_lig_eff, temp_smiles)

            ### Then append that object to the list
            list.append(data)
            i = i + 1
            j = j + 1

        else:
            print "Mismatch occurred: Dock count " + str(i) + "; Moe count " + str(j)
            j = j + 1

    print " done making combined csv "

    ### Declare some additional empty lists for sorting purposes
    sorted_fps_vdw    = []
    sorted_fps_es     = []
    sorted_fps_sum    = []
    sorted_dce_sum    = []
    sorted_totalScore = []

    ### Sort by the vdw footprint
    list.sort(by_fps_vdw)
    for i in range(len(list)):
        sorted_fps_vdw.append(list[i])
    print "Done sorting by fps_vdw"

    ### Sort by the es footprint
    list.sort(by_fps_es)
    for i in range(len(list)):
        sorted_fps_es.append(list[i])
    print "Done sorting by fps_es"

    ### Sort by the footprint sum
    list.sort(by_fps_sum)
    for i in range(len(list)):
        sorted_fps_sum.append(list[i])
    print "Done sorting by fps_sum"

    ### Sort by the dock cartesian energy sum
    list.sort(by_dce_sum)
    for i in range(len(list)):
        sorted_dce_sum.append(list[i])
    print "Done sorting by dce_sum"

    ### Sort by the total score (fps+dce)
    list.sort(by_totalScore)
    for i in range(len(list)):
        sorted_totalScore.append(list[i])
    print "Done sorting by totalScore"

    return sorted_fps_vdw, sorted_fps_es, sorted_fps_sum, sorted_dce_sum, sorted_totalScore

####################################################################################################

def choose_cluster_heads_and_sort_by_family(data,outfileprefix,outfilesuffix):

    ### Create some empty arrays and declare some variables
    families = []
    clusterheads = []
    size = []
    clusters_visited = []
    i = 0
    j = 0
    number_clusterheads = 0

    
    ### For every element of the list of LIG_DATA object
    for i in range(len(data)):

        ### If the cluster number has been seen before
        if (data[i].cluster in clusters_visited):
            skip = "yes"

        ### If the cluster number has NOT been seen before
        else:
             ### Add the cluster number to the list of "seen" clusters 
             clusters_visited.append(data[i].cluster)
             j = 0
             flag = 0

             ### Then iterate over every element in the sorted list
             for j in range(len(data)):

                 ### If that element has the same cluster number as the current cluster number
                 if ((data[j].cluster == data[i].cluster)):

                     ### Then add that obect to the families list and increment flag
                     temp = data[j]
                     families.append(temp)
                     flag = flag + 1

                     ### If this is the clusterhead, add it to the clusterheads list, and increment the count
                     if (j == i):
                        number_clusterheads = number_clusterheads + 1
                        clusterheads.append(temp)

             k = 1
             no = flag + 1

             ### Populate the size array with the clustersize for each member of the families array
             while (k < no) :
                 size.append(flag)
                 k = k + 1

    p = 0
    q = 0

    ### For every member of families, update the cluster size field
    for p in range(len(families)):
        families[p].clustersize = size[p]

    ### For every member of clusterheads, update the cluster size field
    col_in_family = 0
    for q in range(len(clusterheads)):
        clusterheads[q].clustersize = families[col_in_family].clustersize
        col_in_family = col_in_family +  int(clusterheads[q].clustersize)


    ### Write the complete sorted csv files
    outputfilename_families = outfileprefix + "_sorted_" + outfilesuffix + "_families.csv"
    write_csv(families, outputfilename_families)

    ### Write the complete sorted csv files
    outputfilename_clusterheads = outfileprefix + "_sorted_" + outfilesuffix + "_clusterheads.csv"
    write_csv(clusterheads, outputfilename_clusterheads)

    return

####################################################################################################

def main():
    #### Get the command line arguments
    dock_inputfilename = sys.argv[1]    # DOCK csv file, top N ligands
    moe_inputfilename  = sys.argv[2]    # MOE csv file, top N ligands
    outfileprefix      = sys.argv[3]


    #### Combine the DOCK csv with the MOE csv and sort by different metrics
    sorted_fps_vdw, sorted_fps_es, sorted_fps_sum, sorted_dce_sum, sorted_totalScore = combine_and_sort(dock_inputfilename,moe_inputfilename)


    ### Update the order of each of the sorted files to group stuff by clusters, and print out
    ### the cluster heads
    choose_cluster_heads_and_sort_by_family(sorted_fps_vdw,outfileprefix,"fps_vdw")
    choose_cluster_heads_and_sort_by_family(sorted_fps_es,outfileprefix,"fps_es")
    choose_cluster_heads_and_sort_by_family(sorted_fps_sum,outfileprefix,"fps_sum")
    choose_cluster_heads_and_sort_by_family(sorted_dce_sum,outfileprefix,"dce_sum")
    choose_cluster_heads_and_sort_by_family(sorted_totalScore,outfileprefix,"totalScore")


    ### Write the complete sorted csv files
    #outputfilename_sorted_fps_vdw = outfileprefix + "_sorted_fps_vdw_bestfirst.csv"
    #write_csv(sorted_fps_vdw, outputfilename_sorted_fps_vdw)
    #outputfilename_sorted_fps_es = outfileprefix + "_sorted_fps_es_bestfirst.csv"
    #write_csv(sorted_fps_es, outputfilename_sorted_fps_es)
    #outputfilename_sorted_fps_sum = outfileprefix + "_sorted_fps_sum_bestfirst.csv"
    #write_csv(sorted_fps_sum, outputfilename_sorted_fps_sum)
    #outputfilename_sorted_dce_sum = outfileprefix + "_sorted_dce_sum_bestfirst.csv"
    #write_csv(sorted_dce_sum, outputfilename_sorted_dce_sum)
    #outputfilename_sorted_totalScore = outfileprefix + "_sorted_totalScore_bestfirst.csv"
    #write_csv(sorted_totalScore, outputfilename_sorted_totalScore)


    return

main()

