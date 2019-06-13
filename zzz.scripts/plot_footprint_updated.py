#!/usr/bin/env python

import sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator, FormatStrFormatter


####################################################################################################

def identify_residues(filename, cutoff, max_res):

    ### Read input file of zinc ids
    inputfile = open(filename,'r')
    zinc_ids = inputfile.readlines()
    inputfile.close()

    ### Read in the first footprint txt file
    temp = open(zinc_ids[0].strip()+'.txt','r')
    templines = temp.readlines()
    temp.close()

    ### Count the number of residues in footprint.txt file
    num_res = 0
    
    for templine in templines:
        templinesplit = templine.split()
        if (len(templinesplit) == 8):
            if (templinesplit[0] != 'resname'):
                num_res += 1

    ### Create an array to analyze the footprint info
    fp_array = [[0 for i in range(2)] for j in range(num_res)]

    for i in range(num_res):
        fp_array[i][0] = i


    ### Go through again and count how many times each residue is above the cutoff
    for zinc_id in zinc_ids:
        footprint = open(zinc_id.strip()+'.txt','r')
        lines = footprint.readlines()
        footprint.close()

        count = -1
        for line in lines:
            linesplit = line.split()
            if (len(linesplit) == 8 and (linesplit[0]) != 'resname'):
                count += 1
                if (float(linesplit[2]) <= -cutoff or float(linesplit[3]) <= -cutoff or float(linesplit[5]) <= -cutoff or float(linesplit[6]) <= -cutoff or float(linesplit[2]) >= cutoff or float(linesplit[3]) >= cutoff or float(linesplit[5]) >= cutoff or float(linesplit[6]) >= cutoff):
                    fp_array[count][1] += 1


    ### Sort the list by number of hits in count
    fp_array.sort(key=lambda x: x[1])
    resindex_selected = []
    resindex_remainder = []

    ### Get the index list for the residues with the highest counts and the remainder
    for i in range(max_res):
        resindex_selected.append(fp_array[(num_res-1)-i][0])

    for i in range(num_res - max_res):
        resindex_remainder.append(fp_array[i][0])

    resindex_selected.sort()
    resindex_remainder.sort()
    del fp_array[:][:]

    return resindex_selected, resindex_remainder

####################################################################################################

def plot_footprints(filename, resindex_selected, resindex_remainder):

    ### Read input file of zinc ids
    inputfile = open(filename,'r')
    zinc_ids = inputfile.readlines()
    inputfile.close()

    ### For each zinc_id, read in the zinc_id.txt one by one
    for zinc_id in zinc_ids:
        footprint = open(zinc_id.strip()+'.txt','r')
        lines = footprint.readlines()
        footprint.close()
	
	value_check_bool=False
        ### Store the resname, resid, and fp information appropriately
        resname = []; resid = []; vdw_ref = []; es_ref = []; vdw_pose = []; es_pose = []
        vdw_score = ""; es_score = ""
        vdw_energy = ""; es_energy = ""
        for line in lines:
	    ### When scores are too large & bleed into previous columns
	    if (len(line.split()) > 3) and (len(line.split()) < 8):
               value_check_bool=True
        
        if value_check_bool == True:
            bad_val_output = (open('bad_val_id_list.txt','w'))
            bad_val_output.write(zinc_id)
            bad_val_output.close()
            print("Bad value in ZINC ID:", zinc_id)
            continue

        for line in lines:
            linesplit = line.split()
            if (len(linesplit) == 3):
                if (linesplit[1] == 'desc_FPS_vdw_fps:'):
                    vdw_score = 'd = '+linesplit[2]
                if (linesplit[1] ==  'desc_FPS_es_fps:'):
                    es_score = 'd = '+linesplit[2]
                if (linesplit[1] == 'desc_FPS_vdw_energy:'):
                    vdw_energy = 'vdw = '+linesplit[2]+' kcal/mol'
                if (linesplit[1] == 'desc_FPS_es_energy:'):
                    es_energy = 'es = '+linesplit[2]+' kcal/mol'
            if (len(linesplit) == 8):
                if (linesplit[0] != 'resname'):
                    resname.append(linesplit[0])
                    resid.append(int(linesplit[1]))
                    vdw_ref.append(float(linesplit[2]))
                    es_ref.append(float(linesplit[3]))
                    vdw_pose.append(float(linesplit[5]))
                    es_pose.append(float(linesplit[6]))

        ### Put the selected residues onto a selected array
        resname_selected = []
        vdw_ref_selected = []; es_ref_selected = []; vdw_pose_selected = []; es_pose_selected = []

        for i in (resindex_selected):
            resname_selected.append(resname[i])
            vdw_ref_selected.append(vdw_ref[i])
            es_ref_selected.append(es_ref[i])
            vdw_pose_selected.append(vdw_pose[i])
            es_pose_selected.append(es_pose[i])

        ### Compute the sums for the remainder residues
        vdw_ref_remainder = 0; es_ref_remainder = 0; vdw_pose_remainder = 0; es_pose_remainder = 0

        for i in (resindex_remainder):
            vdw_ref_remainder += vdw_ref[i]
            es_ref_remainder += es_ref[i]
            vdw_pose_remainder += vdw_pose[i]
            es_pose_remainder += es_pose[i]

        ### Append the remainders to the end of the selected arrays
        resname_selected.append('REMAIN')
        vdw_ref_selected.append(vdw_ref_remainder)
        es_ref_selected.append(es_ref_remainder)
        vdw_pose_selected.append(vdw_pose_remainder)
        es_pose_selected.append(es_pose_remainder)
        
        ### Create an index for plotting
        residue = []
        
        for i in range(len(resname_selected)):
            residue.append(i)

        ### Check for self consistency
        #print (sum(vdw_pose_selected) - sum(vdw_pose)) + (sum(vdw_ref_selected) - sum(vdw_ref)) + (sum(es_pose_selected) - sum(es_pose)) + (sum(es_ref_selected) - sum(es_ref))

        ### Plot the figure
        fig = plt.figure(figsize=(12, 11))
        ax1 = fig.add_subplot(2,1,1)
        ax1.set_title(zinc_id.strip())
        plt.plot(residue, vdw_ref_selected, 'b', linewidth=3)
        plt.plot(residue, vdw_pose_selected, 'r', linewidth=3)
        ax1.set_ylabel('VDW Energy')
        ax1.set_ylim(-10, 5)
        ax1.set_xlim(0, len(resname_selected))
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax1.set_xticks(residue)
        ax1.xaxis.grid(which='major', color='black', linestyle='solid')
        ax1.set_xticklabels(resname_selected, rotation=90)
        ax1.legend(['Reference', 'Pose'])
        ax1.annotate(vdw_score, xy=(37,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax1.annotate(vdw_energy, xy=(37,-9), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        
        ax2 = fig.add_subplot(2,1,2)
        plt.plot(residue, es_ref_selected, 'b', linewidth=3)
        plt.plot(residue, es_pose_selected, 'r', linewidth=3)
        ax2.set_ylabel('ES Energy')
        ax2.set_ylim(-10, 5)
        ax2.set_xlim(0, len(resname_selected))
        ax2.xaxis.set_major_locator(MultipleLocator(1))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%s'))
        ax2.set_xticks(residue)
        ax2.xaxis.grid(which='major', color='black', linestyle='solid')
        ax2.set_xticklabels(resname_selected, rotation=90)
        ax2.legend()
        ax2.annotate(es_score, xy=(37,-8), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        ax2.annotate(es_energy, xy=(37,-9), backgroundcolor='white', bbox={'facecolor':'white', 'alpha':1.0, 'pad':10})
        
        plt.savefig(zinc_id.strip()+'.pdf')
        plt.close()

        del resname[:]
        del resid[:]
        del vdw_ref[:]
        del es_ref[:]
        del vdw_pose[:]
        del es_pose[:]
        del resname_selected[:]
        del vdw_ref_selected[:]
        del es_ref_selected[:]
        del vdw_pose_selected[:]
        del es_pose_selected[:]
        del residue[:]

    return

def checkexist(filename):
  if not os.path.isfile(filename):
    print " The file [ %s ] does not exist" % filename
    quit()


####################################################################################################

def main():

    ### Get the command line arguments
    if len(sys.argv)==4:
      filename = sys.argv[1]
      checkexist(filename) # make sure the file exists
      cutoff = float(sys.argv[2])
      max_res = int(sys.argv[3])
    else:
      print "The correct number of argument were not provided"
      print "\t%s [filename] [cutoff] [max_residues]" % sys.argv[0]
      quit()
    ### Go through the first time to identify interactions above the threshold
    (resindex_selected, resindex_remainder) = identify_residues(filename, cutoff, max_res)

    #print resindex_selected; print "\n"; print resindex_remainder
    #print "\n"; print len(resindex_selected); print len(resindex_remainder)

    ### Go through a second time to write plots to file
    plot_footprints(filename, resindex_selected, resindex_remainder)


    return

####################################################################################################

main()

