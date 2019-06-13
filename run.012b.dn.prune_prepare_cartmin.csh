#!/bin/tcsh -fe

### This script takes a file of denovo results that have already been scored with tanimoto 
### and calls DOCK iteratively to create a new file with only the top scoring copy of each 
### molecule in a separate file. 
### This script requires a multimol2 file containing tanimoto results, and split_on_tanimoto.py
### WARNING removes all files starting with "unique_"  and "temp" so please run once per directory

set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}/bin"
set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"
set fragdir   = "${FRAGLIB}"
set ancdir    = "${ANCLIB}"

### Set some variables manually
set attractive   = "6"
set repulsive    = "12"

### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn


set wcl   = 48:00:00
set nodes = 4
set ppn   = 28
set queue = "long"
@ numprocs = (${nodes} * ${ppn})

### name of file you want to make unique, output is called unique_${FIL}

rm -f unique*

### this assumes that $FIL already has been scored with Tanimoto to an arbitrary reference
python ${scriptdir}/split_on_tanimoto_new.py ${system}.final.denovo_build.mol2

set FNUM = "0"

### now we go through all these individual files one by one
### if you wanted to make the script faster by calling the python script 
### after each rescoring this logic would still work
while ( -e temp.split.${FNUM}.mol2 )

	### this is used as a checklist of molecules in each file to keep track of which we have 
	### copied a conformer of
	set USED=()

### make DOCK input file for rescoring
#####################################################################################################################
#####################################################################################################################
cat <<EOFB > temp_rescore.${FNUM}.in

conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             temp.split.${FNUM}.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           no
footprint_similarity_score_secondary                         no
pharmacophore_score_primary                                  no
pharmacophore_score_secondary                                no
descriptor_score_primary                                     yes
descriptor_score_secondary                                   no
descriptor_use_grid_score                                    no
descriptor_use_multigrid_score                               no
descriptor_use_continuous_score                              no
descriptor_use_footprint_similarity                          no
descriptor_use_pharmacophore_score                           no
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     no
descriptor_use_volume_overlap                                no
descriptor_fingerprint_ref_filename                          ref.mol2
descriptor_weight_fingerprint_tanimoto                       1
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${DOCKHOMEWORK}/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               ${DOCKHOMEWORK}/parameters/flex.defn
flex_drive_file                                              ${DOCKHOMEWORK}/parameters/flex_drive.tbl
ligand_outfile_prefix                                        temp
write_orientations                                           no

EOFB
#####################################################################################################################
#####################################################################################################################

	@ MOLS ="`grep --count MOL temp.split.${FNUM}.mol2`"
	set COUNT = "1"

	### if theres only one mol we know its unique so copy it over
	if ( ${MOLS} == 1 ) then
		cat temp.split.${FNUM}.mol2 >> unique_${system}.final.denovo_build.mol2
		@ FNUM++
		continue 
	endif

###########################################################################################

	while ( ${COUNT} <= ${MOLS} ) # we go through each molecule in the file
		set BOOL = "true"                
		### check to see if we have "USED" the molecule
		foreach ITEM ($USED)
			if ( $COUNT ==  $ITEM ) then
				echo $COUNT
				echo $ITEM
				set BOOL = "false"
				echo $BOOL
			else
				continue
			endif
		end #for ITEM in USED
   		echo $BOOL
		if ($BOOL == true) then
		echo $BOOL

		### take the "COUNT" molecule of file and use it as new reference
		### so on the first iteration we use the first molecule, on the second we use the second...
			cat temp.split.${FNUM}.mol2 | awk -v count=${COUNT} 'BEGIN{num = 0} /Name/{num += 1} {if (num == count) print}' > ref.mol2

		### run dock to get tanimoto of all molecules in file to the new reference
			${dockdir}/dock6 -i temp_rescore.${FNUM}.in -o temp_rescore.${FNUM}.out

		### we make a list containing the multigrid score and tanimoto (to new reference) of each mol
			grep "MultiGrid_Score:" temp.split.${FNUM}.mol2 > temp_scores.txt
			grep "Tanimoto_Score:" temp_scored.mol2 >  temp_tan.txt
			paste temp_scores.txt temp_tan.txt > temp_merge.txt

		### find conformer with best score that has a tanimoto of 1 to the new ref
			set BEST=` cat temp_merge.txt | awk 'BEGIN{score = 9999; position = 0} {if (($3 < score)&&($6==1)) score = $3} {if (($3 == score)&&($6==1)) position = NR} END{print position}' `
                	echo $BEST

		### mark all the molecules that match new ref so we don't consider them again
			foreach VAL (` cat temp_merge.txt | awk '{if ($6==1) print NR}'`) 
				echo $VAL
				set USED=($USED $VAL)
			end #for VAL in USED

		### assuming your molecule entries start with 
		####### Name: copies molecule number $BEST into the output file
		### we copy from temp.split.${FNUM} to preserve original score data
			cat temp.split.${FNUM}.mol2 | awk -v best=${BEST} 'BEGIN{num = 0} /Name/{num += 1} {if (num == best) print}' >> unique_${system}.final.denovo_build.mol2

		endif #if we haven't already USED the molecule

	@ COUNT++ 

	end  # while $COUNT
#unset USED
@ FNUM++ 
end  # while FN
unset MOLS
rm -f temp* ref.mol2

####################################################################################################################
### Done with make unique - now combine, rename, and prepare the molecules. 

cd ${rootdir}/${system}/012.denovo/${vendor}
rm -rf anchors_all

set max_anchors = `ls -l ${rootdir}/${system}/012.denovo/${vendor}/ | grep anchor | wc -l`
echo "max_anchors = ${max_anchors}"

set a = "1"
while ($a <= ${max_anchors})

        cd ${rootdir}/${system}/012.denovo/${vendor}
        mkdir -p anchors_all
        cd anchors_all
        touch combined_denovo_anchors.mol2
        cd ../
        cd anchor_$a
        set date = `date +%Y%m%d_%H%M%S`
        ### The denovo name will be in this format:
        ### YearMonthDay_HourMinuteSecond_anchor#_rank#
        set denovo_name = ${date}'_a'${a}
        echo ${denovo_name}

### De novo molecules must be given a unique identifier
        python ${scriptdir}/replace_denovo_names.py ${system}.final.denovo_build.mol2 ${system}.renamed.denovo_build.mol2 ${denovo_name}

### The residues must be renumbered in column 7 for it to be compatible with Moe
        python ${scriptdir}/renumber_residues_in_mol2_for_moe.py ${system}.renamed.denovo_build.mol2 ${system}.renamed.renumbered.denovo_build.mol2

### Concatenate all anchor denovo files together to begin minimizing, rescoring, and clustering.
        cat ${system}.renamed.renumbered.denovo_build.mol2 >> ../anchors_all/combined_denovo_anchors.mol2
echo "Concatenating all of the anchors together..."

        @ a++
end

cd ${rootdir}/${system}/012.denovo/${vendor}/anchors_all

echo "Minmizing de novo molecules ..."
### Write the dock.in file
##################################################
cat <<EOF >${system}.denovo.cartesian_min.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             combined_denovo_anchors.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     yes
continuous_score_secondary                                   no
cont_score_rec_filename                                      ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
cont_score_att_exp                                           ${attractive}
cont_score_rep_exp                                           ${repulsive}
cont_score_rep_rad_scale                                     1
cont_score_use_dist_dep_dielectric                           yes
cont_score_dielectric                                        4.0
cont_score_vdw_scale                                         1
cont_score_es_scale                                          1
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        yes
simplex_coefficient_restraint                                5.0
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        ${system}.denovo.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################


### Write the Cluster submit file to maui
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.denovo.cartesian_min.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.denovo.cartmin
#PBS -V
#PBS -q ${queue}

cd ${rootdir}/${system}/012.denovo/${vendor}/anchors_all

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.denovo.cartesian_min.in \
-o ${system}.denovo.cartesian_min.out

EOF

##################################################
### Write the Cluster submit file to slurm
##################################################
###if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
else
cat <<EOF >${system}.denovo.cartesian_min.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${system}.denovo.cartesian_min.qsub.csh
#SBATCH --output=${system}.denovo.cartesian_min.qsub.csh
#SBATCH -p rn-long

cd ${rootdir}/${system}/012.denovo/${vendor}/anchors_all

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.denovo.cartesian_min.in \
-o ${system}.denovo.cartesian_min.out

EOF
endif
###################################################
## Submit the job
echo "Submitting ${system}.denovo.cartesian_min "
qsub ${system}.denovo.cartesian_min.qsub.csh > & ${system}.denovo.cartesian_min.qsub.log
date


exit

