#!/bin/tcsh -fe

### This scripts submits a de novo job with multiple anchors, which is user-defined 

### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}/bin"
set fragdir   = "${FRAGLIB}"
set ancdir    = "${ANCLIB}"

set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"

### Declare some denovo variables 
### max_anchors can be changed to however many anchors you want to seed dn growth with
set max_anchors = "30"          #default is 30 anchors  
set dn_heur_unmatched_num="1"
set dn_heur_matched_rmsd="2.0"
set dn_max_grow_layers="9"
set dn_max_current_aps="5"
set dn_max_root_size="25"	#standard sampling is 25; enhanced is 100
set dn_max_layer_size="25"	#standard sampling is 25; enhanced is 100
set simplex_restraint_min="no"
set dn_prun_cut="100.0"
set dn_prun_scale="2.0"
set dn_fc="2.0"          	#default is 2.0
set max_scaf_lyr="1"

### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn

set wcl   = 120:00:00
set nodes = 1
set ppn   = 24
set queue = "extended"
@ numprocs = (${nodes} * ${ppn})

### Make the appropriate directory. If it already exists, remove previous dock results from only
### the same vendor.
if (! -e ${rootdir}/${system}/012.denovo) then
        mkdir -p ${rootdir}/${system}/012.denovo/
endif
rm -rf ${rootdir}/${system}/012.denovo/${vendor}
mkdir -p ${rootdir}/${system}/012.denovo/${vendor}
cd ${rootdir}/${system}/012.denovo/${vendor}

### Loop through each anchor(a) to create input files, populate multigrids
### and create PBS scripts to submit. 
set a = "1"
while ($a <= ${max_anchors})
	
	cd ${rootdir}/${system}/012.denovo/${vendor}
	rm -fr anchor_$a
	mkdir anchor_$a
	cd anchor_$a
	
	set num_of_grids = `wc -l ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.primary_residues.dat | awk '{print $1+1}'`

###############################################################################################################
cat << eof > denovo.anchor_${a}.in
conformer_search_type                                        denovo
dn_fraglib_scaffold_file                                     ${fragdir}/fraglib_scaffold.mol2
dn_fraglib_linker_file                                       ${fragdir}/fraglib_linker.mol2
dn_fraglib_sidechain_file                                    ${fragdir}/fraglib_sidechain.mol2
dn_fraglib_rigid_file                                        ${fragdir}/fraglib_rigid.mol2
dn_user_specified_anchor                                     yes
dn_fraglib_anchor_file                                       ${ancdir}/anchor_${a}.mol2 
dn_use_torenv_table                                          yes
dn_torenv_table                                              ${fragdir}/fraglib_torenv.dat
dn_sampling_method                                           graph
dn_graph_max_picks                                           30
dn_graph_breadth                                             3
dn_graph_depth                                               2
dn_graph_temperature                                         100
dn_pruning_conformer_score_cutoff                            ${dn_prun_cut}
dn_pruning_conformer_score_scaling_factor                    ${dn_prun_scale}
dn_pruning_clustering_cutoff                                 100.0
dn_constraint_mol_wt                                         750
dn_constraint_rot_bon                                        15
dn_constraint_formal_charge                                  ${dn_fc}
dn_heur_unmatched_num                                        ${dn_heur_unmatched_num}
dn_heur_matched_rmsd                                         ${dn_heur_matched_rmsd}
dn_unique_anchors                                            3
dn_max_grow_layers                                           ${dn_max_grow_layers}
dn_max_root_size                                             ${dn_max_root_size}
dn_max_layer_size                                            ${dn_max_layer_size}
dn_max_current_aps                                           ${dn_max_current_aps}
dn_max_scaffolds_per_layer                                   ${max_scaf_lyr}
dn_write_checkpoints                                         yes
dn_write_prune_dump                                          yes
dn_write_orients                                             no
dn_write_growth_trees                                        no
dn_output_prefix                                             ${system}.final
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ${rootdir}/${system}/003.spheres/${system}.rec.clust.close.sph
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
gist_score_primary                                           no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           no
footprint_similarity_score_secondary                         no
descriptor_score_primary                                     yes
descriptor_score_secondary                                   no
descriptor_use_grid_score                                    no
descriptor_use_multigrid_score                               yes
descriptor_use_tanimoto                                      no
descriptor_use_hungarian                                     no
descriptor_multigrid_score_rep_rad_scale                     1
descriptor_multigrid_score_vdw_scale                         1
descriptor_multigrid_score_es_scale                          1
descriptor_multigrid_score_number_of_grids                   ${num_of_grids}
eof

set i = "0"
while ($i < ${num_of_grids})

set file = `ls -l ${rootdir}/${system}/004.grid | grep resid | grep nrg | awk -v c="$i" '{if (NR == (c+1))  print $9}' | sed 's/.nrg//'`
cat << eof >> denovo.anchor_${a}.in
descriptor_multigrid_score_grid_prefix${i}                   ${rootdir}/${system}/004.grid/${file}
eof

@ i++
end
cat << eof >> denovo.anchor_${a}.in
descriptor_multigrid_score_fp_ref_mol                        yes
descriptor_multigrid_score_footprint_ref                     ${rootdir}/${system}/004.grid/${system}.lig.multigridmin.mol2
descriptor_multigrid_score_use_euc                           yes
descriptor_multigrid_score_use_norm_euc                      no
descriptor_multigrid_score_use_cor                           no
descriptor_multigrid_vdw_euc_scale                           1
descriptor_multigrid_es_euc_scale                            1
descriptor_weight_multigrid_score                            1
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        ${simplex_restraint_min}
simplex_coefficient_restraint                                10.0
atom_model                                                   all
vdw_defn_file                                                ${DOCKHOMEWORK}/parameters/vdw_de_novo.defn
flex_defn_file                                               ${DOCKHOMEWORK}/parameters/flex.defn
flex_drive_file                                              ${DOCKHOMEWORK}/parameters/flex_drive.tbl
eof

###############################################################################################################
### Rescore each compound with Tanimoto to prepare for make_unique 
### in the next step. 
###############################################################################################################

cat <<EOF >rescore.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ${system}.final.denovo_build.mol2
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
gist_score_primary                                           no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     no
continuous_score_secondary                                   no
footprint_similarity_score_primary                           no
footprint_similarity_score_secondary                         no
ph4_score_primary                                            no
ph4_score_secondary                                          no
descriptor_score_primary                                     yes
descriptor_score_secondary                                   no
descriptor_use_grid_score                                    no
descriptor_use_multigrid_score                               yes
descriptor_use_continuous_energy                             no
descriptor_use_footprint_similarity                          no
descriptor_use_pharmacophore_score                           no
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_multigrid_score_rep_rad_scale                     1
descriptor_multigrid_score_vdw_scale                         1
descriptor_multigrid_score_es_scale                          1
descriptor_multigrid_score_number_of_grids                   ${num_of_grids}
EOF

set i = "0"

while ($i < ${num_of_grids})

set file = `ls -l ${rootdir}/${system}/004.grid | grep resid | grep nrg | awk -v c="$i" '{if (NR == (c+1))  print $9}' | sed 's/.nrg//'`
cat <<EOF >> rescore.in
descriptor_multigrid_score_grid_prefix${i}                   ${rootdir}/${system}/004.grid/${file}
EOF

@ i++

end

cat <<EOF >> rescore.in
descriptor_multigrid_score_fp_ref_mol                        yes
descriptor_multigrid_score_footprint_ref                     ${rootdir}/${system}/004.grid/${system}.lig.multigridmin.mol2
descriptor_multigrid_score_foot_compare_type                 Euclidean
descriptor_multigrid_score_normalize_foot                    no
descriptor_multigrid_score_use_cor                           no
descriptor_multigrid_vdw_euc_scale                           1
descriptor_multigrid_es_euc_scale                            1
descriptor_weight_multigrid_score                            1
descriptor_fingerprint_ref_filename                          ${rootdir}/${system}/004.grid/${system}.lig.multigridmin.mol2
descriptor_hms_score_ref_filename                            ${rootdir}/${system}/004.grid/${system}.lig.multigridmin.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_weight_fingerprint_tanimoto                       1
descriptor_weight_hms_score                                  1
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${DOCKHOMEWORK}/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               ${DOCKHOMEWORK}/parameters/flex.defn
flex_drive_file                                              ${DOCKHOMEWORK}/parameters/flex_drive.tbl
ligand_outfile_prefix                                        rescore
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF

###############################################################################################################
### Create PBS script for each anchor and submit
###############################################################################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat << eof > run.anchor${a}.qsub.sh
#/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N anchor${a}
#PBS -V
#PBS -q ${queue}


echo "Job started on"
date
cd ${rootdir}/${system}/012.denovo/${vendor}/anchor_${a}

${dockdir}/dock6 -i denovo.anchor_${a}.in -o denovo.anchor_${a}.out
${dockdir}/dock6 -i rescore.in -o rescore.out

echo "Job finished on"
date
eof
###############################################################################################################

### Write the Cluster submit file to slurm
##################################################
###if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
else
cat << EOF > run.anchor${a}.qsub.sh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${score}.${system}.${vendor}.descriptor
#SBATCH --output=${score}.${system}.${vendor}.descriptor
#SBATCH -p ${queue}


echo "Job started on"
date
cd ${rootdir}/${system}/012.denovo/${vendor}/anchor_${a}

${dockdir}/dock6 -i denovo.anchor_${a}.in -o denovo.anchor_${a}.out
${dockdir}/dock6 -i rescore.in -o rescore.out

echo "Job finished on"
date

EOF
endif
##################################################




### Consider submitting individually so as not to overwhelm the queue
#qsub run.anchor${a}.qsub.sh

	@ a++
	echo $a
end
