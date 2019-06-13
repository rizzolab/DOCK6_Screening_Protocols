#!/bin/tcsh -fe

#
# Once the ligands have been minimized in Cartesian space, rescore each ligand with a footprint
# reference. The resulting output mol2 files will be transferred back to Cluster and fed into
# the MOE sorting / clustering protocol.
#


### Set some variables manually
set attractive   = "6"
set repulsive    = "12"


### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set mpidir    = "${VS_MPIDIR}/bin"

set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set descriptdir = "${rootdir}/zzz.descriptor"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"

### Compile with intel 2018.3 compiler in case this is not default for some system
### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn

#-------!!!!!!!!!!! DO NOT SUBMIT THIS IN PARALLEL it will mess up the rank order mol2-----#
# This scripts ask for the entire node but ONLY USES one proc
# That is why we used dock6 and dock6.mpi

echo "\t ATTN: This script SHOULD not be run PARALLEL "

set wcl   = 12:00:00
set nodes = 1
set ppn   = 24
set queue = "long-24core"
@ numprocs = (${nodes} * ${ppn})



### Make the appropriate directory. If it already exists, remove previous dock results from only
### the same vendor.
if (! -e ${rootdir}/${system}/009.descriptor-rescore) then
	mkdir -p ${rootdir}/${system}/009.descriptor-rescore/
endif

if (! -e ${rootdir}/${system}/009.descriptor-rescore/${vendor}) then
	mkdir -p ${rootdir}/${system}/009.descriptor-rescore/${vendor}
endif

### Compute descriptor scores for the minimized crystal pose
echo "Computing Crystal pose descriptor scores..."

if (! -e ${rootdir}/${system}/009.descriptor-rescore/${vendor}/reference) then
        mkdir -p ${rootdir}/${system}/009.descriptor-rescore/${vendor}/reference
endif

cd ${rootdir}/${system}/009.descriptor-rescore/${vendor}/reference

### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.reference.descriptor_rescore.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
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
descriptor_use_continuous_score                              yes
descriptor_use_footprint_similarity                          yes
descriptor_use_pharmacophore_score                           yes
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_use_volume_overlap                                yes
descriptor_cont_score_rec_filename                           ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
descriptor_cont_score_att_exp                                6
descriptor_cont_score_rep_exp                                12
descriptor_cont_score_rep_rad_scale                          1
descriptor_cont_score_use_dist_dep_dielectric                yes
descriptor_cont_score_dielectric                             4.0
descriptor_cont_score_vdw_scale                              1
descriptor_cont_score_es_scale                               1
descriptor_fps_score_use_footprint_reference_mol2            yes
descriptor_fps_score_footprint_reference_mol2_filename       ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_fps_score_foot_compare_type                       Euclidean
descriptor_fps_score_normalize_foot                          no
descriptor_fps_score_foot_comp_all_residue                   yes
descriptor_fps_score_receptor_filename                       ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
descriptor_fps_score_vdw_att_exp                             6
descriptor_fps_score_vdw_rep_exp                             12
descriptor_fps_score_vdw_rep_rad_scale                       1
descriptor_fps_score_use_distance_dependent_dielectric       yes
descriptor_fps_score_dielectric                              4.0
descriptor_fps_score_vdw_fp_scale                            1
descriptor_fps_score_es_fp_scale                             1
descriptor_fps_score_hb_fp_scale                             0
descriptor_fms_score_use_ref_mol2                            yes
descriptor_fms_score_ref_mol2_filename                       ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_fms_score_write_reference_pharmacophore_mol2      no
descriptor_fms_score_write_reference_pharmacophore_txt       no
descriptor_fms_score_write_candidate_pharmacophore           no
descriptor_fms_score_write_matched_pharmacophore             no
descriptor_fms_score_compare_type                            overlap
descriptor_fms_score_full_match                              yes
descriptor_fms_score_match_rate_weight                       5.0
descriptor_fms_score_match_dist_cutoff                       1.0
descriptor_fms_score_match_proj_cutoff                       0.7071
descriptor_fms_score_max_score                               20
descriptor_fingerprint_ref_filename                          ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_hms_score_ref_filename                            ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_volume_score_reference_mol2_filename              ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_volume_score_overlap_compute_method               analytical
descriptor_weight_cont_score                                 1
descriptor_weight_fps_score                                  1
descriptor_weight_pharmacophore_score                        2
descriptor_weight_fingerprint_tanimoto                       -1
descriptor_weight_hms_score                                  10
descriptor_weight_volume_overlap_score                       -20
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
chem_defn_file                                               ${paramdir}/chem.defn
pharmacophore_defn_file                                      ${paramdir}/ph4.defn
ligand_outfile_prefix                                        ${system}.reference.output
write_footprints                                             yes
write_hbonds                                                 yes
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################

${dockdir}/dock6 -i ${system}.${vendor}.reference.descriptor_rescore.in -o ${system}.${vendor}.reference.descriptor_rescore.out

cd ../

foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
  echo ${score}

  rm -rf ${rootdir}/${system}/009.descriptor-rescore/${vendor}/${score}_rank
  mkdir -p ${rootdir}/${system}/009.descriptor-rescore/${vendor}/${score}_rank
  cd ${rootdir}/${system}/009.descriptor-rescore/${vendor}/${score}_rank



### Compute descriptor scores for the minimized poses
echo "Computing ${score} ranked descriptor scores..."


### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.${score}.descriptor_rescore.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${rootdir}/${system}/008.rank_mol/${vendor}/${score}_rank/${vendor}.${score}.output_ranked.mol2
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
descriptor_use_continuous_score                              yes
descriptor_use_footprint_similarity                          yes
descriptor_use_pharmacophore_score                           yes
descriptor_use_tanimoto                                      yes
descriptor_use_hungarian                                     yes
descriptor_use_volume_overlap                                yes
descriptor_cont_score_rec_filename                           ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
descriptor_cont_score_att_exp                                6
descriptor_cont_score_rep_exp                                12
descriptor_cont_score_rep_rad_scale                          1
descriptor_cont_score_use_dist_dep_dielectric                yes
descriptor_cont_score_dielectric                             4.0
descriptor_cont_score_vdw_scale                              1
descriptor_cont_score_es_scale                               1
descriptor_fps_score_use_footprint_reference_mol2            yes
descriptor_fps_score_footprint_reference_mol2_filename       ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_fps_score_foot_compare_type                       Euclidean
descriptor_fps_score_normalize_foot                          no
descriptor_fps_score_foot_comp_all_residue                   yes
descriptor_fps_score_receptor_filename                       ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
descriptor_fps_score_vdw_att_exp                             6
descriptor_fps_score_vdw_rep_exp                             12
descriptor_fps_score_vdw_rep_rad_scale                       1
descriptor_fps_score_use_distance_dependent_dielectric       yes
descriptor_fps_score_dielectric                              4.0
descriptor_fps_score_vdw_fp_scale                            1
descriptor_fps_score_es_fp_scale                             1
descriptor_fps_score_hb_fp_scale                             0                      
descriptor_fms_score_use_ref_mol2                            yes
descriptor_fms_score_ref_mol2_filename                       ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_fms_score_write_reference_pharmacophore_mol2      no
descriptor_fms_score_write_reference_pharmacophore_txt       no
descriptor_fms_score_write_candidate_pharmacophore           no
descriptor_fms_score_write_matched_pharmacophore             no
descriptor_fms_score_compare_type                            overlap
descriptor_fms_score_full_match                              yes
descriptor_fms_score_match_rate_weight                       5.0
descriptor_fms_score_match_dist_cutoff                       1.0
descriptor_fms_score_match_proj_cutoff                       0.7071
descriptor_fms_score_max_score                               20
descriptor_fingerprint_ref_filename                          ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_hms_score_ref_filename                            ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_hms_score_matching_coeff                          -5
descriptor_hms_score_rmsd_coeff                              1
descriptor_volume_score_reference_mol2_filename              ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2
descriptor_volume_score_overlap_compute_method               analytical
descriptor_weight_cont_score                                 1
descriptor_weight_fps_score                                  1
descriptor_weight_pharmacophore_score                        2
descriptor_weight_fingerprint_tanimoto                       -1
descriptor_weight_hms_score                                  10
descriptor_weight_volume_overlap_score                       -20
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
chem_defn_file                                               ${paramdir}/chem.defn
pharmacophore_defn_file                                      ${paramdir}/ph4.defn
ligand_outfile_prefix                                        ${vendor}.${score}_rank.output
write_footprints                                             yes
write_hbonds                                                 yes
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################


### Write the Cluster submit file to maui
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.${vendor}.${score}.descriptor_rescore.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${score}.${system}.${vendor}.descriptor
#PBS -q ${queue}
#PBS -V
#PBS -j oe


cd ${rootdir}/${system}/009.descriptor-rescore/${vendor}/${score}_rank

echo "Job  Started"
date

${dockdir}/dock6 -v \
-i ${system}.${vendor}.${score}.descriptor_rescore.in \
-o ${system}.${vendor}.${score}.descriptor_rescore.out

echo "Job  Finished"
date

EOF

##################################################

### Write the Cluster submit file to slurm
##################################################
###if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
else
cat <<EOF >${system}.${vendor}.${score}.descriptor_rescore.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${score}.${system}.${vendor}.descriptor
#SBATCH --output=${score}.${system}.${vendor}.descriptor
#SBATCH -p rn-long


cd ${rootdir}/${system}/009.descriptor-rescore/${vendor}/${score}_rank

echo "Job  Started"
date

${dockdir}/dock6 -v \
-i ${system}.${vendor}.${score}.descriptor_rescore.in \
-o ${system}.${vendor}.${score}.descriptor_rescore.out

echo "Job  Finished"
date

EOF
endif
##################################################


### Submit the job
echo "Submitting ${system}.${vendor}.${score}.descriptor_rescore"
qsub ${system}.${vendor}.${score}.descriptor_rescore.qsub.csh > & ${system}.${vendor}.${score}.descriptor_rescore.qsub.log
date

end
exit

