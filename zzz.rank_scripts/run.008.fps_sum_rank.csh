#!/bin/tcsh -fe

#
# Once the ligands have been minimized in Cartesian space, rescore each ligand with a footprint
# reference. The resulting output mol2 files will be transferred back to Cluster and fed into
# the MOE sorting / clustering protocol.
#


### Set some variables manually
set attractive   = "6"
set repulsive    = "12"
set max_num      = "${MAX_NUM_MOL}"

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

### Compile with intel 2013 compiler in case this is not default for some system
### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn

set wcl   = 48:00:00
set nodes = 6
set ppn   = 28
set queue = "long"
@ numprocs = (${nodes} * ${ppn})



### Make the appropriate directory. If it already exists, remove previous dock results from only
### the same vendor.
if (! -e ${rootdir}/${system}/008.rank_mol) then
	mkdir -p ${rootdir}/${system}/008.rank_mol/
endif

if (! -e ${rootdir}/${system}/008.rank_mol/${vendor}) then
       mkdir -p ${rootdir}/${system}/008.rank_mol/${vendor}/
endif

rm -rf ${rootdir}/${system}/008.rank_mol/${vendor}/fps_sum_rank
mkdir -p ${rootdir}/${system}/008.rank_mol/${vendor}/fps_sum_rank
cd ${rootdir}/${system}/008.rank_mol/${vendor}/fps_sum_rank


### Compute descriptor scores for the minimized poses
echo "Ranking scores..."


### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.fps_sum_rank_score.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${rootdir}/${system}/007.cartesian-min/${vendor}/${vendor}.output_scored.mol2
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
descriptor_use_footprint_similarity                          yes
descriptor_use_pharmacophore_score                           no
descriptor_use_tanimoto                                      no
descriptor_use_hungarian                                     no
descriptor_use_volume_overlap                                no
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
descriptor_weight_fps_score                                  1
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        ${vendor}.fps_sum.output
write_footprints                                             no
write_hbonds                                                 no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 yes
max_ranked_ligands                                           ${max_num}
EOF
##################################################


### Write the Cluster submit file to maui
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.${vendor}.fps_sum_rank.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.vs.fps_sum_rank
#PBS -q ${queue}
#PBS -V
#PBS -j oe


cd ${rootdir}/${system}/008.rank_mol/${vendor}/fps_sum_rank

echo "Job  Started"
date

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.fps_sum_rank_score.in \
-o ${system}.${vendor}.fps_sum_rank_score.out

echo "Job  Finished"
date

EOF

##################################################

### Write the Cluster submit file to slurm
##################################################
###if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
else
cat <<EOF >${system}.${vendor}.fps_sum_rank.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=24
#SBATCH --job-name=${system}.vs.fps_sum_rank
#SBATCH --output=${system}.vs.fps_sum_rank
#SBATCH -p rn-long


cd ${rootdir}/${system}/008.rank_mol/${vendor}/fps_sum_rank

echo "Job  Started"
date

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.fps_sum_rank_score.in \
-o ${system}.${vendor}.fps_sum_rank_score.out

echo "Job  Finished"
date

EOF
endif
##################################################


### Submit the job
echo "Submitting ${system}.${vendor}.fps_sum_rank.qsub.csh"
qsub ${system}.${vendor}.fps_sum_rank.qsub.csh > & ${system}.${vendor}.fps_sum_rank.qsub.log
date


exit

