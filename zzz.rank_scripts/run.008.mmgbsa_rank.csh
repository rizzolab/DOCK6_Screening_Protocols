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
set nodes = 8
set ppn   = 24
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

rm -rf ${rootdir}/${system}/008.rank_mol/${vendor}/mmgbsa_rank
mkdir -p ${rootdir}/${system}/008.rank_mol/${vendor}/mmgbsa_rank
cd ${rootdir}/${system}/008.rank_mol/${vendor}/mmgbsa_rank


### Compute descriptor scores for the minimized poses
echo "Ranking scores..."


### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.mmgbsa_rank_score.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100
ligand_atom_file                                             ${rootdir}/${system}/007.cartesian-min/${vendor}/${vendor}.output_scored.mol2
limit_max_ligands                                            no
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
descriptor_score_primary                                     no
descriptor_score_secondary                                   no
gbsa_zou_score_primary                                       no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_primary                                   yes
gbsa_hawkins_score_secondary                                 no
gbsa_hawkins_score_rec_filename                              ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
gbsa_hawkins_score_solvent_dielectric                        78.5
gbsa_hawkins_use_salt_screen                                 no
gbsa_hawkins_score_gb_offset                                 0.09
gbsa_hawkins_score_cont_vdw_and_es                           yes
gbsa_hawkins_score_vdw_att_exp                               ${attractive}
gbsa_hawkins_score_vdw_rep_exp                               ${repulsive}
grid_score_rep_rad_scale                                     1
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 yes
max_ranked_ligands                                           ${max_num}
EOF



### Write the Cluster submit file
##################################################
if (`hostname -f` == "lired.cm.cluster" || `hostname -f` == "login.cm.cluster" ) then
cat <<EOF >${system}.${vendor}.mmgbsa_rank.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.${vendor}.mmgbsa_rank
#PBS -q ${queue}
#PBS -V
#PBS -j oe


cd ${rootdir}/${system}/008.rank_mol/${vendor}/mmgbsa_rank

echo "Job  Started"
date

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.mmgbsa_rank_score.in \
-o ${system}.${vendor}.mmgbsa_rank_score.out

echo "Job  Finished"
date

EOF
endif
##################################################



### Submit the job
echo "Submitting ${system}.${vendor}.mmgbsa_rank.qsub.csh"
qsub ${system}.${vendor}.mmgbsa_rank.qsub.csh > & ${system}.${vendor}.mmgbsa_rank.qsub.log
date


exit

