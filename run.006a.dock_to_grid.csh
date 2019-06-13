#!/bin/tcsh -fe

#
# Run this script on Cluster when you are ready to dock to the grid. Assumes receptor / grids
# and ligand library have already been prepared.
#


### Set some paths
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

### Compile with intel 2013 compiler in case this is not default for some system


### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn

set wcl   = 48:00:00
set nodes = 8
set ppn   = 28
set queue = "long"
@ numprocs = (${nodes} * ${ppn})


### Make the appropriate directory. If it already exists, remove previous dock results from only
### the same vendor.
if (! -e ${rootdir}/${system}/006.dock-to-grid) then
	mkdir -p ${rootdir}/${system}/006.dock-to-grid/
endif
rm -rf ${rootdir}/${system}/006.dock-to-grid/${vendor}
mkdir -p ${rootdir}/${system}/006.dock-to-grid/${vendor}
cd ${rootdir}/${system}/006.dock-to-grid/${vendor}


### Count the number of chunks
set num_chunks = `ls -l ${rootdir}/${system}/005.zinclibs/${vendor}/ | grep chunk | grep mol2 | wc -l`
echo "num_chunks = ${num_chunks}"
set chunk = "0"


### Iterate over each chunk
while (${chunk} < ${num_chunks})

	mkdir chunk${chunk}/
	cd chunk${chunk}/
	echo "Docking chunk named chunk${chunk}_scored.mol2"


### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.${chunk}.dock_to_grid.in
conformer_search_type                                        flex
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${rootdir}/${system}/005.zinclibs/${vendor}/chunk${chunk}_scored.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
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
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ${rootdir}/${system}/004.grid/${system}.rec
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
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
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        ${vendor}.${chunk}.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################


### Write the Cluster submit file for Maui
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.${vendor}.${chunk}.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.${vendor}.${chunk}
#PBS -V
#PBS -q ${queue}

echo "Job started on"
date
cd ${rootdir}/${system}/006.dock-to-grid/${vendor}/chunk${chunk}/

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.${chunk}.dock_to_grid.in \
-o ${system}.${vendor}.${chunk}.dock_to_grid.out
echo "Job finished on"
date

EOF

##################################################

### Write the Cluster submit file for slurm
##################################################
###(`hostname -f` == "rizzo.cm.cluster" ) 
else
cat <<EOF >${system}.${vendor}.${chunk}.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${system}.${vendor}.${chunk}
#SBATCH --output=${system}.${vendor}.${chunk}.out
#SBATCH -p ${queue}

echo "Job started on"
date
cd ${rootdir}/${system}/006.dock-to-grid/${vendor}/chunk${chunk}/

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.${chunk}.dock_to_grid.in \
-o ${system}.${vendor}.${chunk}.dock_to_grid.out
echo "Job finished on"
date

EOF

endif
##################################################


	### Submit the job
	echo "Submitting ${system}.${vendor}.${chunk}.dock_to_grid " >> ../zzz.submit.log
	qsub ${system}.${vendor}.${chunk}.qsub.csh > & ${system}.${vendor}.${chunk}.qsub.log
        sleep 2 # so the queue is not overwhelmed
	date >> ../zzz.submit.log

	@ chunk++
	cd ../
end


exit

