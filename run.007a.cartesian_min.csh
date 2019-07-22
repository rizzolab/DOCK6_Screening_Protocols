#!/bin/tcsh -fe

#
# After docking all of the ligands from the library to the grid, minimize and rescore each of them
# in Cartesian space.
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
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"



### Choose parameters for cluster
### LIRED    24 ppn
### SeaWulf  28 ppn
### Rizzo    24 ppn


set wcl   = 48:00:00
set nodes = 4
set ppn   = 28
set queue = "long"
@ numprocs = (${nodes} * ${ppn})

### Make the appropriate directory. If it already exists, remove previous dock results from only
### the same vendor.
if (! -e ${rootdir}/${system}/007.cartesian-min) then
	mkdir -p ${rootdir}/${system}/007.cartesian-min/
endif
rm -rf ${rootdir}/${system}/007.cartesian-min/${vendor}
mkdir -p ${rootdir}/${system}/007.cartesian-min/${vendor}
cd ${rootdir}/${system}/007.cartesian-min/${vendor}


### Count the number of chunks
set num_chunks = `ls -l ${rootdir}/${system}/005.zinclibs/${vendor}/ | grep chunk | grep mol2 | wc -l`
echo "num_chunks = ${num_chunks}"
set chunk = "0"
echo "Concatenating all of the chunks together..."

### Iterate over each chunk, concatenate them together
while (${chunk} < ${num_chunks})

	if ( -e ${rootdir}/${system}/006.dock-to-grid/${vendor}/chunk${chunk}/${vendor}.${chunk}.output_scored.mol2 ) then
		cat  ${rootdir}/${system}/006.dock-to-grid/${vendor}/chunk${chunk}/${vendor}.${chunk}.output_scored.mol2 >> input.mol2
	else
		echo "Could not find an output mol2 file in ${rootdir}/${system}/006.dock-to-grid/${vendor}/${chunk}"
	endif

	@ chunk++

end

### Cat on the leftover chunk, if it exists
if ( -e ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/${vendor}.leftover.output_scored.mol2 ) then
	cat  ${rootdir}/${system}/006.dock-to-grid/${vendor}/leftover/${vendor}.leftover.output_scored.mol2 >> input.mol2
endif

### Minimize the reference ligand position in the context of the receptor
echo "Minimizing reference ligand..."

##################################################
cat <<EOF >${system}.${vendor}.reference_minimization.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100
ligand_atom_file                                             ${rootdir}/${system}/001.lig-prep/${system}.lig.python.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ${rootdir}/${system}/001.lig-prep/${system}.lig.python.mol2
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
pharmacophore_score_score_secondary                          no
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
ligand_outfile_prefix                                        output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################

### Execute dock on the headnode
${dockdir}/dock6 -i ${system}.${vendor}.reference_minimization.in -o ${system}.${vendor}.reference_minimization.out
mv output_scored.mol2 ${system}.lig.python.min.mol2

### Write footprints for multigrids later
###########################################################################################
cat <<EOF >${system}.footprint_rescore.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             ${system}.lig.python.min.mol2
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
footprint_similarity_score_primary                           yes
footprint_similarity_score_secondary                         no
fps_score_use_footprint_reference_mol2                       yes
fps_score_footprint_reference_mol2_filename                  ${system}.lig.python.min.mol2
fps_score_foot_compare_type                                  Euclidean
fps_score_normalize_foot                                     no
fps_score_foot_comp_all_residue                              no
fps_score_choose_foot_range_type                             threshold
fps_score_vdw_threshold                                      0.5
fps_score_es_threshold                                       0.2
fps_score_hb_threshold                                       0.5
fps_score_use_remainder                                      yes
fps_score_receptor_filename                                  ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2
fps_score_vdw_att_exp                                        6
fps_score_vdw_rep_exp                                        12
fps_score_vdw_rep_rad_scale                                  1
fps_score_use_distance_dependent_dielectric                  yes
fps_score_dielectric                                         4.0
fps_score_vdw_fp_scale                                       1
fps_score_es_fp_scale                                        1
fps_score_hb_fp_scale                                        0
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              no
atom_model                                                   all
vdw_defn_file                                                ${paramdir}/vdw_AMBER_parm99.defn
flex_defn_file                                               ${paramdir}/flex.defn
flex_drive_file                                              ${paramdir}/flex_drive.tbl
ligand_outfile_prefix                                        output
write_footprints                                             yes
write_hbonds                                                 no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
#################################################
### Create primary residues dat file

${dockdir}/dock6 -i ${system}.footprint_rescore.in -o ${system}.footprint_rescore.out
cp ${system}.footprint_rescore.out footprint_rescore.out

###########################################################################################

grep -A 1 "range_union" footprint_rescore.out | grep -v "range_union" | grep -v "\-" | sed -e '{s/,/\n/g}' | sed -e '{s/ //g}' |  sed '/^$/d' | sort -n | uniq > temp.dat

foreach i ("`cat temp.dat`")  
	printf "%0*d\n" 3 $i >> primary_residues.dat
end

foreach r ("`cat temp.dat`")
        grep " $r " output_footprint_scored.txt  | awk -v temp=$r '{if ($2 == temp) print $0;}' | awk '{print $1 "  " $3 "  " $4}' >> reference.txt
end

grep "remainder" output_footprint_scored.txt | sed -e '{s/,/  /g}' | tr -d '\n' | awk '{print $2 "  " $3 "  " $6}' >> reference.txt

rm temp.dat
cp output_footprint_scored.txt ${system}.footprint.txt
#EOF
###########################################################################################

cp primary_residues.dat ${system}.primary_residues.dat
mv reference.txt ${system}.reference.txt

### Write the dock.in file
##################################################
cat <<EOF >${system}.${vendor}.cartesian_min.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             input.mol2
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
ligand_outfile_prefix                                        ${vendor}.output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
##################################################


### Write the Cluster submit file
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.${vendor}.cartesian_min.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.${vendor}.cartmin
#PBS -V
#PBS -q ${queue}

cd ${rootdir}/${system}/007.cartesian-min/${vendor}/

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.cartesian_min.in \
-o ${system}.${vendor}.cartesian_min.out

EOF

##################################################

### Write the Cluster submit file
##################################################
###if (`hostname -f` == "rizzo.cm.cluster")
else
cat <<EOF >${system}.${vendor}.cartesian_min.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=24
#SBATCH --job-name=${system}.${vendor}.cartmin
#SBATCH --output=${system}.${vendor}.cartmin
#SBATCH -p rn-long

cd ${rootdir}/${system}/007.cartesian-min/${vendor}/

${mpidir}/mpirun -np ${numprocs} \
${dockdir}/dock6.mpi -v \
-i ${system}.${vendor}.cartesian_min.in \
-o ${system}.${vendor}.cartesian_min.out

EOF

endif
##################################################

### Submit the job
echo "Submitting ${system}.${vendor}.cartesian_min "
qsub ${system}.${vendor}.cartesian_min.qsub.csh > & ${system}.${vendor}.cartesian_min.qsub.log
date


exit

