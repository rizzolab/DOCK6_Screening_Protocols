#!/bin/tcsh -fe
#While the cartmin job is on the queue or running:
#This script will create multigrids based off of the primary residue dat file created in 007.
#The reference is first minimized in multigrid space and then multigrids are created. 

### Set some variables manually
set attractive   = "6"
set repulsive    = "9"

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
set ppn   = 24
set queue = "long-24core"
@ numprocs = (${nodes} * ${ppn})

### cd into the appriopriate directory
cd ${rootdir}/${system}/004.grid


echo "Minimizing reference ligand in gridspace..."
### DOCK input file for minimizing on grid
###########################################################################################
cat <<EOF >${system}.reference_gridmin.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
ligand_atom_file                                             ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ${system}.rec
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
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
ligand_outfile_prefix                                        grid_output
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF
###########################################################################################

### Execute dock on the headnode
${dockdir}/dock6 -i ${system}.reference_gridmin.in -o ${system}.reference_gridmin.out
mv grid_output_scored.mol2 ${system}.lig.gridmin.mol2

### Grid input file for multi-grid
###########################################################################################
cat <<EOF >${system}.multigrid.in
compute_grids                  yes
grid_spacing                   0.4
output_molecule                yes
contact_score                  no
chemical_score                 no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            ${attractive}
repulsive_exponent             ${repulsive}
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  temp.mol2
box_file                       ./box.pdb
vdw_definition_file            ./vdw.defn
chemical_definition_file       ./chem.defn
score_grid_prefix              temp.rec
receptor_out_file              temp.rec.grid.mol2
EOF
###########################################################################################

setenv num_of_grids `wc -l ../007.cartesian-min/${vendor}/${system}.primary_residues.dat | awk '{print $1+1}'`

### DOCK input file for minimizing on multi-grids
###########################################################################################
cat <<EOF >${system}.reference_multigridmin.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
ligand_atom_file                                             ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2 
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      yes
multigrid_score_secondary                                    no
multigrid_score_rep_rad_scale                                1
multigrid_score_vdw_scale                                    1
multigrid_score_es_scale                                     1
multigrid_score_number_of_grids                              ${num_of_grids}
EOF
############################################################################################

        setenv counter 0

foreach residue (`cat ../007.cartesian-min/${vendor}/${system}.primary_residues.dat`)
        
###########################################################################################
cat <<EOF >>${system}.reference_multigridmin.in
multigrid_score_grid_prefix${counter}                                 ${rootdir}/${system}/004.grid/${system}.resid_${residue}
EOF
###########################################################################################

	@ counter++ 
	end


###########################################################################################
cat <<EOF >>${system}.reference_multigridmin.in
multigrid_score_grid_prefix${counter}                        ${rootdir}/${system}/004.grid/${system}.resid_remaining
multigrid_score_fp_ref_mol                                   yes
multigrid_score_fp_ref_text                                  no
multigrid_score_footprint_ref                                ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2
multigrid_score_foot_compare_type                            Euclidean
multigrid_score_use_euc                                      yes
multigrid_score_normalize_foot                               no
multigrid_score_use_cor                                      no
multigrid_score_vdw_euc_scale                                1
multigrid_score_es_euc_scale                                 1
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
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
###########################################################################################

setenv primary_res ` cat ../007.cartesian-min/${vendor}/${system}.primary_residues.dat | sed -e 's/\n/ /g'`
echo ${primary_res}

python ${dockdir}/multigrid_fp_gen.py ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2 ${system}.resid ${system}.multigrid.in ${primary_res}

rm temp.mol2
rm ${system}.resid*.rec.grid.mol2
echo "Minimizing reference ligand in multigrid space..."
${dockdir}/dock6 -i ${system}.reference_multigridmin.in -o ${system}.reference_multigridmin.out
mv output_scored.mol2 ${system}.lig.multigridmin.mol2

