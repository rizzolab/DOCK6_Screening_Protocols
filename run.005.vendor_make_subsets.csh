#!/bin/tcsh -fe

#
# This script will sift through the vendor catalogues to create subsets suitable for virtual
# screening. The environment variable DOCKHOME must be set prior to running this script. The
# min_rot and max_rot variables should be set to indicate the minimum and maximum number of 
# rotatable bonds that you will want to consider. Likewise, the charge limit should be set to 
# prune out molecules with absolute net charge greater than or less than the value specified.
#
# The chunk limit (chunklim) specifies the size of the largest chunk - the one that should have
# mostly the 0-rotatable bond molecules. Chunks get progressively smaller after the first one. 
# So if the chunklim is set too small for the size of the library used, it will never converge.
# The fix is to either increase the chunklim or use a smaller library. A good rule of thumb is
# chunklim >= 10% of total library.
#


### Set some variables manually
set min_rot    = "0"
set max_rot    = "15"
set chargelim  = "2.0"
set chunklim   = "60000"
set multiplier = "95"


### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"


### Make a directory for linking vendor library, don't remove previously-linked vendor libraries
### except from the same vendor
if (! -e ${rootdir}/${system}/005.zinclibs) then
	mkdir -p ${rootdir}/${system}/005.zinclibs/
endif
rm -f ${rootdir}/${system}/005.zinclibs/${vendor}
cd ${rootdir}/${system}/005.zinclibs/
ln -s ../../zzz.zinclibs/${vendor}/ ./${vendor}
cd ${zincdir}/${vendor}/


### Here, begin writing a tcsh script that can be run on the queue
##################################################
cat <<EOFA >${vendor}.prep.csh
#!/bin/tcsh -fe


### This is to preserve the variable 'chunklim' within the tcsh script
set chunklim = ${chunklim}


### cd into the vendor catalog (which should already be downloaded from zinc)
cd ${zincdir}/${vendor}/


### Concatenate the vendor downloaded subsets together
zcat ${zincdir}/${vendor}/*.mol2.gz > input.mol2
set num_mol = \`grep -c "MOL" input.mol2\`
echo "Total number of molecules = \${num_mol}"


### Precompute the approximate number of chunks that will be created. Kill the program if it 
### exceeds 100 chunks. This is changeable, but right now it is hard-coded so that the scaling
### -factor does not try to make an infinite number of chunks.
set temp = "0"
set temp_chunklim = "${chunklim}"
set counter = "0"

while (\${temp} < \${num_mol})
	@ temp = \${temp} + \${temp_chunklim}
	@ temp_chunklim = \${temp_chunklim} * ${multiplier} / 100
	@ counter++

	if (\${counter} > 100) then
		echo "The current settings will divide the library into more than 100 chunks. Increase the chunklim or"
		echo "decrease the library size. Exiting."
		exit
	endif
end

echo "Total number of chunks = \${counter}"
echo "Max chunk size = \${chunklim}"
echo "Min chunk size = <\${temp_chunklim}"


### Iteratively read in the input.mol2, and spit out molecules to individual files by number of
### rotatable bonds
set num_rot = "${max_rot} - ${min_rot} + 1"
set counter = "0"

while (\${counter} < \${num_rot})

##################################################
cat <<EOFB >chunk.\${counter}.dock.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             input.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          yes
dbfilter_max_heavy_atoms                                     9999
dbfilter_min_heavy_atoms                                     0
dbfilter_max_rot_bonds                                       \${counter}
dbfilter_min_rot_bonds                                       \${counter}
dbfilter_max_molwt                                           9999.0
dbfilter_min_molwt                                           0.0
dbfilter_max_formal_charge                                   ${chargelim}
dbfilter_min_formal_charge                                   -${chargelim}
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              no
atom_model                                                   all
vdw_defn_file                                                ${DOCKHOMEWORK}/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               ${DOCKHOMEWORK}/parameters/flex.defn
flex_drive_file                                              ${DOCKHOMEWORK}/parameters/flex_drive.tbl
ligand_outfile_prefix                                        output.\${counter}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOFB
##################################################

	${dockdir}/dock6 -i chunk.\${counter}.dock.in -o chunk.\${counter}.dock.out
	set tempnum = \`grep -c "MOL" output.\${counter}_scored.mol2\`
	echo "number of mols with \${counter} rotatable bond(s) = \${tempnum}"
	cat output.\${counter}_scored.mol2 >> temp.mol2

	### Make a csv file for num of rotatable bonds
	grep -B 2 "Elapsed time for docking" chunk.\${counter}.dock.out | grep ZINC | awk -v c=\${counter} '{print \$2 "," c}' >> num_rot_bonds.dat

	@ counter++
end

rm chunk.*.dock.in chunk.*.dock.out output.*_scored.mol2



### Count and report how many molecules were pruned
set new_num_mol = \`grep -c "MOL" temp.mol2\`
set difference = "\${num_mol} - \${new_num_mol}"
echo "Number of molecules pruned = \${difference}"
echo "Total number of molecules after sorting and pruning = \${new_num_mol}"


### Save some space:
rm input.mol2


### Divide the combined mol2 (temp.mol2) into chunks. Larger chunks first (fewer rotatable bonds)
### and smaller chunks last to optimize BlueGene performance.
set mol_output = "0"
set counter = "0"
set max_lig = "0"

while (\${mol_output} < \${new_num_mol})

	@ max_lig = \${chunklim} + \${mol_output}

##################################################
cat <<EOFB >dock_chunk\${counter}.in
conformer_search_type                                        rigid
use_internal_energy                                          no
ligand_atom_file                                             temp.mol2
limit_max_ligands                                            yes
max_ligands                                                  \${max_lig}
EOFB
##################################################

	if (\${counter} == 0) then
		echo "skip_molecule                                                no" >> dock_chunk\${counter}.in
	else
		echo "skip_molecule                                                yes" >> dock_chunk\${counter}.in
		echo "initial_skip                                                 \${mol_output}" >> dock_chunk\${counter}.in
	endif

##################################################
cat <<EOFB >>dock_chunk\${counter}.in
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
use_internal_energy                                          no
flexible_ligand                                              no
bump_filter                                                  no
score_molecules                                              no
ligand_outfile_prefix                                        chunk\${counter}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOFB
##################################################

	echo "Chunk \${counter}:"
	echo "Molecules output this round = \${chunklim}"
	@ mol_output = \${mol_output} + \${chunklim}
	echo "Molecules output so far = \${mol_output}"
	@ chunklim = \${chunklim} * ${multiplier} / 100

	${dockdir}/dock6 -i dock_chunk\${counter}.in -o dock_chunk\${counter}.out
	@ counter++
end

rm dock_chunk*.in



### Save some space:
rm temp.mol2


EOFA
##################################################
### Here is where that giant tcsh script ends


### Write the qsub script
##################################################
#Create Jobfile for LIRED or SeaWulf
if (`hostname -f` == "lired.cm.cluster" || `hostname -f` == "login.cm.cluster" ) then
cat <<EOFA >${vendor}.qsub.csh
#!/bin/tcsh
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=1
#PBS -q long
#PBS -N ${vendor}.prep
#PBS -V

cd ${zincdir}/${vendor}/
tcsh ${vendor}.prep.csh

EOFA
endif

##################################################


### Submit the script
qsub ${vendor}.qsub.csh


exit


