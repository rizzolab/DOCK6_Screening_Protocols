#!/bin/tcsh -fe

#
# This script prepares the receptor. Input required is ${system}.rec.foramber.pdb, which was
# prepared in MOE by removing everything except the protein and saving as a PDB. (no hydrogens or
# charges needed).
#
# Note that the environment variable $AMBERHOME should be set externally to this script.
#


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

### Check to see if the ligand file exists
if ( ! -e ${rootdir}/${system}/001.lig-prep/${system}.lig.am1bcc.mol2 ) then
	echo "You have to prepare the ligand first. Exiting."
	exit
endif


### Check to see that if a cofactor exists, a corresponding am1bcc.mol2 also exists
if ( -e ${masterdir}/${system}.cof.moe.mol2 ) then
	if ( ! -e ${rootdir}/${system}/001.lig-prep/${system}.cof.am1bcc.mol2 ) then
		echo "Cofactor expected for ${system} could not be found. Exiting."
		exit
	endif
endif


### Make the rec-prep directory
rm -fr ${rootdir}/${system}/002.rec-prep/
mkdir -p ${rootdir}/${system}/002.rec-prep/
cd ${rootdir}/${system}/002.rec-prep/


### Check to see if the protein file is present, and prepare it for amber
if ( ! -e ${masterdir}/${system}.rec.foramber.pdb ) then

	if ( -e ${masterdir}/${system}.rec.moe.pdb ) then
		# If there is not a hand-made rec.foramber.pdb file, make it here from rec.moe.pdb
		echo "Convert sed ${system}.rec.moe.pdb -> ${system}.rec.foramber.pdb" | tee ${system}.report.txt
		${scriptdir}/foramber.sed ${masterdir}/${system}.rec.moe.pdb > ./${system}.rec.foramber.pdb
	else 
		echo "Receptor file does not seem to exist. Exiting."
		exit
	endif
else
	# Hand-made rec.foramber.pdb file exists, so copy it over
	cp ${masterdir}/${system}.rec.foramber.pdb ./${system}.rec.foramber.pdb
endif


### Remove unusual newlines from receptor
perl -pi -e 's/\r\n/\n/g' ${system}.rec.foramber.pdb


### Track the number of alpha carbons and TERs in the protein, warn if it changes
set num_CA  = `grep -c " CA " ${system}.rec.foramber.pdb`
set num_TER = `grep -c "TER" ${system}.rec.foramber.pdb`
echo "${system}.rec.foramber.pdb has ${num_CA} alpha carbons and ${num_TER} TERS" | tee ${system}.report.txt


### Deal with non-standard residues
${scriptdir}/fix_nonstandard_res.pl ${system}.rec.foramber.pdb ../001.lig-prep/${system}.lig.am1bcc.mol2 out.pdb | tee ${system}.report.txt
mv ${system}.rec.foramber.pdb ${system}.rec.foramber.pdb.old
mv out.pdb ${system}.rec.foramber.pdb


### Copy some parameter files over that leap will need
cp -fr ${paramdir}/ions.frcmod ./
cp -fr ${paramdir}/ions.lib ./
cp -fr ${paramdir}/gaff_cz_mass_fix.frcmod ./
cp -fr ${paramdir}/heme.frcmod ./
cp -fr ${paramdir}/heme.prep ./
cp -fr ${paramdir}/y2p.frcmod ./
cp -fr ${paramdir}/y2p.off ./
cp -fr ${paramdir}/vdw_AMBER_parm99.defn ./vdw.defn
cp -fr ${paramdir}/chem.defn ./
touch box.pdb    # Needed for dock grid program


###
### Naming convention from here on is PRO = protein; COF = cofactor; REC = receptor = PRO + COF
###


### Read protein in with leap to renumber residues from 1
echo "Using AMBER version located at $AMBERHOME"
echo "------------ LEAP RUN_000 SUMMARY -------------" 
echo "Purpose: Split residues, renumber residues from 1"
echo "Removing hydrogens from rec.foramber.pdb --> pro.noH.pdb"
${scriptdir}/remove_hydrogens.pl ${system}.rec.foramber.pdb > pro.noH.pdb
echo -n "atoms in pro.noH.pdb = "
grep -c ATOM pro.noH.pdb
echo "Running leap.000 (Renumbering): pro.noH.pdb -> ${system}.pro.parm+pro.crd = pro.000.pdb"

##################################################
cat << EOF > pro.leap.in
set default PBradii mbondi2
source oldff/leaprc.ff14SB
loadoff ions.lib
loadamberparams ions.frcmod
loadamberparams heme.frcmod
loadamberparams frcmod.ions234lm_126_tip3p
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams frcmod.tip3p
loadamberprep heme.prep
loadoff y2p.off
loadamberparams y2p.frcmod
REC = loadpdb pro.noH.pdb
saveamberparm REC ${system}.pro.parm pro.crd
charge REC
quit
EOF
##################################################

${amberdir}/tleap -s -f pro.leap.in >& ${system}.pro.000.leap.log
${amberdir}/ambpdb -p ${system}.pro.parm -tit "pdb.000" -c pro.crd > pro.000.pdb
echo -n "atoms in pro.000.pdb (protonated) = "
grep -c ATOM pro.000.pdb


### Check to see if there were any errors with leap
set num_CA  = `grep -c " CA " pro.000.pdb`
set num_TER = `grep -c "TER" pro.000.pdb`
echo "pro.000.pdb has ${num_CA} alpha carbons and ${num_TER} TERs"
echo -n "Number of missing heavy atoms added : "
grep -c "Added missing heavy atom" ${system}.pro.000.leap.log  
grep -A1 "WARNING"              ${system}.pro.000.leap.log 
grep -A3 "Splitting"            ${system}.pro.000.leap.log  
if ( { grep "FATAL" ${system}.pro.000.leap.log  >> ${system}.report.txt } ) then
        echo "Leap run had fatal errors. Exiting." >> ${system}.report.txt
        exit
endif


### Remake the pro.noH.pdb with the fixed receptor structure
echo "Removing hydrogens pro.000.pdb --> pro.noH.pdb"
${scriptdir}/remove_hydrogens.pl pro.000.pdb > pro.noH.pdb
echo -n "atoms in pro.noH.pdb = "
grep -c ATOM pro.noH.pdb


### Run leap an extra time to get correct residue numbers for long bonds
echo "------------ LEAP RUN_001 SUMMARY -------------"
echo "Purpose: Get correct residue numbers for long bonds, ss bonds" 
${amberdir}/tleap -s -f pro.leap.in >& ${system}.pro.001.leap.log
${amberdir}/ambpdb -p ${system}.pro.parm -tit "pdb.001" -c pro.crd > pro.001.pdb
echo -n "atoms in pro.001.pdb = "
grep -c ATOM pro.001.pdb
set num_CA  = `grep -c " CA " pro.001.pdb`
set num_TER = `grep -c "TER" pro.001.pdb`
echo "pro.001.pdb has ${num_CA} alpha carbons and ${num_TER} TERs"


### Generate SS bonds, fix HIE/HID, remove hydrogens
echo "Running fix_long_bonds.pl: ssbonds.txt, hie/hid, del hydrogens" 
${scriptdir}/fix_long_bonds.pl pro.001.pdb ${system}.pro.001.leap.log pro.noH.pdb
echo -n "atoms in pro.noH.pdb = "
grep -c ATOM pro.noH.pdb
if ( -s ssbonds.txt ) then
	set numss = `wc -l ssbonds.txt | awk '{print $1}'`
	echo "${numss} disulphide bonds are required"
else
	echo "0 disulphide bonds are required"
endif
echo "--------------------------------------------------------"


### Preapre the ligand file with antechamber
echo "Creating ligand prep file with antechamber"
${amberdir}/antechamber -i ../001.lig-prep/${system}.lig.am1bcc.mol2 -fi mol2  -o ${system}.lig.ante.mol2 -fo mol2 -dr n
#echo "Creating ligand pdb file with antechamber"
#${amberdir}/antechamber -i ../001.lig-prep/${system}.lig.am1bcc.mol2 -fi mol2  -o ${system}.lig.ante.pdb -fo pdb -dr n
echo "Creating ligand prep file with parmchk2"
${amberdir}/parmchk2 -i ${system}.lig.ante.mol2 -f mol2 -o ${system}.lig.ante.frcmod


### Prepare the cofactor file with antechamber, if it exists
if ( -e ../001.lig-prep/${system}.cof.am1bcc.mol2 ) then
	echo "Creating cofactor prep file with antechamber" 
	${amberdir}/antechamber -i ../001.lig-prep/${system}.cof.am1bcc.mol2 -fi mol2  -o ${system}.cof.ante.prep -fo prepi
	echo "Creating cofactor pdb file with antechamber" 
	${amberdir}/antechamber -i ../001.lig-prep/${system}.cof.am1bcc.mol2 -fi mol2  -o ${system}.cof.ante.pdb -fo pdb
	echo "Creating cofactor frcmod file with parmchk2" 
	${amberdir}/parmchk2 -i ${system}.cof.ante.prep -f prepi -o ${system}.cof.ante.frcmod
endif


### Make tleap input for Complex
##################################################
cat << EOF > com.leap.in
set default PBradii mbondi2
source oldff/leaprc.ff14SB
source leaprc.gaff
loadoff ions.lib
loadamberparams ions.frcmod
loadamberparams frcmod.ions234lm_126_tip3p
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams frcmod.tip3p
loadamberparams heme.frcmod
loadamberprep heme.prep
loadoff y2p.off
loadamberparams y2p.frcmod
loadamberparams gaff_cz_mass_fix.frcmod
PRO = loadpdb pro.noH.pdb
EOF
cat ssbonds.txt >> com.leap.in
cat << EOF >> com.leap.in
loadamberparams ${system}.lig.ante.frcmod
LIG = loadmol2 ${system}.lig.ante.mol2
EOF
##################################################

if ( -e ${masterdir}/${system}.cof.moe.mol2 ) then
	echo "Generating complex = pro+lig+cof"
	echo "loadamberparams ${system}.cof.ante.frcmod" >> com.leap.in
	echo "loadamberprep ${system}.cof.ante.prep" >> com.leap.in
	echo "COF = loadpdb ${system}.cof.ante.pdb" >> com.leap.in
	echo "REC = combine { PRO COF }" >> com.leap.in
	echo "saveamberparm COF ${system}.cof.parm ${system}.cof.ori.crd"  >> com.leap.in
else
	echo "Generating complex = pro+lig (no cof)"
	echo "REC = combine { PRO }" >> com.leap.in
endif
echo "COM = combine { REC LIG }" >> com.leap.in
echo "saveamberparm LIG ${system}.lig.parm ${system}.lig.ori.crd"  >> com.leap.in
echo "saveamberparm PRO ${system}.pro.parm ${system}.pro.ori.crd"  >> com.leap.in
echo "saveamberparm REC ${system}.rec.parm ${system}.rec.ori.crd"  >> com.leap.in
echo "saveamberparm COM ${system}.com.parm ${system}.com.ori.crd"  >> com.leap.in
echo "quit" >> com.leap.in


### Use leap to generate complex
echo "------------ LEAP RUN_002 SUMMARY -------------"
echo "Purpose: Generate complex with ssbonds"
${amberdir}/tleap -s -f com.leap.in >& ${system}.com.leap.log
${amberdir}/ambpdb -p ${system}.lig.parm -tit "lig" -c ${system}.lig.ori.crd > ${system}.lig.ori.pdb
${amberdir}/ambpdb -p ${system}.pro.parm -tit "pro" -c ${system}.pro.ori.crd > ${system}.pro.ori.pdb
${amberdir}/ambpdb -p ${system}.rec.parm -tit "rec" -c ${system}.rec.ori.crd > ${system}.rec.ori.pdb
${amberdir}/ambpdb -p ${system}.com.parm -tit "com" -c ${system}.com.ori.crd > ${system}.com.ori.pdb
echo -n "atoms in ${system}.lig.ori.pdb = "
grep -c ATOM ${system}.lig.ori.pdb
echo -n "atoms in ${system}.pro.ori.pdb = "
grep -c ATOM ${system}.pro.ori.pdb
echo -n "atoms in ${system}.rec.ori.pdb = "
grep -c ATOM ${system}.rec.ori.pdb
echo -n "atoms in ${system}.com.ori.pdb = "
grep -c ATOM ${system}.com.ori.pdb


### Run sander to minimize hydrogen positions
echo "Creating ori.mol2 files before minimization"

${amberdir}/ambpdb -p ${system}.lig.parm -c ${system}.lig.ori.crd -mol2 -sybyl > ${system}.lig.ori.mol2 
${amberdir}/ambpdb -p ${system}.pro.parm -c ${system}.pro.ori.crd -mol2 -sybyl > ${system}.pro.ori.mol2 
${amberdir}/ambpdb -p ${system}.rec.parm -c ${system}.rec.ori.crd -mol2 -sybyl > ${system}.rec.ori.mol2 
${amberdir}/ambpdb -p ${system}.com.parm -c ${system}.com.ori.crd -mol2 -sybyl > ${system}.com.ori.mol2

#${amberdir}/antechamber -i ${system}.lig.ori.0.mol2 -fi mol2 -at sybyl -o ${system}.lig.ori.mol2 -fo mol2 -dr n
#${amberdir}/antechamber -i ${system}.pro.ori.0.mol2 -fi mol2 -at sybyl -o ${system}.pro.ori.mol2 -fo mol2 -dr n
#${amberdir}/antechamber -i ${system}.rec.ori.0.mol2 -fi mol2 -at sybyl -o ${system}.rec.ori.mol2 -fo mol2 -dr n
#${amberdir}/antechamber -i ${system}.com.ori.0.mol2 -fi mol2 -at sybyl -o ${system}.com.ori.mol2 -fo mol2 -dr n
 
if ( -e ${masterdir}/${system}.cof.moe.mol2 ) then
	${amberdir}/ambpdb -p ${system}.cof.parm -c ${system}.cof.ori.crd -mol2 -sybyl ${system}.cof.ori.mol2 

# top2mol2 has been deprecated since Amber15 but can be rebuilt in the antechamber makefile if uncommented
#${amberdir}/top2mol2 -p ${system}.lig.parm -c ${system}.lig.ori.crd -o ${system}.lig.ori.mol2 -at sybyl -bt sybyl 
#${amberdir}/top2mol2 -p ${system}.pro.parm -c ${system}.pro.ori.crd -o ${system}.pro.ori.mol2 -at sybyl -bt sybyl
#${amberdir}/top2mol2 -p ${system}.rec.parm -c ${system}.rec.ori.crd -o ${system}.rec.ori.mol2 -at sybyl -bt sybyl
#${amberdir}/top2mol2 -p ${system}.com.parm -c ${system}.com.ori.crd -o ${system}.com.ori.mol2 -at sybyl -bt sybyl
#if ( -e ${masterdir}/${system}.cof.moe.mol2 ) then
#	${amberdir}/top2mol2 -p ${system}.cof.parm -c ${system}.cof.ori.crd -o ${system}.cof.ori.mol2 -at sybyl -bt sybyl
endif

##################################################
cat <<EOF >sander.in
01mi minimization
 &cntrl
    imin = 1, maxcyc = 100,
    ntpr = 10, ntx=1,
    ntb = 0, cut = 10.0,
    ntr = 1, drms=0.1,
    restraintmask = "!@H=",
    restraint_wt  = 1000.0
&end
EOF
##################################################

echo "---------------------------------------------------------"
echo "Minimizing complex with sander"
${amberdir}/sander -O -i sander.in -o sander.out -p ${system}.com.parm -c ${system}.com.ori.crd -ref ${system}.com.ori.crd -r ${system}.com.min.rst
${amberdir}/ambpdb -p ${system}.com.parm -tit "${system}.com.min" -c ${system}.com.min.rst > ${system}.com.min.pdb
grep "SANDER BOMB" sander.out  
grep -A1 NSTEP sander.out | tail -2

if (! -s ${system}.com.min.rst) then
	echo "Complex minimizaton failed! Terminating."
	exit
endif


### Run sander on ligand alone to see if gaff screwed up anything
echo "---------------------------------------------------------"

##################################################
cat <<EOF1 >sander.lig.in
01mi minimization
 &cntrl
    imin = 1, maxcyc = 1000,
    ntpr = 10, ntx=1,
    ntb = 0, cut = 10.0,
    ntr = 0, drms=0.1,
&end
EOF1
##################################################

echo "Minimizing unrestrained gas-phase ligand alone with sander"
${amberdir}/sander -O -i sander.lig.in -o sander.lig.out -p ${system}.lig.parm -c ${system}.lig.ori.crd -r ${system}.lig.only.min.rst
${amberdir}/ambpdb -p ${system}.lig.parm -c ${system}.lig.only.min.rst -mol2 -sybyl > ${system}.lig.only.min.mol2
#${amberdir}/antechamber -i ${system}.lig.only.min.0.mol2 -fi mol2 -at sybyl -o ${system}.lig.only.min.mol2 -fo mol2 -dr n 
#${amberdir}/top2mol2 -p ${system}.lig.parm -c ${system}.lig.only.min.rst -o ${system}.lig.only.min.mol2 -at sybyl -bt sybyl
grep "SANDER BOMB" sander.lig.out
grep -A1 NSTEP sander.lig.out | tail -2
echo -n "Minimizing Ligand 1000 steps alone rmsd "
python ${scriptdir}/calc_rmsd_mol2.py ${system}.lig.ori.mol2 ${system}.lig.only.min.mol2


### Run sander on cofactor alone (if it exists) to see if gaff screwed up anything
if ( -e ${masterdir}/${system}.cof.moe.mol2 ) then
	echo "---------------------------------------------------------"
	echo "Minimizing unrestrained gas-phase cofactor alone with sander"
	cp sander.lig.in sander.cof.in
	${amberdir}/sander -O -i sander.cof.in -o sander.cof.out -p ${system}.cof.parm -c ${system}.cof.ori.crd -r ${system}.cof.only.min.rst
	${amberdir}/ambpdb -p ${system}.cof.parm -c ${system}.cof.only.min.rst -mol2 -sybyl > ${system}.cof.only.min.mol2
	#${amberdir}/top2mol2 -p ${system}.cof.parm -c ${system}.cof.only.min.rst -o ${system}.cof.only.min.mol2 -at sybyl -bt sybyl
	grep "SANDER BOMB" sander.cof.out
	grep -A1 NSTEP sander.cof.out | tail -2
	echo -n "Minimizing Cofactor 1000 steps alone rmsd "
	python ${scriptdir}/calc_rmsd_mol2.py ${system}.cof.ori.mol2 ${system}.cof.only.min.mol2
else
	echo "No cofactor present to minimize"
endif


### Extract some files from the minimized complex
echo "---------------------------------------------------------"
echo "Extracting receptor with cpptraj"
echo "trajin ${system}.com.min.rst" > rec.ptraj.in
echo "strip :LIG" >> rec.ptraj.in
echo "trajout ${system}.rec.min.rst restart"  >> rec.ptraj.in
${amberdir}/cpptraj ${system}.com.parm rec.ptraj.in >& rec.ptraj.out
grep STRIP rec.ptraj.out 
echo "Writing receptor mol2"

#Yuzhang modification made for multiple ligands aka waters etc.
${amberdir}/ambpdb -p ${system}.rec.parm -c ${system}.rec.min.rst -mol2 -sybyl > ${system}.rec.min.mol2
#${amberdir}/antechamber -i ${system}.rec.min.0.mol2 -fi mol2 -at sybyl -o ${system}.rec.min.mol2 -fo mol2 -dr n 

#${amberdir}/top2mol2 -p ${system}.rec.parm -c ${system}.rec.min.rst -o ${system}.rec.min.mol2 -at sybyl -bt sybyl
echo "Creating ligand mol2 file"
echo "trajin ${system}.com.min.rst" > lig.ptraj.in
echo "strip !(:LIG)" >> lig.ptraj.in
echo "trajout ${system}.lig.min.rst restart"  >> lig.ptraj.in
${amberdir}/cpptraj ${system}.com.parm lig.ptraj.in >& lig.ptraj.out
grep STRIP lig.ptraj.out
${amberdir}/ambpdb -p ${system}.lig.parm -c ${system}.lig.min.rst -mol2 > ${system}.lig.min.0.mol2 
${amberdir}/ambpdb -p ${system}.com.parm -c ${system}.com.min.rst -mol2 -sybyl > ${system}.com.min.mol2
${amberdir}/antechamber -i ${system}.lig.min.0.mol2 -fi mol2 -at sybyl -o ${system}.lig.min.mol2 -fo mol2 -dr n
#${amberdir}/antechamber -i ${system}.com.min.0.mol2 -fi mol2 -at sybyl -o ${system}.com.min.mol2 -fo mol2 -dr n 
#${amberdir}/top2mol2 -p ${system}.lig.parm -c ${system}.lig.min.rst -o ${system}.lig.min.mol2 -at sybyl -bt sybyl 
#${amberdir}/top2mol2 -p ${system}.com.parm -c ${system}.com.min.rst -o ${system}.com.min.mol2 -at sybyl -bt sybyl


### Compute some RMSDs to check for consistency
echo -n "Minimized Ligand rmsd "
python ${scriptdir}/calc_rmsd_mol2.py ${system}.lig.ori.mol2 ${system}.lig.min.mol2
echo -n "Minimized Receptor rmsd "
python ${scriptdir}/calc_rmsd_mol2.py ${system}.rec.ori.mol2 ${system}.rec.min.mol2
echo -n "Minimized Complex rmsd "
python ${scriptdir}/calc_rmsd_mol2.py ${system}.com.ori.mol2 ${system}.com.min.mol2
python ${scriptdir}/clean_mol2.py ${system}.rec.min.mol2 ${system}.rec.python.mol2 
python ${scriptdir}/clean_mol2.py ${system}.lig.min.mol2 ${system}.lig.python.mol2

### Create the pdb that will be used for MD simulations 
### Now this pdb matches the mol2 used for docking
${amberdir}/ambpdb -p ${system}.rec.parm -c ${system}.rec.min.rst > ${system}.rec.clean.pdb



### Run check grid
##################################################
cat <<EOF >grid.in
compute_grids                  no
output_molecule                yes
box_file                       box.pdb
receptor_file                  ${system}.rec.python.mol2
receptor_out_file              ${system}.rec.clean.mol2
EOF
##################################################

${dockdir}/grid -v -i grid.in -o grid.out
echo -n "CHECK GRID: " 
grep "Total charge on" grid.out  


### Remove some extra files
rm -f showbox vdw.defn chem.defn box.pdb
rm -f antechamber tleap teLeap parmchk ambpdb sander
rm -f ANTE* ATOMTYPE.INF NEWPDB.PDB PREP.INF
rm -f ions.frcmod ions.lib parm.e16.dat gaff*frcmod y2p.* heme.*
#rm -f ${system}.rec.min.mol2 ${system}.rec.nomin.mol2 ${system}.rec.foramber.pdb 
#rm -f ${system}.com.* ${system}.lig.* ${system}.rec.leap* ${system}.rec.gas*
rm -f mdinfo grid.in sander.* ssbonds.txt


exit
