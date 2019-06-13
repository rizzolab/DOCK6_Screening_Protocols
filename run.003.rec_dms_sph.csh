#!/bin/tcsh -fe

#
# This script will run dms and sphgen on the receptor. Prior to execution, confirm sphere cut-off
# distance (sphcut) and the maximum number of spheres you want to keep (maxkeep).
#


### Set some variables manually
set sphcut    = "8"
set maxkeep   = "75" 


### Set some paths
set dockdir   = "${DOCKHOMEWORK}/bin"
set amberdir  = "${AMBERHOMEWORK}/bin"
set moedir    = "${MOEHOMEWORK}/bin"
set rootdir   = "${VS_ROOTDIR}"
set dmsdir    = "${DMSHOMEWORK}/bin"
set masterdir = "${rootdir}/zzz.master"
set paramdir  = "${rootdir}/zzz.parameters"
set scriptdir = "${rootdir}/zzz.scripts"
set zincdir   = "${rootdir}/zzz.zinclibs"
set system    = "${VS_SYSTEM}"
set vendor    = "${VS_VENDOR}"


### Make the sphere directory
rm -fr ${rootdir}/${system}/003.spheres
mkdir -p ${rootdir}/${system}/003.spheres
cd ${rootdir}/${system}/003.spheres


### Copy some executables to the current location
cp -f ${scriptdir}/radii ./


### Make a small version of the pdb file that only includes residues with 20 A from the ligand
${scriptdir}/keep_close_atoms.pl ${rootdir}/${system}/001.lig-prep/${system}.lig.processed.mol2 ${rootdir}/${system}/002.rec-prep/pro.noH.pdb 20 > dms.pdb


### Run the program DMS 
time ${dmsdir}/dms dms.pdb -a -g ${system}.rec.dms.log -n -o ${system}.rec.dms.out


### Write an input file and generate sphere clusters on the molecular surface with sphgen
##################################################
cat <<EOF >INSPH
${system}.rec.dms.out
R
X
0.0
4.0
1.4
${system}.rec.sph
EOF
##################################################
${dockdir}/sphgen


### Convert the clusters (pruned shperes, not cluster 0) to a PDB file for viewing
##################################################
cat <<EOF >showsphere.in
${system}.rec.sph
-1
N
clustertemp
N
EOF
##################################################
${dockdir}/showsphere < showsphere.in


### Make a PDB file that containes all the pruned clusters for viewing
cat *clustertemp* >> temp.file
mv temp.file all.clust.pdb

if ( ! -s all.clust.pdb ) then
  echo "WARNING:   Spheres were not successfully generate!"
endif

### Make a PDB and a SPH file that contains only spheres close to the ligand
${scriptdir}/keep_close_spheres.pl ${rootdir}/${system}/001.lig-prep/${system}.lig.processed.mol2 all.clust.pdb ${sphcut} ${maxkeep}
mv temp.pdb ${system}.rec.clust.close.pdb
mv temp.sph ${system}.rec.clust.close.sph
mv OUTSPH ${system}.rec.sph.out


### Remove some files
rm -f INSPH OUTSPH showshpere.in sphgen showsphere *clustertemp* temp.pdb temp.sph radii dms


exit

