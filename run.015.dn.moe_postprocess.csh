#!/bin/tcsh -fe

###
### This script will postprocess the output files from virtual screening on cluster. Briefly, it
### concatenates everything together, generates a ranked csv file consisting of ZINC ids and other
### descriptors, calculates some additional descriptors in MOE, clusters by fingerprint, and writes
### some final output files ranked by different scoring functions. It is these files that should be
### visually inspected for purchasing.
###
### Prior to running this script, set the max number of molecules that should be clustered. If the
### number is set to around 100,000 then you can expect the calculation to take about a week. (?)
###


set max_num = "${MAX_NUM_MOL}"


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

set wcl   = 12:00:00
set nodes = 1
set ppn   = 24
set queue = "long-24core"
@ numprocs = (${nodes} * ${ppn})

### If doing best first clustering, this should be set to "_bestfit", else ""
set type = "_bestfit"
set sim  = "0.95"
### or
# set type = ""
# set sim  = "0.75"

### Make a directory for compiling all of the docked results for a given vendor. If the top
### directory already exists, don't remove other vendor results.
if (! -e ${rootdir}/${system}/015.moe-denovo/) then
	mkdir -p ${rootdir}/${system}/015.moe-denovo/
endif

if (! -e ${rootdir}/${system}/015.moe-denovo/${vendor}) then
	mkdir -p ${rootdir}/${system}/015.moe-denovo/${vendor}
endif

### There are only three licenses for moe currently - only run three at a time

foreach score (dce_sum fms_score fps_es)
#foreach score (fps_sum fps_vdw totalScore)
#foreach score (vo_score hms_score descriptor_score)

  echo ${score}

  rm -rf ${rootdir}/${system}/015.moe-denovo/${vendor}/${score}_rank
  mkdir -p ${rootdir}/${system}/015.moe-denovo/${vendor}/${score}_rank
  cd ${rootdir}/${system}/015.moe-denovo/${vendor}/${score}_rank

echo "Submitting to queue..."


### Create a large csh script to be submitted to the queue
##################################################
cat <<EOFA >${system}.denovo.${score}.postprocess.csh 
#!/bin/tcsh -fe


### Copy the relevant files over
cp ${rootdir}/${system}/014.denovo-desc-rescore/${vendor}/${score}_rank/denovo.${score}_rank.output_scored.mol2 ${system}.denovo.${score}_rank.total.mol2
cp ${rootdir}/${system}/014.denovo-desc-rescore/${vendor}/${score}_rank/denovo.${score}_rank.output_footprint_scored.txt ${system}.denovo.${score}_rank.total_fp.txt

### 2. Get scores and descriptors from docked mol2, sort by Descriptor_Score, and save top "max_num".
### (duplicates are also removed at this stage)
echo "Getting descriptors from DOCK for the combined set..."
python ${scriptdir}/dock_csv_updated_descriptor.dn.py ${system}.denovo.${score}_rank.total.mol2 ${system}.denovo.${score}_rank ${max_num} 
mv ${system}.denovo.${score}_rank_sorted_totScore_${max_num}.csv ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_dock.csv
rm -f ${system}.denovo.${score}_rank_sorted_totScore_all.csv 


### 3. Generate ranked list of ZINC ids with at most max_num molecules
echo "Generating ranked csv file..."
cat ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_dock.csv | awk -F "," '{print \$1}' > ${system}.denovo.${score}_rank.${max_num}.codes_withheading.txt
sed '1d' ${system}.denovo.${score}_rank.${max_num}.codes_withheading.txt > ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_codes.txt
rm -f ${system}.denovo.${score}_rank.${max_num}.codes_withheading.txt


### 4. Chop up the unique mol2 into individual files, then concatenate them in the order of the ranked list
echo "Generating ranked mol2 file..."
mkdir temp/
cd temp/
python ${scriptdir}/break_into_mol.py ../${system}.denovo.${score}_rank.total.mol2

foreach mol2 ( \` cat ../${system}.denovo.${score}_rank.sorted_${score}_${max_num}_codes.txt \` )
	cat \${mol2}.mol2 >> ../${system}.denovo.${score}_rank.sorted_${score}_${max_num}_dock.mol2
end

cd ../

##################################################
#rm -f ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_codes.txt
rm -rf temp/


### 5. Calculate descriptors in MOE and cluster by fingerprints
### Make sure the three extra MOE function files are available:
if ( ! -e ${MOEHOMEWORK}/svl/db_extractmoleculename.svl || ! -e ${MOEHOMEWORK}/svl/db_sm_Extract.svl || ! -e ${MOEHOMEWORK}/svl/ph4.svl/ph4clust_jj.svl ) then
        echo "Warning: You have to make sure that three files (1)db_extractmoleculename.svl, (2)db_sm_Extract.svl and (3)ph4clust_jj.svl are copied from zzz.scripts to ${MOEHOMEWORK}/lib/svl/"
        exit
endif
###############################################################

echo "Calculating descriptors and clustering in MOE ..."

###############################################################
${moedir}/moebatch -exec "\
mdb1 = db_Open ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'create'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'mol', 'molecule'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'Weight', 'float'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'b_rotN', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'lip_don', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'lip_acc', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'lip_druglike', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'lip_violation', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'SlogP', 'float'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'FCharge', 'int'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'logS', 'float'];\
db_CreateField ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'chiral', 'int'];\
db_ImportMOL2 ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_dock.mol2', '${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'mol'];\
db_ExtractMoleculeName [ mdb1, 'name'];\
db_sm_Extract [mdb1, 'SMILES'];\
codes1 = [ 'Weight', 'b_rotN', 'lip_don', 'lip_acc', 'lip_druglike', 'lip_violation', 'SlogP', 'FCharge', 'logS', 'chiral'];\
codes2 = [ 'name', 'SMILES', 'CLUSTER_NO', 'Weight', 'b_rotN', 'lip_don', 'lip_acc', 'lip_druglike', 'lip_violation', 'SlogP', 'FCharge', 'logS', 'chiral'];\
QuaSAR_DescriptorMDB ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'mol', codes1 ];\
ph4_FingerprintMDB ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'mol', 'FP:BIT_MACCS', 0];\
ph4_ClusterMDB${type} ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'FP:BIT_MACCS', 'tanimoto', [ overlap:0.75, sim:${sim}, cfield:'CLUSTER_NO', esel:0, mfield:'mol'] ];\
db_Close mdb1;\
mdb2 = db_Open ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', 'read'];\
entries = db_Entries '${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb';\
db_ExportASCII ['${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.mdb', '${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.csv', codes2, entries, [delimiter:',',quotes:0,titles:1]];\
db_Close mdb2\
"
###############################################################


### 6. Combine the csvs containing MOE descriptors and DOCK scores, and print out csv files of just
###    the clusterheads, and of everything grouped by cluster into families 
echo "Combining DOCK and MOE descriptors..."
python ${scriptdir}/sort_and_save_updated_descriptor.py ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_dock.csv ${system}.denovo.${score}_rank.sorted_${score}_${max_num}_moe.csv ${system}.denovo.${score}_rank.final

### 7. Prepare SMILE strings file (only_smiles_and_id.out) for ZINC15 for purchasing analogs
sed '1d' ${system}.denovo.${score}_rank.final_sorted_${score}_families.csv > tmpfile; 
mv tmpfile ${system}.denovo.${score}_rank.final_sorted_${score}_families_no_first_line.csv

awk -F, '{print \$35" "\$1}' ${system}.denovo.${score}_rank.final_sorted_${score}_families_no_first_line.csv > only_smiles_and_id.out

exit

EOFA
################################################################
### Here is where that giant tcsh script ends


### Write the qsub script
################################################################
################################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOF >${system}.denovo.${score}.postprocess.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${score}.${system}.denovo.postprocess
#PBS -V
#PBS -q ${queue}

echo "Job started"
date
cd ${rootdir}/${system}/015.moe-denovo/${vendor}/${score}_rank
tcsh ${system}.denovo.${score}.postprocess.csh 
echo "Job finished"
date

EOF
#################################################################

### Write the Cluster submit file to slurm
##################################################
###if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
else
cat <<EOF >${system}.denovo.${score}.postprocess.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${score}.${system}.denovo.postprocess
#SBATCH --output=${score}.${system}.denovo.postprocess
#SBATCH -p rn-long

echo "Job started"
date
cd ${rootdir}/${system}/015.moe-denovo/${vendor}/${score}_rank
tcsh ${system}.denovo.${score}.postprocess.csh
echo "Job finished"
date

EOF
endif
####################################################

### Submit the script
qsub ${system}.denovo.${score}.postprocess.qsub.csh

end
exit


