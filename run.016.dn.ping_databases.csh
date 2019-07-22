#!/bin/tcsh -fe

#
# This script takes the final clustered csv files and converts things into mol2 files. For each
# scoring method, it will write a mol2 of the top 'max_size' mol2 files from the clusterheads
# csv and from the families csv, to seperate files. All header information is taken from one of
# the csv files (they are equivalent), and it is written over the header information in the
# original mol2 file. In addition, it will create footprint plots for all molecules that appear
# in any of the output mol2 files. Finally, it will also copy a few extra files into the root 
# directory for visualization and analysis purposes.
#
# You will want to tune max_size so that you don't have multi-mol2 files which are exceedingly 
# large. The current version of Chimera (circa March 2013) has trouble opening more than about
# 13,000 - 14,000 mol2 objects in ViewDock at one time, so plan accordingly.
#


### Set some variables manually 
set max_size = "1000"
set max_num  = "100000"
set cutoff   = "0.2"
set max_res  = "50"


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

set wcl   = 120:00:00
set nodes = 1
set ppn   = 24
set queue = "extended"
@ numprocs = (${nodes} * ${ppn})

### Make a directory for compiling all of the docked results for a given vendor. If the top
### directory already exists, don't remove other vendor results.
if (! -e ${rootdir}/${system}/016.denovo-final-results/) then
	mkdir -p ${rootdir}/${system}/016.denovo-final-results/
endif

if (! -e ${rootdir}/${system}/016.denovo-final-results/${vendor}) then
	mkdir -p ${rootdir}/${system}/016.denovo-final-results/${vendor}
endif

if (! -e ${rootdir}/${system}/016.denovo-final-results/system-files/) then
        echo "Creating the system-files directory and copying the corresponding files\n"
	mkdir -p ${rootdir}/${system}/016.denovo-final-results/system-files/
	cd ${rootdir}/${system}/016.denovo-final-results/system-files/
	cp ${rootdir}/${system}/001.lig-prep/${system}.lig.am1bcc.mol2 ./
	cp ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.python.min.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.pdb ./
	cp ${rootdir}/${system}/002.rec-prep/pro.noH.pdb ./${system}.rec.noH.pdb
	cp ${rootdir}/${system}/003.spheres/${system}.rec.clust.close.sph ./
	cp ${rootdir}/${system}/004.grid/box.pdb ./${system}.box.pdb
endif

foreach dockscore (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)

  echo ${dockscore}
#qsub here

rm -rf ${rootdir}/${system}/016.denovo-final-results/${vendor}/${dockscore}_rank
mkdir -p ${rootdir}/${system}/016.denovo-final-results/${vendor}/${dockscore}_rank
cd ${rootdir}/${system}/016.denovo-final-results/${vendor}/${dockscore}_rank





### Write out mol2 files for each different sorting method and for each group - families and 
### clusterheads. Each resulting mol2 file will contain $max_size molecules.
if( -e ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 ) then
   echo " Writing the families and clusterhead mol2 files for each scoring metric\n"
   mkdir temp/
   cd temp/
   cp ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 ./
   python ${scriptdir}/break_into_mol.py ${system}.denovo.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2
   perl ${scriptdir}/concatenate_mol2_new_headers_descriptor.pl ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.final_sorted_dce_sum_families.csv
  cd ../
else
  echo " The file  ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 does not exist"
  exit
endif

set head_size = ${max_size}
@ head_size++
foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
	foreach group (clusterheads families)

		head -n ${head_size} ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.final_sorted_${score}_${group}.csv | awk -F "," '{print $1}' | sed '1d' > ${score}_${group}_zinc_codes.txt
                set listrank=1 #reset rank to 1 for a new clusterhead+scoring metric

                if(${group} == "clusterheads" )then
		   foreach mol2 (` cat ${score}_${group}_zinc_codes.txt `)
			rm -f temp/temp.mol2
			cp temp/${mol2}_new.mol2 temp/temp.mol2
			echo "##########                            List_Rank:      ${listrank}" | cat - temp/temp.mol2 > temp/temp2.mol2
			echo "##########                            From_List:      ${score}" | cat - temp/temp2.mol2 > temp/temp3.mol2
			cat temp/temp3.mol2 >> ${system}.denovo.${dockscore}_rank.final_sorted_${score}_${group}_${max_size}.mol2
                        @ listrank++ # increase rank
		   end
                else
                   foreach mol2 (` cat ${score}_${group}_zinc_codes.txt `)
                        rm -f temp/temp.mol2
                        cp temp/${mol2}_new.mol2 temp/temp.mol2
			echo "##########                            From_List:      ${score}" | cat - temp/temp.mol2 > temp/temp2.mol2
                        cat temp/temp2.mol2 >> ${system}.denovo.${dockscore}_rank.final_sorted_${score}_${group}_${max_size}.mol2
                   end 
                endif  
		cat ${score}_${group}_zinc_codes.txt >> used_zinc_codes.txt
	end
end
rm -rf temp/
cat used_zinc_codes.txt | sort | uniq > zinc_codes.txt
rm -f used_zinc_codes.txt


### Write a final footprint.txt file for each score-group combo
echo "Writing footprint txt files for each scoring+group combo\n"

mkdir temp/
cd temp/
python ${scriptdir}/break_into_fp.py ${rootdir}/${system}/015.moe-denovo/${vendor}/${dockscore}_rank/${system}.denovo.${dockscore}_rank.unique_fp.txt
cd ../

foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
        foreach group (clusterheads families)
                foreach zincid (`cat ${score}_${group}_zinc_codes.txt `)
			cat temp/${zincid}.txt >> ${system}.denovo.${dockscore}_rank.footprints_${score}_${group}_${max_size}.txt
                end
        end
end

rm -rf temp/


### Write footprint plots for all molecules that show up in any of the output mol2 files.
echo "Create the footprint pdfs for each score+group combo"
mkdir temp/
cd temp/

foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
	foreach group (clusterheads families)

		python ${scriptdir}/break_into_fp.py ../${system}.denovo.${dockscore}_rank.footprints_${score}_${group}_${max_size}.txt
	end
end

python -W ignore ${scriptdir}/plot_footprint_updated.py ../zinc_codes.txt ${cutoff} ${max_res}
cd ../


foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
	foreach group (clusterheads families)
		set filenames = ""

		foreach zincid (`cat ${score}_${group}_zinc_codes.txt `)
			set filenames = "${filenames} temp/${zincid}.pdf"
		end

		gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=${system}.denovo.${dockscore}_rank.footprint_plots_${score}_${group}_${max_size}.pdf ${filenames}
	end
end

rm -rf temp/
rm -f *zinc_codes.txt


### Write the qsub script for zinc15 and pubchem 
##################################################
### For Maui
##################################################
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" ) then
cat <<EOFB >${system}.denovo.pingzinc.qsub.csh
#!/bin/tcsh
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1
#PBS -N ${system}.dn.pingzinc
#PBS -V
#PBS -q ${queue}

echo "Job started"
date
cd ${rootdir}/${system}/016.denovo-final-results/${vendor}
python ${scriptdir}/updated_similarity_mols_in_zinc15.py ${rootdir}/${system}/015.moe-denovo/only_smiles_and_id.out 0
.75 10 output.pingzinc_results

echo "Job finished"
date

EOFB
##################################################
cat <<EOFA >${system}.denovo.pingpubchem.qsub.csh
#!/bin/tcsh
#PBS -l walltime=${wcl}
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -N ${system}.dn.pubchem
#PBS -V
#PBS -q ${queue}

echo "Job started"
date
cd ${rootdir}/${system}/016.denovo-final-results/${vendor}
python ${scriptdir}/pubchem_ping_denovo.py ${rootdir}/${system}/015.moe-denovo/only_smiles_and_id.out 0.75 10 output.pingpubchem_results

echo "Job finished"
date

EOFA
##################################################
### For slurm
##################################################
else
cat <<EOF >${system}.denovo.pingzinc.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${system}.dn.pingzinc
#SBATCH --output=${system}.dn.pingzinc
#SBATCH -p rn-long
echo "Job started"
date
cd ${rootdir}/${system}/016.denovo-final-results/${vendor}
python ${scriptdir}/updated_similarity_mols_in_zinc15.py ${rootdir}/${system}/015.moe-denovo/only_smiles_and_id.out 0
.75 10 output.pingzinc_results

echo "Job finished"
date

EOF
##################################################
cat <<EOFB >${system}.denovo.pingpubchem.qsub.csh
#!/bin/tcsh
#SBATCH --time=${wcl}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ppn}
#SBATCH --job-name=${system}.dn.pingpubchem
#SBATCH --output=${system}.dn.pingpubchem
#SBATCH -p rn-long

echo "Job started"
date
cd ${rootdir}/${system}/016.denovo-final-results/${vendor}
python ${scriptdir}/pubchem_ping_denovo.py ${rootdir}/${system}/015.moe-denovo/only_smiles_and_id.out 0.75 10 output.pingpubchem_results

echo "Job finished"
date

EOFB
endif
##################################################



### Submit the scripts
qsub ${system}.denovo.pingzinc.qsub.csh
qsub ${system}.denovo.pingpubchem.qsub.csh

end


exit



