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


### Make a directory for compiling all of the docked results for a given vendor. If the top
### directory already exists, don't remove other vendor results.
if (! -e ${rootdir}/${system}/011.final-results/) then
	mkdir -p ${rootdir}/${system}/011.final-results/
endif

if (! -e ${rootdir}/${system}/011.final-results/${vendor}) then
	mkdir -p ${rootdir}/${system}/011.final-results/${vendor}
endif

if (! -e ${rootdir}/${system}/011.final-results/system-files/) then
        echo "Creating the system-files directory and copying the corresponding files\n"
	mkdir -p ${rootdir}/${system}/011.final-results/system-files/
	cd ${rootdir}/${system}/011.final-results/system-files/
	cp ${rootdir}/${system}/001.lig-prep/${system}.lig.am1bcc.mol2 ./
	cp ${rootdir}/${system}/007.cartesian-min/${vendor}/${system}.lig.am1bcc.min.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.mol2 ./
	cp ${rootdir}/${system}/002.rec-prep/${system}.rec.clean.pdb ./
	cp ${rootdir}/${system}/002.rec-prep/pro.noH.pdb ./${system}.rec.noH.pdb
	cp ${rootdir}/${system}/003.spheres/${system}.rec.clust.close.sph ./
	cp ${rootdir}/${system}/004.grid/box.pdb ./${system}.box.pdb
endif

foreach dockscore (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)

  echo ${dockscore}
#qsub here

rm -rf ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
mkdir -p ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank
cd ${rootdir}/${system}/011.final-results/${vendor}/${dockscore}_rank





### Write out mol2 files for each different sorting method and for each group - families and 
### clusterheads. Each resulting mol2 file will contain $max_size molecules.
if( -e ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 ) then
   echo " Writing the families and clusterhead mol2 files for each scoring metric\n"
   mkdir temp/
   cd temp/
   cp ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 ./
   python ${scriptdir}/break_into_mol.py ${system}.${vendor}.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2
   perl ${scriptdir}/concatenate_mol2_new_headers_descriptor.pl ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.final_sorted_dce_sum_families.csv
  cd ../
else
  echo " The file  ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.sorted_${dockscore}_${max_num}_dock.mol2 does not exist"
  exit
endif

set head_size = ${max_size}
@ head_size++
foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
	foreach group (clusterheads families)

		head -n ${head_size} ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.final_sorted_${score}_${group}.csv | awk -F "," '{print $1}' | sed '1d' > ${score}_${group}_zinc_codes.txt
                set listrank=1 #reset rank to 1 for a new clusterhead+scoring metric

                if(${group} == "clusterheads" )then
		   foreach mol2 (` cat ${score}_${group}_zinc_codes.txt `)
			rm -f temp/temp.mol2
			cp temp/${mol2}_new.mol2 temp/temp.mol2
			echo "##########                            List_Rank:      ${listrank}" | cat - temp/temp.mol2 > temp/temp2.mol2
			echo "##########                            From_List:      ${score}" | cat - temp/temp2.mol2 > temp/temp3.mol2
			cat temp/temp3.mol2 >> ${system}.${vendor}.${dockscore}_rank.final_sorted_${score}_${group}_${max_size}.mol2
                        @ listrank++ # increase rank
		   end
                else
                   foreach mol2 (` cat ${score}_${group}_zinc_codes.txt `)
                        rm -f temp/temp.mol2
                        cp temp/${mol2}_new.mol2 temp/temp.mol2
			echo "##########                            From_List:      ${score}" | cat - temp/temp.mol2 > temp/temp2.mol2
                        cat temp/temp2.mol2 >> ${system}.${vendor}.${dockscore}_rank.final_sorted_${score}_${group}_${max_size}.mol2
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
python ${scriptdir}/break_into_fp.py ${rootdir}/${system}/010.moe-postprocess/${vendor}/${dockscore}_rank/${system}.${vendor}.${dockscore}_rank.total_fp.txt
cd ../

foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
        foreach group (clusterheads families)
                foreach zincid (`cat ${score}_${group}_zinc_codes.txt `)
			cat temp/${zincid}.txt >> ${system}.${vendor}.${dockscore}_rank.footprints_${score}_${group}_${max_size}.txt
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

		python ${scriptdir}/break_into_fp.py ../${system}.${vendor}.${dockscore}_rank.footprints_${score}_${group}_${max_size}.txt
	end
end

#removed the ignore flag
#python -W ignore ${scriptdir}/plot_footprint_updated.py ../zinc_codes.txt ${cutoff} ${max_res}
python ${scriptdir}/plot_footprint_updated.py ../zinc_codes.txt ${cutoff} ${max_res}
cd ../


foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score)
	foreach group (clusterheads families)
		set filenames = ""

		foreach zincid (`cat ${score}_${group}_zinc_codes.txt `)
			set filenames = "${filenames} temp/${zincid}.pdf"
		end

		gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=${system}.${vendor}.${dockscore}_rank.footprint_plots_${score}_${group}_${max_size}.pdf ${filenames}
	end
end

mv ./temp/bad_val_id_list.txt .
rm -rf temp/
rm -f *zinc_codes.txt

end

exit



