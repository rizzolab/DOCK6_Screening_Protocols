#!/bin/tcsh -fe



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
set rankdir   = "${rootdir}/zzz.rank_scripts"

#1st set
foreach score (dce_sum fps_es fps_sum fps_vdw totalScore fms_score vo_score hms_score descriptor_score ) 
  echo ${score}
  tcsh ${rankdir}/run.008.${score}_rank.csh
  cd ${VS_ROOTDIR}
end

#Depending on which queue you use chose for the scripts
#you might have to submit in batches hence the two sets


#2nd set
#foreach score (vo_score hms_score descriptor_score)
#  echo ${score}
#  tcsh ${rankdir}/run.008.${score}_rank.csh
#  cd ${VS_ROOTDIR}
#end
