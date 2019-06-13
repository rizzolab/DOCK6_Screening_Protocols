
#!/bin/tcsh -fe

#
# Edit the lines below to match your desired environment variables. Source this script prior to
# running the others.
#
# Remeber to define VS_ROOTDIR and VS_SYSTEM

### For running the scripts on Cluster
if (`hostname -f` == "login1.cm.cluster" || `hostname -f` == "login2.cm.cluster" || `hostname -f ` == "rizzo.cm.cluster") then

	### DOCK home directory
        setenv DOCKHOMEWORK /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3

	### AMBER home directory
        setenv AMBERHOMEWORK /gpfs/software/amber/16_gpu/amber16/

	### MOE home directory
	setenv MOEHOMEWORK /gpfs/projects/rizzo/zzz.programs/moe_2016.0801/

        ## DMS home directory
        setenv DMSHOMEWORK /gpfs/projects/rizzo/zzz.programs/dms

	### Root directory (the directory where all run.xxx.csh scripts are located)
	setenv VS_ROOTDIR  /gpfs/projects/rizzo/leprentis/DJ1/VS_2or3

        ### MPI directory (this is where mpirun is located, compatible with dock6.mpi)
        #setenv VS_MPIDIR /gpfs/software/intel/parallel-studio-xe/2019/compilers_and_libraries/linux/mpi/intel64
        setenv VS_MPIDIR /gpfs/software/intel/parallel-studio-xe/2018_3/compilers_and_libraries/linux/mpi/intel64

	### System name
	setenv VS_SYSTEM 2or3 #example - gp41.outerpocket

	### Vendor name
	setenv VS_VENDOR zinc15 #example - cdiv, chbr nt1105 sp100309
	
	### Fragment Library (only needed for DN and GA)
	setenv FRAGLIB /gpfs/projects/rizzo/zzz.programs/dock6.9_mpiv2018.0.3/parameters

	### Anchor Library (only needed for DN)
	setenv ANCLIB /gpfs/projects/rizzo/leprentis/zinc1_ancs_freq

        ### Max number of molecules to pass to Moe
        setenv MAX_NUM_MOL 100000 #standard is 100K 

        echo " Running the VS protocol on the SeaWulf Cluster"
        exit
endif

if (`hostname -f` == "lired.cm.cluster") then

        ### DOCK home directory
        setenv DOCKHOMEWORK /gpfs/projects/rizzo/zzz.programs/dock6

        ### AMBER home directory
        setenv AMBERHOMEWORK /gpfs/projects/rizzo/zzz.programs/amber16

        ### MOE home directory
        setenv MOEHOMEWORK /gpfs/projects/rizzo/zzz.programs/moe_2016.0801/

        ## DMS home directory
        setenv DMSHOMEWORK /gpfs/projects/rizzo/zzz.programs/dms

        ### Root directory (the directory where all run.xxx.csh scripts are located)
        setenv VS_ROOTDIR /gpfs/scratch/tmcgee/VS_zika_postfusion_stem/vs-protocol

        ### MPI directory (this is where mpirun is located, compatible with dock6.mpi)
        setenv VS_MPIDIR /gpfs/software/intel_2016_update1/impi/5.1.2.150/intel64/

        ### System name
        setenv VS_SYSTEM zika.postfusion.trimer.stem
        ### Vendor name
	setenv VS_VENDOR zinc15 #example - cdiv, chbr

        echo " Running the VS protocol on the LIRED Cluster"
          



else 
	echo "There was a problem assigning env variables needed for the VS protocol; you may have to do it manually."
        echo "This error is a result of the hostname not being recognized"
        echo " Type the commmand  'hostname -f ' and add that result to this script"
endif

