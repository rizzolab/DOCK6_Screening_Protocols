#!/usr/bin/perl
use strict;

#
# This program reads in a csv file of ZINC names and descriptors. It will go through each line,
# writing the contents of the line to the header of a mol2 file, then it will find the corresponding
# molecule and write the molecule information to the same mol2 file. Header information in the
# original mol2 file is ignored. The program assumes the mol2 files were already made and are in
# the current directory.
#


unless (@ARGV){die "Usage: perl $0 <input csv file>\n";}


open (INPUT, $ARGV[0]);
while (<INPUT>){

	my @fields = split(',',$_);

	unless ($fields[0] =~ m/name_dock/){

		# Open a new output file
		open (OUTPUT, ">", $fields[0]."_new.mol2");
		
		# Print all the header information into it
		print OUTPUT "##########             Name_DOCK:      ".$fields[0]."\n";
		print OUTPUT "##########              Name_MOE:      ".$fields[1]."\n";
		print OUTPUT "##########               Cluster:      ".$fields[2]."\n";
		print OUTPUT "##########          Cluster_size:      ".$fields[3]."\n";
		print OUTPUT "##########  TotalScore_(FPS+DCE):      ".$fields[4]."\n";
		print OUTPUT "##########      DCE_sum_(vdw+es):      ".$fields[5]."\n";
		print OUTPUT "##########               DCE_vdw:      ".$fields[6]."\n";
		print OUTPUT "##########                DCE_es:      ".$fields[7]."\n";
		print OUTPUT "##########       Internal_energy:      ".$fields[8]."\n";
		print OUTPUT "##########      FPS_sum_(vdw+es):      ".$fields[9]."\n";
		print OUTPUT "##########               FPS_vdw:      ".$fields[10]."\n";
		print OUTPUT "##########                FPS_es:      ".$fields[11]."\n";
		print OUTPUT "##########                FPS_hb:      ".$fields[12]."\n";
		print OUTPUT "##########        FPS_vdw_numres:      ".$fields[13]."\n";
		print OUTPUT "##########         FPS_es_numres:      ".$fields[14]."\n";
		print OUTPUT "##########         FPS_hb_numres:      ".$fields[15]."\n";
		print OUTPUT "##########           Num_H-bonds:      ".$fields[16]."\n";
		print OUTPUT "##########        DOCK_rot_bonds:      ".$fields[17]."\n";
		print OUTPUT "##########         MOE_rot_bonds:      ".$fields[19]."\n";
		print OUTPUT "##########      Molecular_weight:      ".$fields[18]."\n";
		print OUTPUT "##########    Num_chiral_centers:      ".$fields[20]."\n";
		print OUTPUT "##########       Lipinski_donors:      ".$fields[21]."\n";
		print OUTPUT "##########    Lipinski_acceptors:      ".$fields[22]."\n";
		print OUTPUT "##########     Lipinski_druglike:      ".$fields[23]."\n";
		print OUTPUT "##########   Lipinski_violations:      ".$fields[24]."\n";
		print OUTPUT "##########                 SlogP:      ".$fields[25]."\n";
		print OUTPUT "##########         Formal_charge:      ".$fields[26]."\n";
		print OUTPUT "##########                  logS:      ".$fields[27]."\n";
		print OUTPUT "##########     Ligand_efficiency:      ".$fields[28]."\n";
		print OUTPUT "##########                SMILES:      ".$fields[29];

		# Print the mol2 information into it
		open (INPUT2, $fields[0].".mol2");
		my @input_mol2 = <INPUT2>;
		my @input_mol2_b = grep !/^#/, @input_mol2;
		close (INPUT2);
		print OUTPUT @input_mol2_b;

		# Close the output file
		close (OUTPUT);

	}

}
close (INPUT);


exit();

