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
		print OUTPUT "##########                            Name_DOCK:      ".$fields[0]."\n";
		print OUTPUT "##########                             Name_MOE:      ".$fields[1]."\n";
		print OUTPUT "##########                              Cluster:      ".$fields[2]."\n";
		print OUTPUT "##########                         Cluster_size:      ".$fields[3]."\n";
		print OUTPUT "##########                 TotalScore_(FPS+DCE):      ".$fields[4]."\n";
		print OUTPUT "##########                     Continuous_Score:      ".$fields[5]."\n";
		print OUTPUT "##########                Continuous_vdw_energy:      ".$fields[6]."\n";
		print OUTPUT "##########                 Continuous_es_energy:      ".$fields[7]."\n";
		print OUTPUT "##########            Internal_energy_repulsive:      ".$fields[8]."\n";
		print OUTPUT "##########           Footprint_Similarity_Score:      ".$fields[9]."\n";
		print OUTPUT "##########                          FPS_vdw_fps:      ".$fields[10]."\n";
		print OUTPUT "##########                           FPS_es_fps:      ".$fields[11]."\n";
		print OUTPUT "##########                           FPS_hb_fps:      ".$fields[12]."\n";
		print OUTPUT "##########                    FPS_vdw_fp_numres:      ".$fields[13]."\n";
		print OUTPUT "##########                     FPS_es_fp_numres:      ".$fields[14]."\n";
		print OUTPUT "##########                     FPS_hb_fp_numres:      ".$fields[15]."\n";
		print OUTPUT "##########                          Num_H-bonds:      ".$fields[16]."\n";
		print OUTPUT "##########                       DOCK_rot_bonds:      ".$fields[17]."\n";
                print OUTPUT "##########                  Pharmacophore_Score:      ".$fields[18]."\n";
                print OUTPUT "##########                Property_Volume_Score:      ".$fields[19]."\n";
                print OUTPUT "##########                       Tanimoto_Score:      ".$fields[20]."\n";
                print OUTPUT "##########  Hungarian_Matching_Similarity_Score:      ".$fields[21]."\n";
                print OUTPUT "##########                     Descriptor_Score:      ".$fields[22]."\n";
		print OUTPUT "##########                        MOE_rot_bonds:      ".$fields[24]."\n";
		print OUTPUT "##########                     Molecular_weight:      ".$fields[23]."\n";
		print OUTPUT "##########                   Num_chiral_centers:      ".$fields[25]."\n";
		print OUTPUT "##########                      Lipinski_donors:      ".$fields[26]."\n";
		print OUTPUT "##########                   Lipinski_acceptors:      ".$fields[27]."\n";
		print OUTPUT "##########                    Lipinski_druglike:      ".$fields[28]."\n";
		print OUTPUT "##########                  Lipinski_violations:      ".$fields[29]."\n";
		print OUTPUT "##########                                SlogP:      ".$fields[30]."\n";
		print OUTPUT "##########                        Formal_charge:      ".$fields[31]."\n";
		print OUTPUT "##########                                 logS:      ".$fields[32]."\n";
		print OUTPUT "##########                    Ligand_efficiency:      ".$fields[33]."\n";
		print OUTPUT "##########                               SMILES:      ".$fields[34];

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

