#!/usr/bin/perl
use strict;

#
# This program will read in a large multi-mol2 that has docked ligands accompanied by energy
# scores. It will go through the mol2 file once to read ZINC names and the DCE score, and in the
# case of duplicate ZINC names, it will only remember the best score. Then, it will go through the
# mol2 file a second time, writing out the mol2 object to a new file if it is the one with the
# best score. Or, if it is a unique ZINC name, it will just write it out. Then, it will open up
# the footprint text file that accompanies the mol2 file, and repeat the second half of the 
# process.
#


unless (@ARGV){die "Usage: perl $0 <input mol2> <input fp txt> <output mol2> <output fp txt>\n"; }


open (INPUT, $ARGV[0]);
my %name_score;
my $name_flag = 0;
my $score_flag = 0;
my ($temp1, $temp2, $this_name, $this_score);


# Go through the first time and create a hash of unique ZINC ids and keep the best corresponding
# descriptor score
while (<INPUT>) {

	if ($_ =~ m/Name/){
		chomp ($_);
		($temp1, $temp2, $this_name) = split(' ', $_); 
		$name_flag = 1;
	}

	if ($_ =~ m/Descriptor_Score/){
		chomp ($_);
		($temp1, $temp2, $this_score) = split(' ', $_);
		$score_flag = 1;
	}

	if ($name_flag == 1 && $score_flag == 1){

		if (exists($name_score{$this_name})){

			if ($name_score{$this_name} > $this_score){
				$name_score{$this_name} = $this_score;
			}

		} else {
			$name_score{$this_name} = $this_score;
		}

		$name_flag = 0;
		$score_flag = 0;
	}
}

close (INPUT);

my @temp = keys(%name_score);

#print "there are ".scalar(@temp)." unique keys\n";
#foreach (@temp){ print $_."\t\t".$name_score{$_}."\n"; }


open (INPUT, $ARGV[0]);
open (OUTPUT, ">>", $ARGV[2]);
my $this_mol2;
my $keep_flag = 0;
my $discard_flag = 0;
$name_flag = 0;
$score_flag = 0;
#my ($temp1, $temp2, $this_name, $this_score);


# Go through a second time to write out to file
while (<INPUT>){

	if ($_ =~ m/Name/){
		chomp ($_);
		($temp1, $temp2, $this_name) = split(' ', $_); 
		$name_flag = 1;
		$_ = $_."\n";

		if ($keep_flag == 1){
			print OUTPUT $this_mol2;
			$this_mol2 = "";
			$keep_flag = 0;
		}
		if ($discard_flag == 1){
			$this_mol2 = "";
			$discard_flag = 0;
		}
	}

	$this_mol2 = $this_mol2.$_;

	if ($_ =~ m/Descriptor_Score/){
		chomp ($_);
		($temp1, $temp2, $this_score) = split(' ', $_);
		$score_flag = 1;
	}

	if ($name_flag == 1 && $score_flag == 1){
		if ($name_score{$this_name} == $this_score){
			$keep_flag = 1;
		} else {
			$discard_flag = 1;
		}

		$name_flag = 0;
		$score_flag = 0;
	}
}

# Print the last one
if ($keep_flag == 1){
	print OUTPUT $this_mol2;
}


close (OUTPUT);
close (INPUT);



open (INPUT, $ARGV[1]);
open (OUTPUT, ">>", $ARGV[3]);
my $this_fptxt;
my $keep_flag = 0;
my $discard_flag = 0;
$name_flag = 0;
$score_flag = 0;
#my ($temp1, $temp2, $this_name, $this_score);


# Go through the footprint file and write out to file
while (<INPUT>){

	if ($_ =~ m/Molecule/){
		chomp ($_);
		($temp1, $temp2, $this_name) = split(' ', $_); 
		$name_flag = 1;
		$_ = $_."\n";

		if ($keep_flag == 1){
			print OUTPUT $this_fptxt;
			$this_fptxt = "";
			$keep_flag = 0;
		}
		if ($discard_flag == 1){
			$this_fptxt = "";
			$discard_flag = 0;
		}
	}

	$this_fptxt = $this_fptxt.$_;

	if ($_ =~ m/Descriptor_Score/){
		chomp ($_);
		($temp1, $temp2, $this_score) = split(' ', $_);
		$score_flag = 1;
	}

	if ($name_flag == 1 && $score_flag == 1){
		if ($name_score{$this_name} == $this_score){
			$keep_flag = 1;
		} else {
			$discard_flag = 1;
		}

		$name_flag = 0;
		$score_flag = 0;
	}
}

# Print the last one
if ($keep_flag == 1){
	print OUTPUT $this_fptxt;
}


close (OUTPUT);
close (INPUT);




exit();





