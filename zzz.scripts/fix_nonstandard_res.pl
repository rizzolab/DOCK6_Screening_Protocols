#!/usr/bin/perl -w

# zzz.fix.nonstandard.res.perl
# reads in a rec.formaber.pdb file and the leap output for this file
# mutates non-standard residues to ALA, and saves a pdb file

die 'Usage: zzz.fix.nonstandard.res.perl RECEPTOR LIGAND OUTPDB' unless (@ARGV==3);
$RECEPTOR = $ARGV[0];  # rec.foramber.pdb
$LIGAND   = $ARGV[1];  # ligand.mol2
$OUTPDB   = $ARGV[2];  # pdb output file

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

# Read in receptor PDB file by breaking on spaces. Deprecated
# Read in pdb file using unpack to allow for cases where there are
# no spaces between the coordinate columns
print "zzz.fix.nonstandard.res.perl: Receptor $RECEPTOR\n";
$i=0;
$num_CA=0;
open (RECEPTOR) or die "can't open $RECEPTOR\n";
while (<RECEPTOR>) {
        if ( (/^ATOM/) ||(/^HETATM/) ) {
                chomp;

                @line = unpack("a6a5xa4xa3xa1a4xxxxa8a8a8",$_);
		next unless (@line > 8); # ignore blank or incomplete records
                #$rec_atonumb[$i]=trim($line[1]);
                $rec_atoname[$i]=trim($line[2]);
                $rec_resname[$i]=trim($line[3]);
                $rec_resnumb[$i]=trim($line[5]);
                $rec_X[$i]=trim($line[6]);
                $rec_Y[$i]=trim($line[7]);
                $rec_Z[$i]=trim($line[8]);
                $lines[$i]=$_."\n";
		#$lines[$i]=~ s/\r/\n/g; #replace CR with LF
                #print "\"$atoname[$i]\" \"$resname[$i]\"\n";
		$num_CA++ if ($rec_atoname[$i] eq "CA");
                $i++;
        }
	if ( /^TER/ ) {
		# Pass through TER's unchanged by appending them to previous line read
		# This should make even multiple TER's print correctly along with the
		# last ATOM line

		$lines[$i-1]=$lines[$i-1].$_."\n" unless ($i==0);

		# When i=0, we have not read in any atoms yet, so we can ignore
		# any TER's that show up here. 
	}
}

$recnatm=$i;
close (RECEPTOR);
print "zzz.fix.nonstandard.res.perl: Receptor has $recnatm atoms and $num_CA alpha carbons\n";

# Read in ligand mol2 file
$i=0;
$yes=0;  #means we are in the @<TRIPOS>ATOM section
print "zzz.fix.nonstandard.res.perl: Reading Ligand $LIGAND \n";
open (LIGAND) or die "zzz.fix.nonstandard.res.perl: can't open $LIGAND\n";
while (<LIGAND>) {
        if (/\@\<TRIPOS\>ATOM/) {
                $_=<LIGAND>; # skip the @<TRIPOS>ATOM line
                $yes=1;
        }
        if (/\@\<TRIPOS\>BOND/) {
                $yes = 0;   # done reading 
		last;       # break out of the loop, done reading molecule
        }
        if ($yes == 1) {
                @line=split;
		next unless (@line > 7); #ignore blank lines or incomplete records
                $lig_sybtype[$i]=$line[5];
		next if ($lig_sybtype[$i] eq "H");  # skip ligand hydrogens
                #$lig_atonumb[$i]=$line[0];
                #$lig_atoname[$i]=$line[1];
                $lig_X[$i]=$line[2];
                $lig_Y[$i]=$line[3];
                $lig_Z[$i]=$line[4];
                #$lig_resnumb[$i]=$line[6];
                #$lig_resname[$i]=$line[7];
                #print "LIGAND 1 $atoname[$i], $X[$i], $Y[$i], $Z[$i]\n";
                $i++;
        }
}
close (LIGAND);
$lignatm = $i;
print "zzz.fix.nonstandard.res.perl: Ligand has $lignatm atoms\n";

# Identify all non-standard residues, mutate to ALA if far from ligand
%amino_acids = ('ALA',1,'ARG',2,'ASN',3,'ASP',4,'CYS',5,'GLU',6,
'GLY',7,'HIS',8,'ILE',9,'LEU',10,'LYS',11,'MET',12,'PHE',13,'PRO',14,
'SER',15,'THR',16,'TRP',17,'TYR',18,'VAL',19,'GLN',20,'HEM',21,'Y2P',22,'GLH',23,
'CAL',24,'MAG',25,'ZIN',26,'CHL',27,'POT',28,'HIE',29,'ACE',30,'NME',31,'NHE',32);
# Y2P = phosphorylated tyrosine

for ($j=0; $j<$recnatm; $j++)
{
	next if exists $amino_acids{$rec_resname[$j]};  #ignore standard amino acids

	# mutate MSE residues to methionine
	# ATOM    208 SE   MSE B  29       9.409  55.334  13.406 25.24 25.24          SE
	# ATOM     96  SD  MET A  16      32.985 -34.967  47.036  1.00 59.46           S
	if ($rec_resname[$j] eq "MSE") {
		$lines[$j]=~ s/SE   MSE/ SD  MET/g;
		$lines[$j]=~ s/MSE/MET/g;
		next; # do not delete any atoms
	}

	# Convert CME->MET. Leave CE and SD unchanged. SG->CG. Delete OH and CZ.
	if ($rec_resname[$j] eq "CME") {
                $lines[$j]=~ s/SG  CME/CG  MET/g;
		$lines[$j]=~ s/ CME / MET /g;
		$lines[$j]="" if (($rec_atoname[$j] eq "OH")or($rec_atoname[$j] eq "CZ"));
		next; # do not delete any more atoms
	}

	# Convert CSD->CYS by deleting two extra oxygens OD1 and OD2 
	if ($rec_resname[$j] eq "CSD") {
		$lines[$j]=~ s/ CSD / CYS /g;
		$lines[$j]="" if (($rec_atoname[$j] eq "OD1")or($rec_atoname[$j] eq "OD2"));
		next; # do not delete any more atoms
	}

        # Convert CCS->CYS by just renaming the residue
        if ($rec_resname[$j] eq "CCS") {
                $lines[$j]=~ s/ CCS / CYS /g;
                next; # do not delete any more atoms
        }

        # Rename  PTR->Y2P by just renaming the residue
        if ($rec_resname[$j] eq "PTR") {
                $lines[$j]=~ s/OH  PTR/OG  Y2P/g;
                $lines[$j]=~ s/ PTR / Y2P /g;
                next; # do not delete any more atoms
        }

        # Convert TYS->TYR by just renaming the residue
        #if ($rec_resname[$j] eq "TYS") {
        #        $lines[$j]=~ s/ OH  TYR / S   TYR /g;
        #        $lines[$j]=~ s/ TYS / TYR /g;
        #        next; # do not delete any more atoms
        #}

	# Measure distance from bad residue to ligand
	$mindist2ligsq = 10000;   # Set min distance is 100A between lig and non-std residue
	for ($i=0; $i<$lignatm; $i++)
	{
		$distsq = ($rec_X[$j]-$lig_X[$i])**2 + ($rec_Y[$j]-$lig_Y[$i])**2 + ($rec_Z[$j]-$lig_Z[$i])**2;
		$mindist2ligsq = $distsq if ($distsq < $mindist2ligsq);
	}	
	print "zzz.fix.nonstandard.res.perl: $rec_resname[$j] $rec_resnumb[$j] $rec_atoname[$j] ";
	printf "%.2f\n", sqrt($mindist2ligsq);

	#next if ($rec_resname[$j] eq "MSE");

	# keep only backbone+CA+CB in non-standard residues
	# lines can be deleted simply by setting them equal to ""
	$lines[$j]="" unless (($rec_atoname[$j] eq "CA")or($rec_atoname[$j] eq "N")or($rec_atoname[$j] eq "CB")or($rec_atoname[$j] eq "C")or($rec_atoname[$j] eq "O"));  
	$lines[$j]=~ s/$rec_resname[$j]/ALA/g;
}

# Write a pdb with fixed non-standard residues
open (OUT, ">".$OUTPDB) or die "fix.nonstd.residues.perl: can't create output pdb file\n";
for ($j=0;$j<$recnatm;$j++) {

	#printf OUT "ATOM %6d %-4s %3s %5d     %7.3f %7.3f %7.3f\n",
        #	$atonumb[$j], $atoname[$j], $resname[$j], $resnumb[$j], $X[$j], $Y[$j], $Z[$j];       
	#print OUT $lines[$j]."\n";
	print OUT $lines[$j];

	# Print TER after each protein chain (Leap guarantees the OXT after each chain)
	print OUT "TER\n" if ( ($rec_atoname[$j] eq "OXT")or($rec_atoname[$j] eq "CAL")or($rec_atoname[$j] eq "MAG")or($rec_atoname[$j] eq "ZIN")or($rec_atoname[$j] eq "CHL"));
        # Print TER card if the PDB has a NME residue
	print OUT "TER\n" if ( ($rec_atoname[$j] eq "HH3") and $rec_resname[$j] eq "NME");
}
close (OUT);

