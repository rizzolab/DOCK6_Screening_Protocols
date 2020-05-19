#!/usr/bin/perl -w

select (STDERR); $|++;
select (STDOUT); $|++;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

$FILE     = $ARGV[0];
$LEAPOUT  = $ARGV[1];
$outfile  = $ARGV[2];

#@existing_ter_res = ();

# Read in receptor PDB file by breaking on spaces. Deprecated
# Read in pdb file using unpack to allow for cases where there are
# no spaces between the coordinate columns
$i=0;
open (FILE) or die "can't open $FILE\n";
while (<FILE>) {
        if ( (/^ATOM/) ||(/^HETATM/) ) {
                chomp;

                @line = unpack("a6a5xa4xa3xa1a4xxxxa8a8a8",$_);
		$atonumb[$i]=trim($line[1]);
                $atoname[$i]=trim($line[2]);
                $resname[$i]=trim($line[3]);
                $resnumb[$i]=trim($line[5]);
                $X[$i]=trim($line[6]);
                $Y[$i]=trim($line[7]);
                $Z[$i]=trim($line[8]);
		$lines[$i]=$_;
		#print "\"$atoname[$i]\" \"$resname[$i]\"\n";
                $i++;
        }

#        push (@existing_ter_res,) if ( (/^TER/ ) {
#                chomp;

}
$natm=$i;
close (FILE);

# Search leap output for long bonds
# WARNING: There is a bond of 6.442116 angstroms between: 
# -------  .R<LEU 37>.A<C 18> and .R<PRO 38>.A<N 1>
# WARNING: There is a bond of 13.371721 angstroms between: 
# -------  .R<PRO 112>.A<C 13> and .R<GLU 113>.A<N 1>
# WARNING: The unperturbed charge of the unit: -2.000000 is not zero.

@longbondres = ();
open (LEAPOUT) or die "can't open $FILE\n";
while (<LEAPOUT>) {   
        if ( /WARNING: There is a bond of / ) {
		$_ = <LEAPOUT>;
		print $_; 
		@line=split(/[ >]/);
		$res1 = $line[3];  # select first residue of long bond
		$res2 = $line[9];  # select second residue of long bond
		if ($res1==$res2) {
			print "WARNING: Long bond inside residue $res1, TER not added\n";
		} else {
			push (@longbondres, $res2);	# draw long bond
			print "Long bond between residues detected: $res1 --- $res2\n";
		}
	}
}
push (@longbondres, -1);			# Mark the end of the long bond queue
close (LEAPOUT);

#fix histidines protonated incorrectly
for ($j=0; $j<$natm; $j++)
{
	next if (!($resname[$j] eq "CAL") and !($resname[$j] eq "MAG") and !($resname[$j] eq "ZIN") and !($resname[$j] eq "ZN") and 
 !($resname[$j] eq "MG") and !($resname[$j] eq "CA"));
	for ($i=0; $i<$natm; $i++)
	{
		next if !($resname[$i] eq "HIE");	# Process only HIE
	
		@line=split(//, $lines[$i]);		# Break into single characters
		#next if !($line[13] eq "H");		# For H atoms only char[13] is H
		next if !($line[77] eq "H");		# Check atom type column for hydrogens
		
		$distsq = ($X[$j]-$X[$i])*($X[$j]-$X[$i]) + ($Y[$j]-$Y[$i])*($Y[$j]-$Y[$i]) + ($Z[$j]-$Z[$i])*($Z[$j]-$Z[$i]);
		
		next if ($distsq > 2.25);		# Ignore if H is > 1.5A from ion

		# Change this whole residue HIE -> HID
		$curresno = $resnumb[$i];
		$k = $i - 1;			
		while ($resnumb[$k] == $curresno) 
		{ $k--; }
		
		$k++;					# k is now the first atom of the HIE residue
		while ($resnumb[$k] == $curresno)	# For the current HIE residue
		{
			$resname[$k]= "HID";		# Change the residue name to HID
			$lines[$k] =~ s/HIE/HID/;	# Replace HIE with HID
			$k++;				# Go to next atom in residue
		}
		#printf "%.2f %s", sqrt($dist)," ".$lines[$i];
	}
}


#store the cysteines with SG in array
for ($j=0;$j<$natm;$j++) 
{                                         
	if ( ($atoname[$j] eq "SG")&& ($resname[$j] eq "CYS" ) ) 
	{
		push(@sg, $j);
	}
}

foreach $cys1 (@sg)
{
	foreach $cys2 (@sg)
	{
		next if ($cys1 == $cys2);
 
		$Xdiff_squared = ($X[$cys1]-$X[$cys2])*($X[$cys1]-$X[$cys2]);
		$Ydiff_squared = ($Y[$cys1]-$Y[$cys2])*($Y[$cys1]-$Y[$cys2]);
		$Zdiff_squared = ($Z[$cys1]-$Z[$cys2])*($Z[$cys1]-$Z[$cys2]);
		$distance_squared= ($Xdiff_squared + $Ydiff_squared + $Zdiff_squared);

		#use SG bond distance squared as 5 ang. squared
		#actual SG bond length is 2.05A
                #we are allowing with sqrt(5)=2.24A  
		if ($distance_squared < 5.0) 
		{
			# Flag residue nums with SS bond
			$flagcode{$resnumb[$cys1]} = 1;
			$flagcode{$resnumb[$cys2]} = 1;

			#check if we already stored this atom
			$stored = 0;
			foreach $atomno (@ssbond)
                        {
                                if (($cys1 eq $atomno) or ($cys2 eq $atomno))
				{
					#set stored and exit the loop
					$stored = 1;
					last;
				}
                        }
			
			if ($stored == 0)
			{
				#save the atom nos. of the SS bond pair
				push(@ssbond, $cys1);
				push(@ssbond, $cys2);
			}	
		}														
	}
}

# write the ssbonds.txt file
open (SSBONDS, ">ssbonds.txt") or die "fix.long.bonds.perl:can't create output pdb file\n";
$len = scalar(@ssbond);
for ($i=0; $i<$len; $i+=2)
{
	#print "REC."."$resnumb[$ssbond[$i]]".".SG-"."REC."."$resnumb[$ssbond[$i+1]]".".SG";
	#print "," if ($i < $len-2);
	print SSBONDS "bond PRO."."$resnumb[$ssbond[$i]]".".SG PRO."."$resnumb[$ssbond[$i+1]]".".SG\n";
}
close(SSBONDS);

# Write a pdb with fixed HID, CYX and long-bond TERs
open (OUT, ">".$outfile) or die "fix.long.bonds.perl :can't create output pdb file\n";
for ($j=0;$j<$natm;$j++) {

	# Rename cysteines flagged as ss-bonded
	if ($flagcode{$resnumb[$j]}) {
		$resname[$j] = "CYX";
		$lines[$j] =~ s/CYS/CYX/;
	}
        
	# Printing file using printf caused atom type to get messed up, although leap reads it correctly
	@chars = split(//, $lines[$j]);  # Break into single characters
        next if ($chars[77] eq "H");     # Throw out all hydrogens

	if ($resnumb[$j] == $longbondres[0]) 
	{
		print OUT "TER\n";
		print "TER added before $resnumb[$j] to fix long bond\n";
		shift(@longbondres);
	} 

	printf OUT "ATOM %6d %-4s %3s %5d     %7.3f %7.3f %7.3f\n",
        $atonumb[$j],   
        $atoname[$j],
        $resname[$j],   
        $resnumb[$j],   
        $X[$j],         
        $Y[$j],         
        $Z[$j];       


	# Print TER after each protein chain (Leap guarantees the OXT after each chain)
	#print OUT "TER\n" if ( ($atoname[$j] eq "OXT")or($atoname[$j] eq "CAL")or($atoname[$j] eq "MAG")or($atoname[$j] eq "ZIN")or($atoname[$j] eq "CHL") );

	if ( ($atoname[$j] eq "OXT")or($atoname[$j] eq "CAL")or($atoname[$j] eq "MAG")or($atoname[$j] eq "ZIN")or($atoname[$j] eq "CHL") ) {
		print "TER added after $atonumb[$j] $atoname[$j] $resname[$j] $resnumb[$j]\n";
		print OUT "TER\n";
	}
         # Print TER card if the residue is NME
        print OUT "TER\n" if ( ($atoname[$j] eq "CH3") and $resname[$j] eq "NME");
}
close (OUT);

