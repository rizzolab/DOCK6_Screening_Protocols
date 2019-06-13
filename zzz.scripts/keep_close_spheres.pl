#!/usr/bin/perl -w

select (STDERR); $|++;
select (STDOUT); $|++;

$file1    = $ARGV[0];   # Ligand MOL2 file
$file2    = $ARGV[1];	# PDB Spheres from sphgen
$cutoff   = $ARGV[2];   # Cutoff distance from ligand 
$maxkeep  = $ARGV[3];	# MAx number of spheres to keep

$cutoffsqrd = $cutoff*$cutoff;

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

$i    = 0;
$yes  = 0;
$natm = 0;

# Read in MOL2 ligand
open (FILE, "$file1") or die "can't open $file1\n";
while (<FILE>) {
	if (/\@\<TRIPOS\>ATOM/) {
		$_=<FILE>;
		$yes   =1;
	}
	if (/\@\<TRIPOS\>BOND/) {
		$yes = 0;
	}
	if ($yes == 1) {
		@line=split; 
		unless ($line[0]) {
			shift(@line);
		}
		$sybtype[$i]=$line[5];
	}
	if (($yes == 1) && ($sybtype[$i] ne "H")) {
		@line=split; 
		unless ($line[0]) {
			shift(@line);
		}
		$atonumb[$i]=$line[0];
		$atoname[$i]=$line[1];
		$X[$i]=$line[2];
		$Y[$i]=$line[3];
		$Z[$i]=$line[4];
		$resnumb[$i]=$line[6];
		$resname[$i]=$line[7];                                    
		#print "FILE 1 $atoname[$i], $X[$i], $Y[$i], $Z[$i]\n";
	        $i++;
	}
}
close (FILE);

$ligand{$resname[0]} = ["$resname[0]"];

# Read in PDB spheres
# Code to read in PDB files by breaking on spaces is now deprecated and causes issues 
# Read in all.clust.pdb by using unpack
open (FILE, "$file2") or die "can't open $file2\n";
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

                $i++;
        }
}
close (FILE);

$natm = $i;
$z    = 0;
for ($j=0;$j<$natm;$j++) {                                         
	if ($ligand{$resname[$j]}) {
  		for ($k=0;$k<$natm;$k++) {
			$Xdiff_squared = ($X[$k]-$X[$j])*($X[$k]-$X[$j]);
			$Ydiff_squared = ($Y[$k]-$Y[$j])*($Y[$k]-$Y[$j]);
			$Zdiff_squared = ($Z[$k]-$Z[$j])*($Z[$k]-$Z[$j]);
			$distancesquared= ($Xdiff_squared + $Ydiff_squared + $Zdiff_squared);
			if ( ($distancesquared <= $cutoffsqrd) && (!($ligand{$resname[$k]})) ) {
                        	if ($closeres{$resnumb[$k]}) {
                                	next;
                                } else {
                                	$closeres{$resnumb[$k]} = $resnumb[$k];
					$map{$z} = ["$atonumb[$k]", "$atoname[$k]", "$resname[$k]", "$resnumb[$k]", "$X[$k]", "$Y[$k]", "$Z[$k]", "$distancesquared" ];
					$z++;
				}
			}
		}
	}
}

#Bubble Sort
for ($k=0; $k<$z; $k++) {
	for ($j=0; $j<$z-1-$k; $j++) {
		if ($map{$j+1}[7] < $map{$j}[7]) {
			$tmp       = $map{$j};
			$map{$j}   = $map{$j+1};
			$map{$j+1} = $tmp;
		}
	}
}

open (PDB, ">temp.pdb");
open (SPH, ">temp.sph");

$count   = 0;
$sphfile = "";
$sphline = "";

# Print out ranked
for ($l=0; $l<$z; $l++) {
	$sphnumber = $l+1;
	$sphtype   = "0.700";
	if ($count >= $maxkeep) {
		last;
	}
 	printf PDB "ATOM %6d %-4s %3s %5d     %7.3f %7.3f %7.3f %7.3f\n",
	$map{$l}[0],
	$map{$l}[1],
	$map{$l}[2],
	$map{$l}[3],
	$map{$l}[4],
	$map{$l}[5],
	$map{$l}[6],
	$map{$l}[7];

# NOTE THAT SPH FILE MUST CONTAIN LINES 80 CHARACTERS LONG
 	$sphline = sprintf "%5d%10.5f%10.5f%10.5f%8.3f%5d",
	$sphnumber,
	$map{$l}[4],
	$map{$l}[5],
	$map{$l}[6],
	$sphtype,
	$sphnumber;

	$sphfile = $sphfile.$sphline."                                \n";

	$count++;
}

$sphmax  = $count;

print SPH "cluster     1   number of spheres in cluster    $sphmax\n";
print SPH $sphfile;

close (PDB);
close (SPH);
exit;

#cluster     1   number of spheres in cluster    72
#    1  11.75500 -35.01400  25.70900   0.700    1
#    2   7.05300 -27.88000  24.84000   0.700    2
