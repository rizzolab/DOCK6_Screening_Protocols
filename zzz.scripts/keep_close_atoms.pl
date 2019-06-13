#!/usr/bin/perl -w

select (STDERR); $|++;
select (STDOUT); $|++;

$ligand    = $ARGV[0];
$receptor  = $ARGV[1];
$cutoff    = $ARGV[2];

$cutoffsqrd = $cutoff*$cutoff;

$i    = 0;
$yes  = 0;
$natm = 0;

# Read ligand mol2 file
open (FILE, "$ligand") or die "can't open $ligand\n";
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

# Read receptor without hydrogens and core out the binding site for dms
# this is needed because dms fails if the receptor is too big
open (FILE, "$receptor") or die "can't open $receptor\n";
while (<FILE>) {             
        if ( (/^ATOM/) ||(/^HETATM/) ) {
		chomp;

		@line = unpack("a6a5xa4xa3xa1a4xxxxa8a8a8",$_);
                $recname[$i]=$line[0];
                $atonumb[$i]=$line[1];
                $atoname[$i]=$line[2];
                $resname[$i]=$line[3];
                $resnumb[$i]=$line[5];
                $X[$i]=$line[6];
                $Y[$i]=$line[7];
                $Z[$i]=$line[8];


                $i++;
        }
}
close (FILE);

$natm = $i;
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
						# print "CLOSE RESIDUE = $closeres{$resnumb[$k]}\n";
                                	}
				}
		}
	}
}

for ($j=0;$j<$natm;$j++) {                                         
	if ($closeres{$resnumb[$j]}) {
		print $recname[$j];
	        print $atonumb[$j];
		print " ";
	        print $atoname[$j];
		print " ";
	        print $resname[$j];
		print "  ";
	        print $resnumb[$j];
		print "    ";
	        print $X[$j];
	        print $Y[$j];
	        print $Z[$j];
	        print "\n";
	}
}



