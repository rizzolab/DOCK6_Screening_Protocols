#!/usr/bin/perl

select (STDERR); $|++;
select (STDOUT); $|++;

$file           = $ARGV[0];
$molname        = $ARGV[1];
$replacemolname = $ARGV[2];

unless   ( ($ARGV[0]) && ($ARGV[1]) && ($ARGV[2]) ) {
        print STDERR "ERROR Please Specify ARGUMENTS\n";
        print STDOUT "ERROR Please Specify ARGUMENTS\n";
        die;
}

$found = 0;
$i     = 0;
$sum   = 0;
$tchrg = 0;
$newcode = "ZZZ";
if (-B $file) {
        open (FILE, "gzip -dc $file | ") or die "can't open $file\n";
} else {
        open (FILE, "$file") or die "can't open $file\n";
}
	BIG:while (<FILE>) {
		@line=split;
		unless ($line[0]) {
			shift(@line);
		}
		$number         = $line[0];
		$oriname        = $line[1];
		$x              = $line[2];
		$y              = $line[3];
		$z              = $line[4];
		$type           = $line[5];
		$resnumb        = $line[6];
		$resname        = $line[7];

		if ($replacemolname eq "yes" ) {
                	$resnumb        = "1";
                	$resname        = "LIG";
		}

		$oricharge      = $line[8];
		@typelist=split(/\./,$type);  #the element before the .
		$atomcode=$typelist[0];

			if (/\@\<TRIPOS\>MOLECULE/) {
				print "@<TRIPOS>MOLECULE\n";
				$_=<FILE>;
				$oldcode = $_;
				chomp($oldcode);
				if ($replacemolname eq "yes" ) {
					$newcode = $molname.".lig";
					print "$newcode\n";
				} else {
					$newcode = $oldcode;
					print "$newcode\n";
				}
				$_=<FILE>;
				@line=split;
				unless ($line[0]) {
					shift(@line);
				}
				$natom  = $line[0];
				$nbonds = $line[1];
				$nset   = $line[2];

				if ($replacemolname eq "yes" ) {
					$nset   = "1";
				}

				$ntemp1 = 0;
				$ntemp2 = 0;
				printf "%4d %4d %4d %4d %4d\n",
				$natom,$nbonds,$nset,$ntemp1,$ntemp2;
				$_=<FILE>;
			}
			if (/\@\<TRIPOS\>ATOM/) {
				print "@<TRIPOS>ATOM\n";
				$i=0;
				$found = 1;
				%counter = ();
			}
			if (/\@\<TRIPOS\>BOND/) {
				print "@<TRIPOS>BOND\n";
				$found = 0;
				next;
			}
                	if (/\@\<TRIPOS\>SUBSTRUCTURE/) {
                        	print $_;
                        	$_=<FILE>;
				print "1     LIG    1";
                        	print " \n";

	                        $tchrg = $sum;

	                        $n = int($tchrg);
	                        if ($tchrg > 0) {
	                          if ($tchrg - $n > 0.5) {
	                           $tchrg = $n + 1;
	                          }
	                          else {
	                           $tchrg = $n
	                          }
	                        }
	                        else {
	                          if ( $n - $tchrg > 0.5) {
	                           $tchrg = $n - 1;
	                          }
	                          else {
	                           $tchrg = $n
	                          }
	                        }  

	                        print STDERR "NEWCODE = $newcode : OLDCODE = $oldcode : CHARGE = $tchrg\n";

				$sum   = 0;
				$tchrg = 0;

				next BIG;
                	}

			if ($found == 1) {
				$counter{$atomcode}++;
				$newname = $atomcode.$counter{$atomcode};
			}
			if ($found == 0) {
				print $_;
			}
			if (($found == 1) && ($i ge 1)) {
				printf "%-5d %-4s %10.4f%10.4f%10.4f %-7s %-3s%-3s %8.4f \n",
				$number,$newname,$x,$y,$z,$type,$resnumb,$resname,$oricharge;
				$sum = $sum + $oricharge;
			}
	$i++;
	}
