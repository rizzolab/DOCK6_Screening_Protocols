#!/usr/bin/perl -w

$INPDB   = "$ARGV[0]";

open (INPDB) or die "can't open input PDB file $INPDB\n";

while (<INPDB>) 
{
	@line=split(//);
        
	if ($line[77])
	{ 
		print $_ unless ($line[77] eq "H");
	}
	else 
	{
		print $_;
	}
}	

close (INPDB);
