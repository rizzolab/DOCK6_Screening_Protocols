#!/usr/bin/perl
use strict;
use POSIX;

unless (@ARGV){ die "Usage: perl $0 <library size> <chunk size> <multiplier>\n"; }

my $library_size = $ARGV[0];
my $chunk_size = $ARGV[1];
my $running_total = 0;
my $counter = 0;

while ($running_total < $library_size && $counter < 100){

	$running_total += $chunk_size;

	print "Chunk: ".$counter."\t\tChunk Size: ".$chunk_size."\t\tRunning Total: ".$running_total."\n";

	$counter++;
	my $temp = $chunk_size * $ARGV[2];
	$chunk_size = ceil($temp);

}

exit();
