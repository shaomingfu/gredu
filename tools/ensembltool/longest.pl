#!/usr/bin/perl

use strict;

if($#ARGV != 0) 
{
	die("usage: ./longest.pl list-file\n");
}

open(FILE, '<', $ARGV[0]) or die("open file $ARGV[0] error.");

my @xline;
while(<FILE>)
{
	chomp;
	my @line = split(',');

	die if($line[4] > $line[5]); 

	if($line[0] ne $xline[0])
	{
		print join(",", @xline) . "\n";
		@xline = @line;
	}
	elsif($line[5] - $line[4] > $xline[5] - $line[4])
	{
		@xline = @line;
	}
}
