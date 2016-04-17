#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: join.pl index human.longest");
}

open(FILE1, '<', $ARGV[0]) or die("open $ARGV[0] error.");
my %fam;
while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	$fam{$line[0]} = $line[1];
}

open(FILE2, '<', $ARGV[1]) or die("open $ARGV[1] error.");
while(<FILE2>)
{
	chomp;
	my @line = split(',');
	if(defined( $fam{$line[6]} ))
	{
		my $t = $fam{$line[6]};
		if($line[3] eq -1)
		{
			$t = 0 - $t;
		}
		print("$line[0] $t $line[2] 1\n");
	}
}
