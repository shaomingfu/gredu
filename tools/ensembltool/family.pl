#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: family.pl human.longest mouse.longest");
}

open(FILE1, '<', $ARGV[0]) or die("open $ARGV[0] error.");
my %fam1;
while(<FILE1>)
{
	chomp;
	my @line = split(',');
	if(defined( $fam1{$line[6]} ))
	{
		$fam1{$line[6]} = $fam1{$line[6]} + 1;
	}
	else
	{
		$fam1{$line[6]} = 1;
	}
}

open(FILE2, '<', $ARGV[1]) or die("open $ARGV[1] error.");
my %fam2;
while(<FILE2>)
{
	chomp;
	my @line = split(',');
	if(defined( $fam2{$line[6]} ))
	{
		$fam2{$line[6]} = $fam2{$line[6]} + 1;
	}
	else
	{
		$fam2{$line[6]} = 1;
	}
}

my $index = 1;
for my $i (keys %fam1)
{
	if(defined($fam2{$i}))
	{
		my $d = abs($fam1{$i} - $fam2{$i});

## TO USE DCJ, PLEASE UNCOMMENT THE FOLLOWING FOUR LINES
#		if(!($fam1{$i} eq $fam2{$i}))
#		{
			next;
#		}

		print("$i $index\n");
		$index = $index + 1;

	}
}
