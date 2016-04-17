#!/usr/bin/perl

if($#ARGV ne 3)
{
	die("usage: build.pl human.list mouse.list human.input mouse.input\n");
}

my $tmp_human_file = "./gredo.tmp.human";
my $tmp_mouse_file = "./gredo.tmp.mouse";
my $tmp_index_file = "./gredo.tmp.index";

`cat $ARGV[0] | grep -v "Strand,Transcript Start" | sort -k1,1 -t "," > $ARGV[2]`;
`cat $ARGV[1] | grep -v "Strand,Transcript Start" | sort -k1,1 -t "," > $ARGV[3]`;
`./longest.pl $ARGV[2] | grep ENS | sort -k3,3 -k5,5n -t "," > $tmp_human_file`;
`./longest.pl $ARGV[3] | grep ENS | sort -k3,3 -k5,5n -t "," > $tmp_mouse_file`;
`./family.pl $tmp_human_file $tmp_mouse_file > $tmp_index_file`;
`./join.pl $tmp_index_file $tmp_human_file > $ARGV[2]`;
`./join.pl $tmp_index_file $tmp_mouse_file > $ARGV[3]`;

`rm -rf $tmp_human_file $tmp_mouse_file $tmp_index_file`;
