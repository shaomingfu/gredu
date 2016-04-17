#!/bin/bash

subs="rmtandem"
dir=`pwd`

for i in `echo $subs`
do
	cur="$dir/$i"
	echo $cur
	cd $cur
	aclocal
	autoconf
	autoheader
	automake -a
	./configure
	make
	cd $dir
done
