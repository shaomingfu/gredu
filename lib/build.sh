#!/bin/bash

libs="gbase"
dir=`pwd`

for i in `echo $libs`
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
