#!/bin/bash

if [ $# -lt 1 ]; then echo "usage: $0 pdb"; exit 1; fi
xdir=$(dirname $0)

for x in $*; do
	n=$(basename $x .pdb)
	[ -f $n.asa ] && continue
	$xdir/pdb2xrn.py $x |$xdir/calASA STDIN > $n.asa
done
