#!/bin/bash

xdir=$(dirname $0)

if [ $# -eq 0 ]; then echo "$0 pdbs.."; exit 1; fi
for m in $*; do
	if [ ! -f $m ]; then echo "file not exists: $m"; exit 1; fi
	n=$(basename $m .pdb)
	[ -f $n.ASA ] && continue
	$xdir/pdb2xrn.py $m |$xdir/calASA STDIN > $n.ASA
done
