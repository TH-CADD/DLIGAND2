#!/usr/bin/env python
from misc_yang import *

def calASA(fn):
	list_rn, list_an = [], []
	for x in open(fn):
		if x[0] == '#':
			print x,
			continue

		ss = x.split()
		ss2 = ss[0].split('_')
		r1 = '_'.join(ss2[1:3] + [ss[1]])
		if len(list_rn)==0 or r1!=list_rn[-1]:
			list_rn.append(r1)
			list_an.append([])
		list_an[-1].append((ss2[0], ss[2]))
	
	for x,yss in zip(list_rn, list_an):
		print x,
		print '%g' % sum([float(ys[1]) for ys in yss]), len(yss),
		ys2 = [float(ys[1]) for ys in yss if not ys[0].endswith("'")]
		print '%g' % sum(ys2), len(ys2)


fn = argv[1]
calASA(fn)
