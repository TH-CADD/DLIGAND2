#!/usr/bin/env python
from misc_yang import *

def rdvdw(fn):
	dict1 = {}
	for x in openfile(fn):
		if x.startswith('TYPE:'): list_type = x.split()[1:]
		elif x.startswith('RADII:'): list_rad = x.split()[1:]
		elif x.startswith('ATOM '):
			ss = x.split()
			dict1[ss[1]+'-'+ss[2]] = ss[3]
	dict_t2r = dict(zip(list_type, list_rad))
#
	dict_a2r = {}
	for x,y in dict1.items():
		dict_a2r[x] = dict_t2r[y]
	return dict_a2r

if __name__ == '__main__':
	if len(argv)<2: die('usage: RUN pdb [vdw]')

	if len(argv) > 2: vdwfile = argv[2]
	else:
		for x in ('.', os.path.expanduser('~/source/lib'), os.path.dirname(argv[1])):
			f1 = os.path.join(x, 'aminorna.mol2')
			if  os.path.isfile(f1):
				vdwfile = f1; break
		else: die('not found vdwfile: aminorna.mol2')
	dict_a2r = rdvdw(vdwfile)
#
	last_str = ''
	for x in openfile(argv[1]):
		if x.startswith('ATOM '):
			if last_str == x[13:27]: continue  # to skip some dup. lines
			else: last_str = x[13:27]
			rn = x[17:20].strip()
			if rn == 'DU': rn = 'U'
			if rn == 'N': rn = 'A'
			an = x[12:16].strip()
			if x[16] not in ' A':
				print >>stderr, 'skipping', x[16], x,
				continue
			ch1 = x[21]
			ridx = x[22:27].strip().replace(' ', '-')
			
			if an[0]=='H' or x[13]=='H': continue
			print ' '.join([x[30+k*8: 38+k*8] for k in range(3)]),
			a1 = '-'.join([rn,an])
			if a1 in dict_a2r: print dict_a2r[a1],
			else:
				print >>stderr, 'not found res/atom type: ', rn, an
				print 0,

			print ch1, '_'.join([an, rn, ridx])

