#!/usr/bin/env python2
from __future__ import division, print_function

import os, re

def build_index(prefix='./'):
	#builds an accession->directory index
	#takes a famXpander directory
	fams = os.listdir(prefix)
	index = {}

	BREAKME = 1000

	for fam in fams:
		f = open('%s/%s/psiblast.tbl' % (prefix, fam))
		x = f.read()
		f.close()

		for l in x.split('\n'):
			if re.match('\s*#', l): continue
			if not l.strip(): continue
			ls = l.split()

			try: index[ls[1]].add(ls[0])
			except KeyError: index[ls[1]] = {ls[0]}

			BREAKME -= 1
			if not BREAKME: break
		if not BREAKME: break

	print(index)
	for i in index: print(i, index[i])
	return index

def rev_iteration(prot2dir, prot1dir):
	f = open(prot2dir + '/report.tbl')
	x = f.read()
	f.close()

	lineno = 0
	STOPAT = 100

	accs, acct = [], []

	for l in x.split('\n'):
		if lineno == 1: 
			lineno += 1
			continue
		if lineno == STOPAT: break
		elif lineno == 0: fams, famt = l.split()[-3], l.split()[-1]
		elif l.strip(): 
			#filter for z-scores or something around here
			accs.append(l.split()[0])
			acct.append(l.split()[1])
		lineno += 1

	fs = open(prot1dir + '/' + fams + '/psiblast.tbl')
	ft = open(prot1dir + '/' + famt + '/psiblast.tbl')
	s, t = fs.read(),  ft.read()
	fs.close(), ft.close()

	raw_input('read the file! continue?')
	for l in s.split('\n'):
		print(l)


def search_gene(acc, prefix='./', preindex=None):
	pass

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser()

	parser.add_argument('-i', help='root directory, contains famXpander and all protocol2 results')

	famxpanderdir = '/ResearchData/Users/amedrano/MFS/famXpander/'
	prot2dir = '/ResearchData/Users/amedrano/MFS/PositiveControl/2.A.1_vs_2.A.100/2.A.1_vs_2.A.100/'
	rev_iteration(prot2dir, famxpanderdir)
