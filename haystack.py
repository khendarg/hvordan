#!/usr/bin/env python2
from __future__ import print_function, division
import subprocess, sys, re, os, tempfile
import cgat

VERBOSITY = 0

def info(*text):
	for line in text: print('[INFO]:', line, file=sys.stderr)
def warn(*text):
	for line in text: print('[WARNING]:', line, file=sys.stderr)
def error(*text):
	for line in text: print('[ERROR]:', line, file=sys.stderr)
	exit()

class Seq(object):
	def __init__(self, fasta): 
		self.header = fasta[1:fasta.find('\n')-1]
		self.seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', fasta[fasta.find('\n')+1:])
		self.fasta = '>' + self.header + '\n' + self.seq

		self.tmss = None

	def hmmtop(self, extend=5):
		if self.tmss: return self.tmss
		else:
			self.tmss = []
			p = subprocess.Popen(['hmmtop', '-if=--', '-sf=FAS', '-pi=spred', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate(input=self.fasta)
			print(err.strip(), file=sys.stderr)

			indices = re.findall('((?:\s*[0-9]\s*)+)$', out)[0].strip().split()
			for i in range(1, len(indices), 2): 
				self.tmss.append([int(indices[i])-extend, int(indices[i+1])+extend])

	def get_tms_sequence(self, n, extend=5): 
		if self.tmss is None: self.hmmtop(extend=extend)
		if self.tmss: return self.seq[self.tmss[n][0]-1:self.tmss[n][1]]
		else: return ''


	def get_tms_fasta(self, n, extend=5):
		if self.tmss is None: self.hmmtop(extend=extend)
		if self.tmss: 
			header = '>' + self.header + ' TMS %d' % (n+1)
			seq = self.seq[self.tmss[n][0]-1:self.tmss[n][1]]
			return header + '\n' + seq
		else: return ''

	def get_tms_seq(self, n, extend=5):
		if self.tmss is None: self.hmmtop(extend=5)
		if self.tmss: 
			header = '>' + self.header + ' TMS %d' % (n+1)
			seq = self.seq[self.tmss[n][0]-1:self.tmss[n][1]]
			return Seq(header + '\n' + seq)
		else: return None

	def cgat(self, other, gapopen=10.0, gapextend=0.5, shuffles=100, silent=True): 
		z = None
		try:
			f1 = tempfile.NamedTemporaryFile(delete=False)
			f2 = tempfile.NamedTemporaryFile(delete=False)

			f1.write('>' + self.header + '\n' + self.seq)
			f2.write('>' + other.header + '\n' + other.seq)
			f1.close(), f2.close()

			z = cgat.cgat(f1.name, f2.name, gapopen=gapopen, gapextend=gapextend, shuffles=shuffles, silent=silent)
		finally:
			os.remove(f1.name)
			os.remove(f2.name)
			return z

def printmatrix(m, places=1):
	out = ''
	for row in m:
		for col in row:
			preout = '\t%%0.%df' % places
			out += preout % col
		out += '\n'
	print(out)
	return out

class SmithWaterman:
	def __init__(self, seq1, seq2):
		self.seq1 = seq1
		self.seq2 = seq2
		self.matrix = None

	#def build_scoring_matrix(self, tmgap=0, gapopen=10., gapextend=0.5, shuffles=2000):
	def build_scoring_matrix(self, tmgap=2, gapopen=10., gapextend=0.5, shuffles=2000, extend=5):
		self.shuffles = shuffles
		self.tmgap = tmgap
		self.gapopen = gapopen
		self.gapextend = gapextend
		self.extend = extend

		if self.matrix is not None: return self.matrix

		self.seq1.hmmtop(extend=extend)
		self.seq2.hmmtop(extend)

		self.simil = []
		for i in range(len(self.seq1.tmss)):
			tms1 = self.seq1.get_tms_seq(i)
			self.simil.append([])
			for j in range(len(self.seq2.tmss)):
				tms2 = self.seq2.get_tms_seq(j)
				z = tms1.cgat(tms2, gapopen=gapopen, gapextend=gapextend, shuffles=shuffles)
				self.simil[-1].append(z)

		scores = []
		for i in range(-1, len(self.seq1.tmss)):
			scores.append([])
			for j in range(-1, len(self.seq2.tmss)):
				if i < 0 or j < 0: scores[-1].append(0)
				#linear gap penalty
				#TODO: implement affine gap penalty
				else: scores[-1].append(max(
					scores[i-1][j-1] + self.simil[i-1][j-1], 
					scores[i-1][j] - tmgap, 
					scores[i][j-1] - tmgap,
					0))
		self.matrix = scores
		return scores

	def align(self, shuffles=200, gapopen=10.0, gapextend=0.5, tmgap=2.0, extend=5): 
		self.build_scoring_matrix(shuffles=shuffles, gapopen=gapopen, gapextend=gapextend, tmgap=tmgap, extend=extend)

		highest = (0, 0, 0) #value, i, j
		for i, row in enumerate(self.matrix):
			for j, col in enumerate(row):
				if col > highest[0]: highest = col, i, j
				elif col == highest[0]:
					if (i+j) > (highest[1]+highest[2]): highest = col, i, j

		run = 1
		i, j = highest[1:3]
		path = []
		while run:
			path.append([i,j])
			gapup = self.matrix[i-1][j], [i-1, j]
			gapleft = self.matrix[i][j-1], [i, j-1]
			back = self.matrix[i-1][j-1], [i-1, j-1]

			best = max(gapup[0], gapleft[0], back[0])

			if best == back[0]: i, j = back[1]
			elif best == gapup[0]: i, j = gapup[1]
			elif best == gapleft[0]: i, j = gapleft[1]

			if i == 0 or j == 0 or self.matrix[i][j] == 0: run = 0
		path.reverse()
		self.path = [x[:] for x in path]

		#take care of gaps by selecting the best-scoring bit
		alignment = []
		for n in range(len(path)):
			if n < 1:
				pair = path[n]
				prev = pair[:]
			elif n >= 1:
				if path[n][0] == prev[0]:
					nov = self.matrix[path[n][0]][path[n][1]]
					vet = self.matrix[prev[0]][prev[1]]
					if nov > vet: 
						path[n-1][0] = None
						pair = path[n]
					elif vet > nov: 
						pair = [None, path[n][1]]
				elif path[n][1] == prev[1]: 
					nov = self.matrix[path[n][0]][path[n][1]]
					vet = self.matrix[prev[0]][prev[1]]
					if nov > vet: 
						path[n-1][1] = None
						pair = path[n]
					elif vet > nov: 
						pair = [path[n][0], None]
				else: pair = path[n]
			alignment.append(pair)
			prev = pair[:]

		#collapse reciprocal gaps
		removeme = []
		for n in range(len(alignment)):
			if n < 1: continue
			else: 
				if alignment[n] == alignment[n-1][::-1]:
					if alignment[n-1][0] is not None: alignment[n-1][1] = alignment[n-1][0]
					elif alignment[n-1][1] is not None: alignment[n-1][0] = alignment[n-1][1]
					removeme.append(n)
		for n in reversed(sorted(removeme)): alignment.pop(n)

		self.alignment = alignment

	def print_path(self, fullthreshold=6., halfthreshold=2.):
		top, middle, bottom = '', '', ''
		for pair in self.alignment:
			if pair[0]: 
				if pair[0] < 10: top += ' '
				top += ' %d' % pair[0]
			else: top += '   '
			if pair[1]: 
				if pair[1] < 10: bottom += ' '
				bottom += ' %d' % pair[1]
			else: bottom += '   '
			if pair[0] and pair[1]: 
				if self.matrix[pair[0]][pair[1]] >= fullthreshold: middle += '  |'
				elif self.matrix[pair[0]][pair[1]] >= halfthreshold: middle += '  :'
				else: middle += '  .'
			else: middle += '   '

		#for pair in self.raw_path: print(self.matrix[pair[0]-1][pair[1]-1], end=', ')
		#print()

		print('########################################')
		print('# Program: haystack')
		print('########################################')
		print()
		print('#=======================================')
		print('#')

		print('# 1: %s' % self.seq1.header)
		print('# 2: %s' % self.seq2.header)
		print('# TMS extension: %d' % self.extend)
		print('# TMS gap penalty: %0.1f' % self.tmgap)
		print('# Shuffles: %d' % self.shuffles)
		print('#')

		print('# needle Gap_penalty: %0.1f' % self.gapopen)
		print('# needle Extend_penalty: %0.1f' % self.gapextend)
		print('#')

		s = 0
		print('# TMSs\tZ-score\tSmith-Waterman score')
		for pair in self.path:
			print('# %d-%d\t%0.3f\t%0.2f' % (tuple(pair) + (self.simil[pair[0]-1][pair[1]-1], self.matrix[pair[0]][pair[1]])))
			s += self.simil[pair[0]-1][pair[1]-1]
		print('# Average Z-score: %0.3f' % (s/len(self.path)))

		print('#')
		print('#')
		print('#=======================================')

		print(top)
		print(middle)
		print(bottom)

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()

	parser.add_argument('fasta1', help='first fasta')
	parser.add_argument('fasta2', help='second fasta')
	parser.add_argument('-s', type=int, default=2000, help='shuffles {default:2000}')
	parser.add_argument('-v', help='verbose output')
	parser.add_argument('-g', type=float, default=10.0, help='gap open cost {default:10.0}')
	parser.add_argument('-e', type=float, default=0.5, help='gap extension cost {default:0.5}')
	parser.add_argument('-t', type=float, default=1.0, help='TMS gap cost {default:1.0}')

	parser.add_argument('-c', type=int, default=5, help='extend TMSs by this many residues {default:5}')

	args = parser.parse_args()

	if args.v: VERBOSITY = 1

	with open(args.fasta1) as f: fas1 = Seq(f.read())
	with open(args.fasta2) as f: fas2 = Seq(f.read())

	SW = SmithWaterman(fas1, fas2)
	SW.align(shuffles=args.s, gapopen=args.g, gapextend=args.e, tmgap=args.t, extend=args.c)
	SW.print_path()
