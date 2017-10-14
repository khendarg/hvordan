#!/usr/bin/env python2
from __future__ import print_function

class Node(object):
	def __init__(self, weight, symbol):
		self.left = None
		self.right = None
		self.parent = None
		self.weight = weight
		self.symbol = symbol

	def __iter__(self):
		out = []
		if self.left is not None: out.append(self.left)
		if self.right is not None: out.append(self.right)

	def pretty_print(self, indent=0):
		out = '    ' * indent
		out += '"%s":%s\n' % (self.symbol, self.weight)
		if self.left is not None: out += self.left.pretty_print(indent+1)
		if self.right is not None: out += self.right.pretty_print(indent+1)

		if indent == 0: return out.strip()
		else: return out

	def __str__(self): return self.pretty_print()

	def __lt__(self, other):
		if self.weight < other.weight: return True
		elif self.weight == other.weight:
			if self.symbol < other.symbol: return True
			else: return False
		else: return False

	def __gt__(self, other):
		if self.weight > other.weight: return True
		elif self.weight == other.weight:
			if self.symbol > other.symbol: return True
			else: return False
		else: return False

	def find(self, target, path=tuple()):
		p = None
		if self.symbol == target: return path

		if self.left is not None:
			p = self.left.find(target, path+(0,))
			if p is not None: return p
		if self.right is not None:
			p = self.right.find(target, path+(1,))
			if p is not None: return p
		if p is None: return None

class Huffman(object):
	def __init__(self, raw):
		self.raw = raw
		self.compressed = None

	def compress(self):
		symbols = {}
		for x in self.raw:
			try: symbols[x] += 1
			except KeyError: symbols[x] = 1

		rawqueue = sorted([Node(symbols[s], s) for s in symbols])
		cpxqueue = []

		def min_raw_cpx(raw, cpx):
			if raw and cpx:
				if raw[0] < cpx[0]: return raw.pop(0)
				else: return cpx.pop(0)
			elif raw: return raw.pop(0)
			elif cpx: return cpx.pop(0)
			return None

		parent = None
		while rawqueue or cpxqueue: 
			#left = minimum between raw and cpx
			#right = minimum between raw and cpx
			#parent = 
			left = min_raw_cpx(rawqueue, cpxqueue)
			right = min_raw_cpx(rawqueue, cpxqueue)

			if right == None: 
				if parent is None: 
					parent = Node(left.weight, left.symbol)
					parent.left = left
					left.parent = parent
				break

			parent = Node(left.weight+right.weight, left.symbol+right.symbol)

			parent.left = left
			parent.right = right
			left.parent = parent
			right.parent = parent

			cpxqueue.append(parent)

		try: self.parent = parent
		except: 
			print(self.raw)
			print(sorted([Node(symbols[s], s) for s in symbols]))

		dictionary = {}
		for x in symbols: dictionary[x] = parent.find(x)

		out = tuple()
		for x in self.raw: out = out + dictionary[x]

		#TODO: finish compressor and implement decompressor

		return 1. * len(out) / len(self.raw)

def seq_entropy(fasta):
	import re
	if type(fasta) is str: f = fasta.split('\n')
	else: f = fasta

	name = None
	entropy = None
	sequence = ''

	entropies = []
	for l in f:
		if l.startswith('>'):
			if sequence:
				h = Huffman(re.sub('[^ACDEFGHIKLMNPQRSTVWXY]', '', sequence))
				#print('%s\t%0.5f\t%d' % (name, entropy, len(sequence)))
				#print('%d\t%0.3f' % (len(sequence), time.time()-last))
				entropies.append(h.compress())
				sequence = ''
				#del h
			name = l.split()[0]
		else: sequence += l.strip()

	if sequence:
		h = Huffman(re.sub('[^ACDEFGHIKLMNPQRSTVWXY]', '', sequence))
		entropies.append(h.compress())
		#print('%s\t%0.5f\t%d' % (name, entropy, len(sequence)))
	return entropies
		
		

if __name__ == '__main__':
	#h = Huffman('A_DEAD_DAD_CEDED_A_BAD_BABE_A_BEADED_ABACA_BED')

	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='fasta to inspect')
	args = parser.parse_args()

	with open(args.infile) as f: seq_entropy(f)
