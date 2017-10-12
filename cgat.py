#!/usr/bin/env python2
#CGAT: The Cludgy Global Alignment Tool
#Making knock-offs of niche software since 2016!
from __future__ import print_function, division, generators

import Bio.Seq
import Bio.SeqIO
import Bio.Alphabet
import subprocess
import random
import tempfile, os, re, sys

VERBOSITY = 0

class Seq(Bio.Seq.Seq):
	def __init__(self, *args, **kwargs):
		Bio.Seq.Seq.__init__(self, *args, **kwargs)


def info(*text):
	for line in text: print('[INFO]: %s' % line, file=sys.stderr)

def cgat(fasta1, fasta2, gapopen=10.0, gapextend=0.5, outfile=None, shuffles=10, shuffleseq=1, silent=False):
	tape1 = Bio.SeqIO.parse(fasta1, 'fasta')
	tape2 = Bio.SeqIO.parse(fasta2, 'fasta')

	run = 1
	def fastafy(title, seq):
		return '>%s\n%s\n' % (title.strip(), seq.strip())
	summary = []
	while run:
		try:
			rec1, rec2 = tape1.next(), tape2.next()
			#seq1 = []
			#seq2 = []
			fas1, fas2 = '', ''
			shuf1, shuf2 = list(rec1.seq), list(rec2.seq)
			testscore = 0.0
			try:
				f1 = tempfile.NamedTemporaryFile(delete=False)
				f2 = tempfile.NamedTemporaryFile(delete=False)
				f1.write(fastafy(rec1.name, rec1.seq)), f1.close()
				f2.write(fastafy(rec2.name, rec2.seq)), f2.close()
				cmd = ['needle', '-asequence', f1.name, '-bsequence', f2.name, '-gapopen', str(gapopen), '-gapextend', str(gapextend)]
				cmd += ['-outfile', 'stdout']
				p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				out, err = p.communicate(fas1 + chr(0) + fas2)
				for l in out.split('\n'):
					if '# Score' in l: testscore = float(l.split()[-1])
				if not silent: print(out.decode('utf-8'))
			finally:
				os.remove(f1.name)
				os.remove(f2.name)
			if VERBOSITY: info('Shuffling...')
			n = 10
			milestones = list(range(n))
			milestones = [x/n for x in milestones]
			if shuffleseq == 0:
				fas2 += fastafy('seq2', (''.join(shuf2)))
			elif shuffleseq == 1:
				fas1 += fastafy('seq1', (''.join(shuf1)))
			for i in range(shuffles):
				if shuffleseq == 0:
					random.shuffle(shuf1)
					fas1 += fastafy('seq1_s%d' % i, (''.join(shuf1)))
				elif shuffleseq == 1:
					random.shuffle(shuf2)
					fas2 += fastafy('seq2_s%d' % i, (''.join(shuf2)))
				else:
					random.shuffle(shuf1)
					fas1 += fastafy('seq1_s%d' % i, (''.join(shuf1)))
					random.shuffle(shuf2)
					fas2 += fastafy('seq2_s%d' % i, (''.join(shuf2)))
				if milestones and i/shuffles > milestones[0] and VERBOSITY: info('%0.1f%% done' % (100*milestones.pop(0)))
			scores = []
			try:
				f1 = tempfile.NamedTemporaryFile(delete=False)
				f2 = tempfile.NamedTemporaryFile(delete=False)
				f1.write(fas1), f1.close()
				f2.write(fas2), f2.close()
				cmd = ['needleall', '-asequence', f1.name, '-bsequence', f2.name, '-gapopen', str(gapopen), '-gapextend', str(gapextend)]
				if outfile: cmd += ['-outfile', oufile]
				else: cmd += ['-outfile', 'stdout']
				p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				if shuffleseq == 0 or shuffleseq == 1: alns = shuffles
				else: alns = shuffles**2
				if VERBOSITY: info('Performing %d alignments...' % (alns))
				out, err = p.communicate(fas1 + chr(0) + fas2)
				if VERBOSITY: info('Done!')
				for l in out.split('\n'):
					if l.strip() and '(' in l: scores.append(float(l.split()[-1][1:-1]))
				if shuffleseq == 0: 
					if len(scores) < fas1.count('>'): scores += [0.] * (fas1.count('>') - len(scores))
				elif shuffleseq == 1: 
					if len(scores) < fas2.count('>'): scores += [0.] * (fas2.count('>') - len(scores))
				#print(out.decode('utf-8'))
			finally:
				os.remove(f1.name)
				os.remove(f2.name)
			mean = sum(scores)/len(scores)
			stdev = (sum([(x - mean)**2 for x in scores])/(len(scores)-1))**.5
			if not silent:
				print('============ FINISHED =============')
				print('Shuffles: %d' % shuffles)
				print('Average Quality (AQ): %0.2f +/- %0.2f' % (mean, stdev))
				print('Test score: %0.1f' % testscore)
				if stdev:
					print('Standard score (Z): %0.1f' % ((testscore - mean)/stdev))
					print('Precise score (Z): %0.3f' % ((testscore - mean)/stdev))
				else:
					print('Standard score (Z): %0.1f' % ((-999.9)))
					print('Precise score (Z): %0.3f' % ((-999.999)))
			summary.append((fas1[:fas1.find('\n')-1], fas2[:fas2.find('\n')-1], stdev))
			try: return ((testscore - mean)/stdev)
			except ZeroDivisionError: return -999.999

		except StopIteration: run = 0

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta1', help='this sequence is not shuffled')
	parser.add_argument('fasta2', help='this sequence is shuffled')
	parser.add_argument('-g', '--gapopen', type=float, default=10., help='gap open penalty {default:10}')
	parser.add_argument('-e', '--gapextend', type=float, default=0.5, help='gap extension penalty {default:0.5}')
	parser.add_argument('-o', '--outfile', default=None, help='where to dump alignment output {default:stdout}')
	parser.add_argument('-s', '--shuffles', type=int, default=10, help='number of shuffles {default:10}')
	parser.add_argument('-v', action='store_true', help='verbose output')

	args = parser.parse_args()

	if args.v: VERBOSITY = 1

	cgat(args.fasta1, args.fasta2, gapopen=args.gapopen, gapextend=args.gapextend, outfile=args.outfile, shuffles=args.shuffles)
