#!/usr/bin/env python2

from __future__ import print_function
import os, sys, subprocess, re 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time

import numpy as np

import warnings

import quod, tcblast

DEBUG = 1
VERBOSITY = 0
def warn(*msgs):
	for l in msgs: print('[WARNING]', l, file=sys.stderr)
def error(*msgs):
	for l in msgs: print('[ERROR]', l, file=sys.stderr)
	exit()
def info(*msgs):
	for l in msgs: print(l, file=sys.stderr)

#
#def retrieve_full(accession):
#	faa = ''
#	if DEBUG:
#		if VERBOSITY:
#			warnings.warn('Running in debug mode, returning predetermined sequence')
#		faa += '>gnl|TC-DB|373567223|9.B.143.4.7 protein of unknown function DUF1275 [Methylobacterium extorquens DSM 13060]\nMSDEPVFRPRTRLPFEERLNDRADPVTRPWQIGFGVVLTALAGFVDALGFIRLGGLYTSL\nMSGNTTQLAVALGHGEPLGAVLPALLIGAFLVGAVSGGAIAALCPPRWVTPAVLGLEAVA\nLTAAVALAAEHAHVGVASLFLALAMGGQNAVLAHVQGFRAGTTFVTGALFAFGQKAALAL\nAGRGPRLGWVGDGSVWLSLLVGAVAGTLAHAHLGIAALAIPAVVAAGLCLAATLFTLVYR\nRAPVTLTKI'
#	if not DEBUG:
#		faa += subprocess.check_output(['blastdbcmd', '-db', 'nr', '-entry', accession, '-target_only'])
#	return faa

def retrieve_full(accession, p2dir='./'):
	#retrieves full Bs or Cs
	#if VERBOSITY: warn('I left a clumsy hack in. If you are /not/ just testing, contact me.')
	foundit = 0
	seq = ''
	f = open(p2dir + '/subjects.faa')
	for l in f: 
		if foundit and l.startswith('>'): break
		elif foundit: seq += l.strip()
		if accession in l: 
			foundit = 1
			seq += l
	f.close()
	if seq: return seq
	f = open(p2dir + '/targets.faa')
	for l in f: 
		if foundit and l.startswith('>'): break
		elif foundit: seq += l.strip()
		if accession in l: 
			foundit = 1
			seq += l
	f.close()
	return seq

def retrieve_mother(accession, p2dir='./'):
	#retrieves full As or Ds

	allgsats = os.listdir(p2dir + '/gsat')

	bestm = ''
	bests = 0
	for fn in allgsats: 
		if accession in fn and re.findall('[0-9]\.[A-Z]\.[0-9]+', fn): 
			m = (fn[fn.find('.')+1:fn.find(accession)-1])
			f = open(p2dir + '/gsat/' + fn)
			score = None
			for l in f: 
				if l.startswith('Precise score'): score = float(l.split()[-1])
			f.close()

			if score > bests:
				bestm = m
				bests = score
	print(accession, bestm)

	#def relevant(x):
	#	#return x[0]
	#	y = x[x.find('.')+1:-5]
	#	child = y.split('.')[-1]
	#	mother = y[:-len(child)-1]
	#	return mother, child

	#seq = '>' + accession + '\n'
	#for fn in os.listdir(gsatdir):
	#	if fn.startswith('SRRSW_'): 
	#		if accession != relevant(fn)[1]: continue
	#		f = open(gsatdir + '/' + fn)
	#		for l in f: 
	#			if l.startswith('B_Sequence'): seq += re.sub('[^A-Z]', '', l[22:71])
	#if seq.endswith('\n'): 
	#	print(accession)
	#	raw_input() #WP_007105673
	#else: return seq

def parse_report(report, cutoff=15., skip=1, p2dir='./'):
	#gets bc
	bc_data = []
	n = 0
	reportlen = len(report.split('\n'))
	for l in report.split('\n'):
		if not l.strip(): continue
		if l.strip().startswith('#'): continue
		if skip: 
			skip -= 1
			continue
		fields = l.split('\t')

		#subjID, targID, SS_Zscore, GSAT_Zscore, Subjlen, Targlen, Sseq, Tseq, TMO score
		#0,      1,      2,         3,           4,       5,       6,    7,    8
		if float(fields[3]) < cutoff: continue

		#OPTIMIZE: grab all sequences in one fell swoop instead of closing and reopening constantly
		full_sseq = retrieve_full(fields[0], p2dir=p2dir)
		full_tseq = retrieve_full(fields[1], p2dir=p2dir)
		aln_sseq = '>%s_aln\n%s' % (fields[0], fields[6])
		aln_tseq = '>%s_aln\n%s' % (fields[1], fields[7])

		#print(fields)
		#raw_input()

		#may not be necessary
		#part_sseq = aln_sseq.sub('-', '')
		#part_tseq = aln_tseq.sub('-', '')

		bc_data.append([full_sseq, full_tseq, aln_sseq, aln_tseq])
		n += 1
		if VERBOSITY: info('Processed %d/%d total records with a cutoff of %0.2f' % (n, reportlen, cutoff))
	return bc_data

def sanitize(putative_filename, tiny=0):
	if tiny: putative_filename = putative_filename.split()[0]
	fn = re.sub('[>/]', '', putative_filename)
	fn = re.sub('[^A-Za-z0-9\.]', '_', fn)
	return fn
def relative(fn):
	if fn.startswith('/'): return fn
	else: return fn[fn.find('/')+1:]

def quod_them(bc_data, outdir='warum_out', overwrite=0, resoln=100):
	if not os.path.isdir(outdir): os.mkdir(outdir)
	if not os.path.isdir(outdir + '/images'): 
		os.mkdir(outdir + '/images')


	#graph A

	#graph B
	title = re.sub('\n.*', '', bc_data[0])
	bfilename = sanitize(title, tiny=1) + '.png'
	labels = [bc_data[0].split('\n')[0]]
	seq1 = bc_data[0].split('\n')[1]
	#os.system('quod -d %s/images -l "%s" -q -t png -r %d -o %s %s' % (outdir, title, resoln, bfilename, seq1))
	quod.what([seq1], title=title, labels=labels, imgfmt='png', directory=(outdir + '/images'), filename=bfilename, dpi=resoln, hide=1)

	#graph C
	title = re.sub('\n.*', '', bc_data[1])
	cfilename = sanitize(title, tiny=1) + '.png'
	labels = [bc_data[1].split('\n')[0]]
	seq1 = bc_data[1].split('\n')[1]
	#os.system('quod -d %s/images -l "%s" -q -t png -r %d -o %s %s' % (outdir, title, resoln, bfilename, seq1))
	quod.what([seq1], title=title, labels=labels, imgfmt='png', directory=(outdir + '/images'), filename=cfilename, dpi=resoln, hide=1)

	#graph D

	#graph BC
	title = re.sub('\n.*', '', bc_data[2]) + '(red) and ' + re.sub('\n.*', '', bc_data[3]) + ' (blue)'
	bcfilename = sanitize(title) + '.png'
	label1 = bc_data[2].split('\n')[0]
	label2 = bc_data[3].split('\n')[0]
	labels = [label1, label2]
	seq1 = bc_data[2].split('\n')[1]
	seq2 = bc_data[3].split('\n')[1]
	#print('quod -d %s -l "%s" -q -t png -r 300 %s %s' % (outdir, title, seq1, seq2))
	#os.system('quod -d %s/images -l "%s" -q -t png -r %d -o %s %s %s' % (outdir, title, resoln, bcfilename, seq1, seq2))
	quod.what([seq1, seq2], title=title, labels=labels, imgfmt='png', directory=(outdir + '/images'), filename=bcfilename, dpi=resoln, hide=1)

	return outdir + '/images/' + bfilename, outdir + '/images/' + cfilename, outdir + '/images/' + bcfilename
	

def rebuild_iter1(gsatdir='./gsat/'):
	def relevant(x):
		#return x[0]
		x = x[x.find('.')+1:-5]
		child = x.split('.')[-1]
		mother = x[:-len(child)-1]
		return mother, child
	#print(os.listdir(gsatdir))
	rawgsats = os.listdir(gsatdir)
	gsats = map(relevant, rawgsats)

	pair2fn = dict(zip(gsats, rawgsats))
	pair2score = {}
	for k in sorted(pair2fn.keys()):
		f = open(gsatdir + '/' + pair2fn[k])
		for l in f: 
			if l.startswith('Precise score'):
				pair2score[k] = float(l.split()[-1])
	child2mother = {}
	finalscores = {}
	
	for l in gsats:
		try:
			if pair2score[(l[1], l[0])] > pair2score[(child2mother[l[1]], l[1])]:
				child2mother[l[1]] = l[0]
				finalscores[l[1]] = pair2score[(l[0], l[1])]
		except KeyError:
			child2mother[l[1]] = l[0]
			finalscores[l[1]] = pair2score[(l[0], l[1])]
	#print(child2mother)
	return child2mother, finalscores

def build_html(bc, graphfns, blasts, outdir='warum_out', outfn='test.html', title='Unnamed report'):
	if not os.path.isfile(outdir + '/openclose.js'):
		f = open(outdir + '/openclose.js', 'w')
		f.write('function toggle_section(sectionid, selfid) {\n\tvar section = document.getElementById(sectionid);\n\tvar me = document.getElementById(selfid);\n\tconsole.log([section, section.style.display]);\n\tif (section.style.display == \'none\') {\n\t\tsection.style.display = \'block\';\n\t\tme.innerHTML = \'Hide\';\n\t} else { \n\t\tsection.style.display = \'none\'; \n\t\tme.innerHTML = \'Show\';\n\t}\n}')
		f.close()
	if not os.path.isfile(outdir + '/nice.css'):
		f = open(outdir + '/nice.css', 'w')
		f.write('body {\n\tfont-family: sans;\n\theight: 100%;\n}\ndiv {\n\tdisplay: block;\n}\ndiv.tcblast {\n\tmax-width: 1500px;\n}\ndiv.fullblast {\n\twidth: 50%;\n\tfloat: left;\n}\ndiv.tabular1 {\n\twidth: 50%;\n\tfloat: left;\n\theight: 100%;\n}\ndiv.tabular2 {\n\twidth: 50%;\n\tfloat: right;\n\theight: 100%;\n}\nimg.bluebarplot {\n\tmax-width: 100%;\n\theight: auto;\n}\n.clear { clear: both; }\n.scrollable {\n\toverflow-y: scroll;\n}\n.resizeable {\n\tresize: vertical;\n\toverflow: auto;\n\tborder-bottom: 1px solid gray;\n\tdisplay: block;\n}\n.bluebars {\n\theight: 25vh;\n}\n.pairwise {\n\theight: 50vh;\n}\n.whatall {\n\theight: 50vh;\n}\n.whataln {\n\twidth: 100%;\n}\n#seqs {\n\tdisplay: none;\n}\n\n\n\n.summtbl {\n\tfont-family: monospace, courier;\n\tfont-size: 75%;\n}\n.oddrow {\n\tbackground-color: #d8d8d8;\n}\ntd {\n\tpadding-right: 1em;\n}\n.red {\n\tcolor: red;\n}\nimg {\n\tborder: 1pt solid black;\n}\n')
		f.close()
	out = '<html><head><title>%s</title>' % title
	out += '\n<link rel="stylesheet" type="text/css" href="nice.css"/>'
	out += '\n<script src="openclose.js"></script>'
	out += '\n</head><body>'

	out += '\n<h1>%s</h1>' % title
	out += '\n<h2>Table of contents</h2>'
	out += '\n<button class="showhide" id="tocsh" onclick="toggle_section(\'toc\', \'tocsh\')">Hide</button>'

	out += '\n<div class="toc" id="toc"> <ol> <li><a href="#summary">Summary</a></li> <li><a href="#pairwise">Pairwise</a></li> <li><a href="#abcd">ABCD hydropathy plots</a></li> <li><a href="#bc">BC hydropathy plot</a></li> </ol> </div>'

	out += '\n<h2>TCBLAST</h2>'

	#bluebars
	out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'
	out += '\n<div class="tcblast" id="tcblast"><a name="summary"><h3>Summary</h3></a>'
	out += '\n<div class="resizeable bluebars"><div class="scrollable tabular1">'
	out += '\n%s' % blasts[0][0]
	out += '\n</div><div class="scrollable tabular2">'
	out += '\n%s' % blasts[1][0]
	out += '\n</div></div>'

	#pairwise
	out += '\n<a name="pairwise"><h3 class="clear">Pairwise</h3></a><div class="resizeable pairwise"><div class="scrollable tabular1">'
	out += '\n%s' % blasts[0][1]
	out += '</div><div class="scrollable tabular2">'
	out += '\n%s' % blasts[1][1]
	out += '\n</div></div></div>'

	#abcd bc
	out += '\n<a name="abcd"><h3 class="clear">ABCD Hydropathy plots</h3></a>'
	out += '\n<button class="showhide" id="abcdsh" onclick="toggle_section(\'abcd\', \'abcdsh\')">Hide</button>'

	out += '\n<div class="whatall" id="abcd">'
	out += '\n<div class="tabular1">'
	out += '\nA<br/><img class="bluebarplot" id="plota" src="%s"/><br/>'
	out += '\nB<br/><img class="bluebarplot" id="plotb" src="%s"/><br/>' % relative(graphfns[0])
	out += '\n</div><div class="tabular2">'
	out += '\nC<br/><img class="bluebarplot" id="plotc" src="%s"/><br/>' % relative(graphfns[1])
	out += '\nD<br/><img class="bluebarplot" id="plotd" src="%s"/><br/>'
	out += '\n</div></div>'

	out += '\n<a name="bc"><h3 class="clear">BC hydropathy plot</h3></a>'
	out += '\n<button class="showhide" id="bcsh" onclick="toggle_section(\'bc\', \'bcsh\')">Hide</button>'

	out += '\n<div class="resizeable whataln" id="bc"><div class="scrollable">'
	out += '<img id="plotbc" src="%s"/><br/>' % relative(graphfns[2])
	out += '\n</div></div>'

	out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'


	out += '\n<br/><div style="height: 10ex"></div>'
	out += '\n</body></html>'

	#print(graphfns)
	f = open(outdir + '/' + outfn, 'w')
	f.write(out)
	f.close()

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='saves your wrist!')

	parser.add_argument('-d', default='./', help='Protocol2 directory {default: ./}')
	parser.add_argument('-r', default=100, type=int, help='Resolution in dpi of graphs {default:100}')

	#tested for 2.A.1 vs 2.A.6, cutoff of 25, B, C, BC
	#approximate details when launching quod as a process: 
	#12dpi: 3K / 1.5s
	#24dpi: 7K / 1.5s
	#50dpi: 20KB / 1.5s
	#100dpi: 40KB / 1.7s
	#300dpi: 180KB / 2.9s
	#600dpi: 350KB / 7s

	#approximate details when calling from quod:
	#12dpi: 3K / 0.7s
	#24dpi: 7K / 0.8s
	#50dpi: 20KB / 0.8s
	#100dpi: 40KB / 0.9s
	#300dpi: 180KB / 2.1s
	#600dpi: 350KB / 62s

	parser.add_argument('-z', default=15., type=float, help='Z-score cutoff {default: 12}')
	parser.add_argument('-o', default='warum_out', help='Output directory {default: warum_out}')

	args = parser.parse_args()

	if not os.path.isdir(args.d): 
		error('Could not find directory %s' % args.d)

	if not os.path.isfile(args.d + '/' + 'report.tbl'): 
		error('Could not find file %s/%s' % (args.d, 'report.tbl'))

	if not os.path.isdir(args.o): os.mkdir(args.o)

	c2m, scores = rebuild_iter1(gsatdir=args.d + '/gsat/')

	f = open(args.d + '/report.tbl')
	report = f.read()
	bc = parse_report(report, cutoff=args.z, p2dir=args.d)

	n = 0
	A = {}
	D = {}
	titles = []
	for l in bc: 
		start = time.time()
		a = (l[0][1:l[0].find(' ')].strip())
		d = (l[1][1:l[1].find(' ')].strip())

		if a not in A: A[a] = retrieve_mother(a, p2dir=args.d)
		if d not in D: D[d] = retrieve_mother(d, p2dir=args.d)

		graphfns = quod_them(l, resoln=args.r)
		n += 1
		info('Built a report (%d/%d) in %0.2fms' % (n, len(bc), (time.time()-start)*1000))

		accs = []

		for x in l:
			#standard FASTA
			if '|' in x.split('\n')[0]: accs.append(x.split('|')[3])
			#simplified TC-FASTA
			elif '-' in x.split('\n')[0]: accs.append(x.split('-')[1].split()[0])
			#worst case
			else: accs.append(sanitize(x.split('\n')[0]))

		if not os.path.isdir(args.o + '/images'): os.mkdir(args.o + '/images')

		blasts = [tcblast.til_warum(l[0], args.o + '/images/' + accs[0] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop'), tcblast.til_warum(l[1], args.o + '/images/' + accs[1] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop')]

		build_html(l, graphfns, blasts, outdir=args.o, outfn='%s__%s.html' % tuple(accs)[0:2], title='%s vs %s' % tuple(accs)[0:2])
		
		#blasts.append(tcblast.til_warum(
		#blasts.append(
		#print(l[0])
		#print(l[1])
		#tcblast.til_warum(l[0], args.o + '/images/' + 

	#graphs

	#sequences!
	#for l in bc: 
	#	print(l)
