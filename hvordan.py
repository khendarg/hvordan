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

import Bio.Entrez

DEBUG = 1
VERBOSITY = 1
def warn(*msgs):
	for l in msgs: print('[WARNING]', l, file=sys.stderr)
def error(*msgs):
	for l in msgs: print('[ERROR]', l, file=sys.stderr)
	exit()
def info(*msgs):
	for l in msgs: print('[INFO]', l, file=sys.stderr)

def fetch(accessions, email=None, db='protein'):
	if not accessions: return ''
	if db == 'tcdb':
		out = ''
		for acc in accessions:
			try: 
				fa = subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-target_only', '-entry', acc])
				out += fa + '\n'
			except ValueError: raise ValueError
		return out
	else:
		if not email:
			if 'ENTREZ_EMAIL' in os.environ: email = os.environ['ENTREZ_EMAIL']
			else: 
				raise TypeError('Missing argument email')
		Bio.Entrez.tool = 'biopython'
		Bio.Entrez.email = email
		
		acclist = ''
		for x in accessions: acclist += ',' + x
		acclist = acclist[1:]

		try: 
			f = Bio.Entrez.efetch(db=db, id=acclist, rettype='fasta', retmode='text')
			out = f.read()
			f.close()
		except Bio.Entrez.urllib2.URLError: 
			out = subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, acclist), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
		return out

def parse_p2report(p2report, minz=15, maxz=None):
	line = 0

	if minz == None: minz = -2**16
	if maxz == None: maxz = 2**16

	bcs = []
	alnregs = {}
	for l in p2report.split('\n'):
		line += 1
		if line == 1:
			fams = [re.split('[, :]+', l)[2], re.split('[, :]+', l)[4]]
			bcs = {fams[0]:[], fams[1]:[]}
		elif line == 2: pass
		else: 
			if not l.strip(): continue
			ls = l.split('\t')
			z = float(ls[3])
			if minz <= z <= maxz: 
				#bcs.append(ls[:2])
				bcs[fams[0]].append(ls[0])
				bcs[fams[1]].append(ls[1])
				try: alnregs[ls[0]][ls[1]] = [ls[6], ls[7]]
				except KeyError: alnregs[ls[0]] = {ls[1]:[ls[6], ls[7]]}
			
	return fams, bcs, alnregs

def seek_initial(p1d, bcs):
	hits = {}
	for fam in sorted(bcs.keys()):
		hits[fam] = {}
		for bc in sorted(bcs[fam]): hits[fam][bc] = []
		f = open(p1d + '/%s.tbl' % fam)
		x = f.read()
		f.close()
		for l in x.split('\n'):
			if not l.strip(): continue
			if l.lstrip().startswith('#'): continue
			ls = l.split('\t')
			try: hits[fam][ls[1]].append((float(ls[4]), ls[0]))
			#except KeyError: hits[fam][ls[1]] = [(float(ls[4]), ls[0])]
			except KeyError: pass
		for bc in sorted(hits[fam]):
			hits[fam][bc] = sorted(hits[fam][bc])[-1]
	return hits

def clean_fetch(accs, outdir, force=False):

	if not os.path.isdir(outdir): os.mkdir(outdir)

	dlme = []
	tcdlme = []
	for acc in accs:
		if os.path.isfile(outdir + '/%s.fa' % acc): continue
		else: 
			if re.match('[0-9]\.[A-Z]\.[0-9]+\.', acc): tcdlme.append(acc)
			else: dlme.append(acc)

	allfaa = ''
	if dlme: 
		#if VERBOSITY: info('Downloading %s' % dlme)
		allfaa += fetch(dlme)
	if tcdlme:
		#if VERBOSITY: info('Loading %s' % tcdlme)
		allfaa += fetch(tcdlme, db='tcdb')
	fastas = {}
	for fa in allfaa.split('\n\n'):
		if not fa.strip(): continue
		for acc in accs:
			if acc in fastas: continue
			if fa.startswith('>' + acc):
				fastas[acc] = fa
	for x in sorted(fastas): 
		f = open(outdir + '/%s.fa' % x, 'w')
		f.write(fastas[x])
		f.close()

def quod_indiv(sequences, indir, outdir, dpi=300, force=False):

	if not os.path.isdir(outdir): os.mkdir(outdir)

	for x in sequences: 
		f = open(indir + '/%s.fa' % x)
		seq = f.read()
		f.close()
		title = x
		seq = seq[seq.find('\n')+1:].replace('\n', '')

		if not force and os.path.isfile(outdir + x + '.png'): continue

		quod.what([seq], title=title, imgfmt='png', directory=outdir, filename=(x+'.png'), dpi=dpi, hide=1)

def build_html(bc, indir, blasts, outdir='hvordan_out/html', filename='test.html'):

	if not os.path.isdir(outdir): os.mkdir(outdir)
	if not os.path.isdir(outdir + '/assets'): os.mkdir(outdir + '/assets')

	if not os.path.isfile(outdir + '/assets/openclose.js'):
		f = open(outdir + '/assets/openclose.js', 'w')
		f.write('function toggle_section(sectionid, selfid) {\n\tvar section = document.getElementById(sectionid);\n\tvar me = document.getElementById(selfid);\n\tconsole.log([section, section.style.display]);\n\tif (section.style.display == \'none\') {\n\t\tsection.style.display = \'block\';\n\t\tme.innerHTML = \'Hide\';\n\t} else { \n\t\tsection.style.display = \'none\'; \n\t\tme.innerHTML = \'Show\';\n\t}\n}')
		f.close()
	if not os.path.isfile(outdir + '/assets/nice.css'):
		f = open(outdir + '/assets/nice.css', 'w')
		f.write('body {\n\tfont-family: sans;\n\theight: 100%;\n}\ndiv {\n\tdisplay: block;\n}\ndiv.tcblast {\n\tmax-width: 1500px;\n}\ndiv.fullblast {\n\twidth: 50%;\n\tfloat: left;\n}\ndiv.tabular1 {\n\twidth: 50%;\n\tfloat: left;\n\theight: 100%;\n}\ndiv.tabular2 {\n\twidth: 50%;\n\tfloat: right;\n\theight: 100%;\n}\nimg.bluebarplot {\n\tmax-width: 100%;\n\theight: auto;\n}\n.clear { clear: both; }\n.scrollable {\n\toverflow-y: scroll;\n}\n.resizeable {\n\tresize: vertical;\n\toverflow: auto;\n\tborder-bottom: 1px solid gray;\n\tdisplay: block;\n}\n.bluebars {\n\theight: 25vh;\n}\n.pairwise {\n\theight: 50vh;\n}\n.whatall {\n\theight: 50vh;\n}\n.whataln {\n\twidth: 100%;\n}\n#seqs {\n\tdisplay: none;\n}\n\n\n\n.summtbl {\n\tfont-family: monospace, courier;\n\tfont-size: 75%;\n}\n.oddrow {\n\tbackground-color: #d8d8d8;\n}\ntd {\n\tpadding-right: 1em;\n}\n.red {\n\tcolor: red;\n}\nimg {\n\tborder: 1pt solid black;\n}\n')
		f.close()
	#bc := [WP_1234567890, AP_1234567890]
	title = 'HVORDAN summary: %s vs %s' % tuple(bc[1:3])

	out = '<html><head><title>%s</title>' % title
	out += '\n<link rel="stylesheet" type="text/css" href="assets/nice.css"/>'
	out += '\n<script src="assets/openclose.js"></script>'
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
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[1]
	out += '\n</div><div class="scrollable tabular2">'
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[2]
	out += '\n</div></div>'

	#pairwise
	out += '\n<div class="clear"></div><a name="pairwise"><h3>Pairwise</h3></a><div class="resizeable pairwise"><div class="scrollable tabular1">'
	warn(blasts[0][1])
	out += '\n%s' % blasts[0][1]
	out += '</div><div class="scrollable tabular2">'
	out += '\n%s' % blasts[1][1]
	out += '\n</div></div></div>'


	#abcd bc
	out += '\n<div class="clear"></div><a name="abcd"><h3>ABCD Hydropathy plots</h3></a>'
	out += '\n<button class="showhide" id="abcdsh" onclick="toggle_section(\'abcd\', \'abcdsh\')">Hide</button>'

	out += '\n<div class="whatall" id="abcd">'
	out += '\n<div class="tabular1">'
	out += '\nA<br/><img class="bluebarplot" id="plota" src="../graphs/%s.png"/><br/>' % (bc[0])
	out += '\nB<br/><img class="bluebarplot" id="plotb" src="../graphs/%s.png"/><br/>' % (bc[1])
	out += '\n</div><div class="tabular2">'
	out += '\nC<br/><img class="bluebarplot" id="plotc" src="../graphs/%s.png"/><br/>' % (bc[2])
	out += '\nD<br/><img class="bluebarplot" id="plotd" src="../graphs/%s.png"/><br/>' % (bc[3])
	out += '\n</div></div>'

	out += '\n<div class="clear"></div><br/><a name="bc"><h3>BC hydropathy plot</h3></a>'
	out += '\n<button class="showhide" id="bcsh" onclick="toggle_section(\'bc\', \'bcsh\')">Hide</button>'
	out += '\n<div class="resizeable whataln" id="bc"><div class="scrollable">'
	out += '<img class="bluebarplot" id="plotbc" src="../graphs/%s_vs_%s.png"/><br/>' % (bc[1], bc[2])
	out += '\n</div></div>'

	out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'

	out += '\n<br/><div style="height: 10ex"></div>'
	out += '\n</body></html>'

	f = open(outdir + '/' + filename, 'w')
	f.write(out)
	f.close()

def get_fulltrans(fams, bcs, abcd):

	pairs = zip(bcs[fams[0]], bcs[fams[1]])
	origs = [abcd[fams[0]], abcd[fams[1]]]

	fulltrans = []
	for p in pairs:
		fulltrans.append([origs[0][p[0]][1], p[0], p[1], origs[1][p[1]][1]])
	return fulltrans

def blastem(acc, indir, outdir, dpi=300):
	f = open(indir + '/sequences/' + acc + '.fa')
	seq= f.read()
	f.close()

	return tcblast.til_warum(seq, outfile='%s/graphs/TCBLAST_%s.png' % (outdir, acc), title=acc, dpi=dpi, outdir='%s/blasts' % outdir)

	#fn = outdir + '/' + filename + '.png'

	#blasts = tcblast.til_warum(seq, fn, dpi=dpi)
	#blasts = [tcblast.til_warum(l[0], args.o + '/images/' + accs[0] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop'), tcblast.til_warum(l[1], args.o + '/images/' + accs[1] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop')]

def summarize(p1d, p2d, outdir, minz=15, maxz=None, dpi=100, force=False):


	if not os.path.isdir(outdir): os.mkdir(outdir)

	if VERBOSITY: info('Reading Protocol2 report')
	f = open(p2d + '/report.tbl')
	p2report = f.read()
	f.close()

	fams, bcs, alnregs = parse_p2report(p2report, minz, maxz)

	if VERBOSITY: info('Selecting best A-B C-D pairs')
	abcd = seek_initial(p1d, bcs)

	allseqs = set()
	for fam in abcd:
		for bc in abcd[fam]:
			allseqs.add(bc)
			allseqs.add(abcd[fam][bc][1])

	#grab all relevant sequences and store them
	if VERBOSITY: info('Retrieving sequences')
	clean_fetch(allseqs, outdir + '/sequences', force=force)

	#make graphs for all individual full-lengthers
	if VERBOSITY: info('Generating QUOD plots')
	quod_indiv(allseqs, outdir + '/sequences', outdir + '/graphs', dpi=dpi, force=force)

	#make graphs for all pairs of sequences
	#%s_vs_%s.png % (B, C)
	for s1 in alnregs: 
		for s2 in alnregs[s1]: #print(s1, s2, alnregs[s1][s2])
			quod.what(alnregs[s1][s2], labels=[s1,s2], title='%s (red) vs %s (blue)' % (s1,s2), imgfmt='png', directory=outdir+'/graphs', filename='%s_vs_%s.png' % (s1,s2), dpi=dpi, hide=1)

	#def build_html(bc, indir, outdir='hvordan_out/html', filename='test.html'):
	#for s1 in alnregs:
	#	for s2 in alnregs[s1]: 
	fulltrans = get_fulltrans(fams, bcs, abcd)
	if VERBOSITY: info('Generating TCBLAST plots')
	blasts = {}
	for pair in fulltrans:
		blasts[tuple(pair)] = [blastem(pair[1], indir=outdir, outdir=outdir, dpi=dpi), blastem(pair[2], indir=outdir, outdir=outdir, dpi=dpi)]
	if VERBOSITY: info('Generating HTML')
	for pair in fulltrans:
		build_html(pair, indir=outdir, blasts=blasts[tuple(pair)], outdir=(outdir + '/html'), filename='%s_vs_%s.html' % tuple(pair[1:3]))
		
if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser(description='HTML Visualization of Reasonable, Decent Alignment Networks')

	parser.add_argument('--p1d', default='.', help='famXpander directory. Note: Running "cut -f1-6" on psiblast.tbl will greatly improve performance')
	parser.add_argument('--p2d', default='.', help='Protocol2 results directory')

	parser.add_argument('--outdir', default='hvordan_out', help='output directory')

	parser.add_argument('-z', '--z-min', default=15, type=int, help='minimum Z score {default:15}')
	parser.add_argument('--z-max', default=None, type=int, help='maximum Z score')

	parser.add_argument('-f', action='store_true', help='force redownloads where applicable')
	parser.add_argument('--dpi', type=int, default=300, help='resolution of graphs {default:300}')

	args = parser.parse_args()

	summarize(args.p1d, args.p2d, args.outdir, minz=args.z_min, maxz=args.z_max, dpi=args.dpi, force=args.f)
