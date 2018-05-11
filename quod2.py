#!/usr/bin/env python2
#QUOD2: Questionable Utility of Doom 2
#most code copied from gblast3.py (which was written by Vamsee Reddy and Vasu Pranav Sai Iddamsetty) except where noted
#-Kevin Hendargo
from __future__ import division, print_function, division
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import os, subprocess, re, sys
import argparse
import numpy as np

from Bio import SeqIO
def sanitize(s): return re.sub('/', '', s)

def roundto(x, n):
	return (x//5)*5

def overlap(span1, span2):
	if span1[0] <= span2[0] <= span1[-1]: return True
	elif span1[0] <= span2[1] <= span1[-1]: return True
	elif span2[0] <= span1[0] <= span2[-1]: return True
	elif span2[0] <= span1[1] <= span2[-1]: return True
	else: return False

def union(span1, span2):
	if not overlap(span1, span2): raise NotImplementedError
	else: return [min(span1[0], span2[0]), max(span1[-1], span2[-1])]

def hex2tuple(s):
	if s.startswith('#'): s = s[1:]
	elif s.startswith('0x'): s = s[2:]

	if len(s) == 3: s = [2*s[i] for i in range(len(s))]
	if len(s) == 1: s *= 6
	
	l = [int(s[i:i+2], 16)/255. for i in range(0, len(s), 2)]
	return tuple(l)

class Plot(object):
	def __init__(self): 
		self.fig = Figure()
		self.canvas = FigureCanvas(self.fig)

		self.ax = self.fig.add_subplot(111)
		self.width, self.height = None, None
		self.xlim = [0, 20]
		self.ylim = [-3, 3]
		self.xticks, self.yticks = None, 1

		self.axeslabels = ['X', 'Y']
		self.titles = []

	def render(self):
		self.ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
		self.ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])
		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()

		maxl = xlim[1] - xlim[0]
		if self.xticks is None:
			self.xticks = roundto(maxl // 5, 5)

		self.ax.set_xticks(np.arange(xlim[0], xlim[1]+1, self.xticks))
		self.ax.set_yticks(np.arange(ylim[0], ylim[1]+1, self.yticks))

		#TODO: allow centimeters or something
		#if self.width is None: self.width = (0.0265) * (maxl) if maxl > 600 else 15.
		if self.width is None: self.width = (0.0265) * (maxl) if maxl > 200 else 15.
		if self.height is None: self.height = 5.5
		self.fig.set_figwidth(self.width)
		self.fig.set_figheight(self.height)
		#self.fig.set_tight_layout(True)
		self.ax.set_xlabel(self.axeslabels[0])
		self.ax.set_ylabel(self.axeslabels[1])

		#only when plotting hydropathy! FIXME
		self.ax.axhline(y=0, color='black', linewidth=0.5)

		title = ''
		for t in self.titles: 
			if t is not None: title += '{}, '.format(t)
		self.ax.set_title(title[:-2])

	def save(self, filename, dpi=80, format='png'):
		self.fig.savefig(filename, dpi=dpi, format=format)

class Entity(object):
	def __init__(self): pass

	def draw(self, plot): raise NotImplementedError

class Curve(Entity):
	def __init__(self, X=[], Y=[], style='auto'):
		Entity.__init__(self)
		self.X = X
		self.Y = Y
		self.style = style

	def draw(self, plot): 
		if len(Y) and not len(X): 
			plot.ax.plot(Y, style)
		else: plot.ax.plot(X, Y, style)

class Vspans(Entity):
	def __init__(self, spans=[], style='orange', alpha=None):
		Entity.__init__(self)
		self.spans = spans
		self.style = style
		self.alpha = alpha

	def get_alpha(self):
		if type(self.style) is tuple and len(self.style) == 4: return self.style[3]
		elif self.alpha is not None: return self.alpha
		else: return 0.25

	def draw(self, plot):
		for span in self.spans:
			plot.ax.axvspan(span[0], span[1], facecolor=self.style, alpha=self.get_alpha())

class Wall(Vspans):
	def __init__(self, spans=[], y=None, ylim=[0,1], style='black', wedge=1, single=False):
		Vspans.__init__(self, spans=spans, style=style)
		self.y = y
		self.ylim = ylim
		self.wedge = wedge
		self.single = single

	def get_y(self):
		if type(self.y) is None: return 2
		else: return self.y

	def get_ylim(self):
		if self.ylim is None: return [0, 1]
		elif type(self.ylim) is int:
			if self.ylim > 0: return [0.5, 1]
			elif self.ylim == 0: return [0, 1]
			else: return [0, 0.5]
		else: return self.ylim

	def draw(self, plot):
		ylim = self.get_ylim()
		for span in self.spans:
			for i, x in enumerate(span): 
				plot.ax.axvline(x=x, color='black', ymin=ylim[0], ymax=ylim[1])
				if self.wedge:
					if single:
						if wedge >= 0: marker = '>'
						else: marker = '<'
					else:
						if i % 2: marker = '<'
						else: marker = '>'

					#wx = x# + (abs(4*self.wedge)**0.5 * np.sign(0.5 - i % 2))
					#FIXME: make the arrows aligned at any scale
					wx = x - 2*(abs(self.wedge) * np.sign((i % 2) - 0.5)) - self.wedge*.5

					if self.y >= 0: wy = 2
					else: wy = -2

					plot.ax.scatter([wx], [wy], marker=marker, color='black', s=25*abs(self.wedge))

class HMMTOP(Vspans):
	def __init__(self, gseq, style='orange'):
		Vspans.__init__(self, style=style)

		fasta = gseq
		if not fasta.startswith('>'): fasta = '>seq\n' + fasta
		p = subprocess.Popen(['hmmtop', '-if=--', '-is=pseudo', '-sf=FAS', '-pi=spred'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate(input=fasta)
		print(err.strip(), file=sys.stderr)

		indices = re.findall('(?:IN|OUT)((?:\s*(?:[0-9]+))+)', out.strip())[0].strip().split()
		indices = [int(i) for i in indices[1:]]

		if not indices: return

		if gseq.startswith('>'): seq = gseq[gseq.find('\n')+1:]
		else: seq = gseq

		seq = re.sub('[^A-Z\-]', '', seq)

		for i, c in enumerate(seq):
			if c in 'ACDEFGHIKLMNPQRSTVWY': pass
			else: 
				for j, n in enumerate(indices):
					if (n - 1) >= i: indices[j] += 1
		self.spans = [[indices[i], indices[i+1]] for i in range(0, len(indices), 2)]

	def add(self, *spans):
		for span in spans: self.spans.append(span)

	def replace(self, *spans):
		removeme = []
		for newspan in spans:
			for j, oldspan in enumerate(self.spans):
				if overlap(newspan, oldspan): 
					removeme.append(j)
		for j in removeme[::-1]: self.spans.pop(j)
		for newspan in spans: self.spans.append(newspan)

	def extend(self, *spans):
		fuseme = []
		appendme = []
		for i, newspan in enumerate(spans):
			for j, oldspan in enumerate(self.spans):
				if overlap(newspan, oldspan): fuseme.append([i, j])
				else: appendme.append(i)
		for pair in fuseme:
			self.spans[pair[1]] = union(self.spans[pair[1]], spans[pair[0]])
		for span in appendme: self.spans.append(span)

class Hydropathy(Curve):

	def __init__(self, gseq, style='r', offset=0, window=19):
		Curve.__init__(self, style=style)
		self.window = window
		#Copied with barely any modification from gblast3.py

 		#sanitize the sequence, which may be a FASTA
		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		seq = re.sub('[X\-]', '', gseq)
		seq = re.sub('[^A-Z]', '', seq)
		prev = 0
		index = {'G':(-0.400,0.48), \
             'I':(4.500,1.38), \
             'S':(-0.800,-0.18), \
			 'Q':(-3.500,-0.85), \
             'E':(-3.500,-0.74), \
             'A':(1.800,0.62), \
             'M':(1.900,0.64), \
             'T':(-0.700,-0.05), \
             'Y':(-1.300,0.26), \
             'H':(-3.200,-0.4), \
             'V':(4.200,1.08), \
             'F':(2.800,1.19), \
             'C':(2.500,0.29), \
             'W':(-0.900,0.81), \
             'K':(-3.900,-1.5), \
             'L':(3.800,1.06), \
             'P':(-1.600,0.12), \
             'N':(-3.500,-0.78), \
             'D':(-3.500,-0.90), \
             'R':(-4.500,-2.53), \
             'U':(0,0), \
             'B':(-3.500,-0.84), \
             'J':(-3.500,-0.80), \
             'Z':(4.150,1.22) \
            }

		midpt = (window+1)//2
		length = len(seq)
		hydro = []
		for i in range(length-window+1):
			total = 0
			for j in range(window): total += index[seq[i+j]][0]
			total /= window
			hydro.append(total)

		if len(seq) == len(gseq): 
			self.Y = np.array(hydro)
			self.X = np.arange(0, len(self.Y))+window//2+offset+1

		else:
			replace = re.finditer('(-+|X+)', gseq)
			inserts = {}
			for i in replace: inserts.setdefault(i.start(), i.end()-i.start())
			first = False
			newhydro = []

			for x, h in enumerate(hydro):
				if x in inserts and not first:
					first = True
					newhydro += [np.nan for y in range(inserts[x])]
					newcount = x + inserts[x]

				elif not first: newhydro.append(h)

				elif first and newcount in inserts:
					newhydro += [np.nan for y in range(inserts[newcount])]
					newcount += inserts[newcount]
				else:
					newhydro.append(h)
					newcount += 1

			self.Y = np.array(newhydro)
			self.X = np.arange(offset+1, len(self.Y)+1)+window//2
	def draw(self, plot):
		plot.xlim[0] = min(plot.xlim[0], self.X[0] - self.window//2)
		plot.xlim[1] = max(plot.xlim[1], self.X[-1] + self.window//2)

		plot.axeslabels = ['Residue number', 'Hydropathy (kcal/mol)']
		plot.ax.plot(self.X, self.Y, color=self.style, linewidth=1)

class What(Entity):
	def __init__(self, seq, style=0, tmscolor=None, linecolor=None):
		Entity.__init__(self)
		self.seq = seq
		self.style = style
		self.tmscolor = tmscolor
		self.linecolor = linecolor
		self.entities = []
		self.entities.append(Hydropathy(seq, style=self.get_curve_color()))
		self.entities.append(HMMTOP(seq, style=self.get_tms_color()))

	def get_title(self, showlength=True):
		if self.seq.startswith('>'): 
			s = self.seq[1:self.seq.find('\n')]
			if showlength: 
				truelength = len(re.sub('[^A-Z]', '', self.seq[self.seq.find('\n')+1:]))
				s += ' ({}aa)'.format(truelength)
			return s
		else: 
			PREVIEW = [8, 3]
			if len(self.seq) <= sum(PREVIEW): return '{} ({}aa)'.format(self.seq, len(self.seq))
			else: return '{}...{} ({}aa)'.format(self.seq[:8], self.seq[-3:], len(self.seq))

	def get_curve_color(self):
		if self.linecolor is None:
			if (self.style % 3) == 0: return 'red'
			elif (self.style % 3) == 1: return 'blue'
			elif (self.style % 3) == 2: return 'green'
		elif self.linecolor.startswith('#') or self.linecolor.startswith('0x'):
			return hex2tuple(self.linecolor)
		else: return self.linecolor 

	def get_tms_color(self):
		if self.tmscolor is None:
			if (self.style % 3) == 0: return 'orange'
			elif (self.style % 3) == 1: return 'cyan'
			elif (self.style % 3) == 2: return 'purple'
		elif self.tmscolor.startswith('#') or self.tmscolor.startswith('0x'):
			return hex2tuple(self.tmscolor)
		else: return self.tmscolor

	def draw(self, plot):
		#for e in self.entities: e.draw(notitle=True)
		for e in self.entities: e.draw(plot)

		#XXX maybe this goes better under Hydropathy or something
		plot.titles.append(self.get_title())

def split_fasta(fastastr):
	sequences = []
	for l in fastastr.split('\n'):
		if l.startswith('>'): sequences.append(l)
		elif not l.strip(): continue
		else: 
			try: sequences[-1] += '\n' + l.strip()
			except IndexError: sequences.append('>untitled sequence\n{}'.format(l))
	return sequences

#def main(infiles, mode='hydropathy', walls=None, ):
#def what(sequences, labels=None, imgfmt='png', directory=None, filename=None, title=False, dpi=80, hide=True, viewer=None, bars=[], color='auto', offset=0, statistics=False, overwrite=False, manual_tms=None, wedges=[], ywedge=2, legend=False, window=19, silent=False, axisfont=None, tickfont=None, xticks=None, mode='hydropathy', width=None, height=None):
def main(infiles, **kwargs):
	plot = Plot()
	entities = []

	n = 0
	if 'force_seq' in kwargs and kwargs['force_seq']: 
		for seq in infiles: 
			entities.append(What(seq, style=n))
			n += 1
	else:
		for fn in infiles:
			with open(fn) as f:
				for seq in split_fasta(f.read().decode('utf-8')):
					entities.append(What(seq, style=n))
					n += 1

	if 'bars' in kwargs and kwargs['bars'] is not None: 
		[entities.append(wall) for wall in parse_walls(kwargs['bars'], wedge=0)]

	if 'dpi' in kwargs: dpi = kwargs['dpi']
	else: dpi = 80

	if 'imgfmt' in kwargs: imgfmt = kwargs['imgfmt']
	else: imgfmt = 'png'

	if 'walls' in kwargs and kwargs['walls'] is not None: 
		[entities.append(wall) for wall in parse_walls(kwargs['walls'])]

	if 'wall' in kwargs and kwargs['wall'] is not None:
		[entities.append(wall) for wall in parse_walls(kwargs['wall'], single=1)]

	if 'outdir' in kwargs and kwargs['outdir']: prefix = kwargs['outdir']
	else: prefix = ''

	if 'outfile' in kwargs and kwargs['outfile']: 
		if prefix: outfile = '{}/{}'.format(prefix, kwargs['outfile'])
		else: outfile = '{}'.format(kwargs['outfile'])
	else: outfile = None

	for e in entities: e.draw(plot)

	plot.render()
	if outfile is None: 
		outfile = sanitize(plot.ax.get_title())

	if len(os.path.splitext(outfile)) == 1: outfile += '.{}'.format(imgfmt)
	elif len(os.path.splitext(outfile)[-1]) not in [3, 4]: outfile += '.{}'.format(imgfmt)

	plot.save(outfile, dpi=dpi, format=imgfmt)
	os.system('xdg-open "{}"'.format(outfile))
		

def parse_walls(strlist, wedge=1, single=False):
	#turns a list of wall specifications into Wall() objects
	if not strlist: return None
	out = []
	for wall in strlist:
		if type(wall) is str:
			coords = [int(x) for x in wall.split(',')]
			if len(coords) == 1: coords.append(5)
			if len(coords) == 2: coords.append(None)

			span = [coords[0], coords[0]+coords[1]]
			out.append(Wall([span], y=coords[2], ylim=coords[2], wedge=wedge, single=single))
		elif type(wall) is int:
			out.append(Wall([[wall, wall]], wedge=wedge, single=single))
	return out

if __name__ == '__main__': 
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='sequence files to read in')
	parser.add_argument('-b', '--bars', nargs='+', type=int, help='Draws vertical bars at these positions')
	parser.add_argument('-d', metavar='outdir', help='Directory to store graphs in (recommended only with autogenerated filenames)')
	parser.add_argument('-o', metavar='outfile', help='Filename of graph, relative to the argument of -d if present and as normally interpreted otherwise')
	parser.add_argument('-r', metavar='resolution', type=int, help='Resolution of graph in dpi. The default is 80dpi, which is suitable for viewing on monitors. Draft-quality images should be about 300dpi, and publication-quality images need to be 600 or 1200dpi depending on the journal.')
	parser.add_argument('-s', action='store_true', help='Force inputs to be interpreted as sequences (this is no longer a default behavior for infile args matching /[A-Z]+/')
	parser.add_argument('-t', metavar='format', default='png', help='Format of graph (\033[1mpng\033[0m, eps, jpeg, jpg, pdf, pgf, ps, raw, rgba, svg, svgz, tif, tiff')
	parser.add_argument('-W', '--wall', metavar='x(,dx(,y))', nargs='+', help='Draws bounds around sequences and such')
	parser.add_argument('-w', '--walls', metavar='x(,dx(,y))', nargs='+', help='Draws bounds around sequences and such')
	parser.add_argument('--mode', default='hydropathy', help='mode to run QUOD in (\033[1mhydropathy\033[0m, entropy)')
	parser.add_argument('--viewer', metavar='viewer', default=None, help='Viewer to be used for opening plots')
	parser.add_argument('--window', metavar='windowsize', default=19, help='Window size for hydropathy')

	args = parser.parse_args()

	main(args.infile, mode=args.mode, walls=args.walls, wall=args.wall, bars=args.bars, dpi=args.r, imgfmt=args.t, force_seq=args.s, outdir=args.d, outfile=args.o)

#with open('gapxmfs_2.A.1.1.1.fa') as f: gapseq = f.read().decode('utf-8')
#with open('mfs_2.A.1.1.1.fa') as f: seq = f.read().decode('utf-8')
#plot = Plot()
#entities = []
#entities.append(What(plot, seq, style=0))
#entities.append(What(plot, gapseq, style=1))
##entities.append(HMMTOP(plot, seq))
#for e in entities: e.draw()
#plot.save('test.png')
#os.system('xdg-open test.png')

