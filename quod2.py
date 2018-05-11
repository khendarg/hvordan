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

def roundto(x, n):
	return (x//5)*5

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
		self.render()

		self.fig.savefig(filename, dpi=dpi, format=format)

class Entity(object):
	def __init__(self, plot): 
		self.plot = plot

	def draw(self): raise NotImplementedError

class Curve(Entity):
	def __init__(self, plot, X=[], Y=[], style='auto'):
		Entity.__init__(self, plot)
		self.X = X
		self.Y = Y
		self.style = style

	def draw(self): 
		if len(Y) and not len(X): 
			self.plot.ax.plot(Y, style)
		else: self.plot.ax.plot(X, Y, style)

class Vspans(Entity):
	def __init__(self, plot, spans=[], style='orange', alpha=None):
		Entity.__init__(self, plot)
		self.spans = spans
		self.style = style
		self.alpha = alpha

	def get_alpha(self):
		if type(self.style) is tuple and len(self.style) == 4: return self.style[3]
		elif self.alpha is not None: return self.alpha
		else: return 0.25

	def draw(self):
		for span in self.spans:
			self.plot.ax.axvspan(span[0], span[1], facecolor=self.style, alpha=self.get_alpha())

class HMMTOP(Vspans):
	def __init__(self, plot, gseq, style='orange'):
		Vspans.__init__(self, plot, style=style)

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

class Hydropathy(Curve):

	def __init__(self, plot, gseq, style='r', offset=0, window=19):
		Curve.__init__(self, plot, style=style)
		self.window = window
		#Copied with barely any modification from gblast3.py

 		#sanitize the sequence, which may be a FASTA
		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		seq = re.sub('[X\-]', '', gseq)
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
	def draw(self):
		self.plot.xlim[0] = min(self.plot.xlim[0], self.X[0] - self.window//2)
		self.plot.xlim[1] = max(self.plot.xlim[1], self.X[-1] + self.window//2)

		self.plot.axeslabels = ['Residue number', 'Hydropathy (kcal/mol)']
		self.plot.ax.plot(self.X, self.Y, color=self.style, linewidth=1)

class What(Entity):
	def __init__(self, plot, seq, style=0, tmscolor=None, linecolor=None):
		Entity.__init__(self, plot)
		self.seq = seq
		self.style = style
		self.tmscolor = tmscolor
		self.linecolor = linecolor
		self.entities = []
		self.entities.append(Hydropathy(self.plot, seq, style=self.get_curve_color()))
		self.entities.append(HMMTOP(self.plot, seq, style=self.get_tms_color()))

	def get_title(self, showlength=True):
		if self.seq.startswith('>'): 
			s = self.seq[1:self.seq.find('\n')]
			if showlength: 
				truelength = len(re.sub('[^A-Z]', '', self.seq[self.seq.find('\n')+1:]))
				s += ' ({}aa)'.format(truelength)
			return s
		else: return '{}...'.format(self.seq[:8])

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

	def draw(self):
		#for e in self.entities: e.draw(notitle=True)
		for e in self.entities: e.draw()

		#XXX maybe this goes better under Hydropathy or something
		self.plot.titles.append(self.get_title())

with open('gapxmfs_2.A.1.1.1.fa') as f: gapseq = f.read().decode('utf-8')
with open('mfs_2.A.1.1.1.fa') as f: seq = f.read().decode('utf-8')
plot = Plot()
entities = []
entities.append(What(plot, seq, style=0))
entities.append(What(plot, gapseq, style=1))
#entities.append(HMMTOP(plot, seq))
for e in entities: e.draw()
plot.save('test.png')
os.system('xdg-open test.png')
#fig = Figure()
## A canvas must be manually attached to the figure (pyplot would automatically
## do it).  This is done by instantiating the canvas with the figure as
## argument.
#FigureCanvas(fig)
#ax = fig.add_subplot(111)
#ax.plot([1, 2, 3])
#ax.set_title('hi mom')
#ax.grid(True)
#ax.set_xlabel('time')
#ax.set_ylabel('volts')
##fig.show()
#
#subprocess.call(['xdg-open', 'test.png'])
#
