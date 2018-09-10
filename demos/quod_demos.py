#!/usr/bin/env python2

from __future__ import print_function, division

import quod

SEQ1 = '''>tr|Q9JHI9 SLC11A3 iron transporter (Ferroportin1) (Solute carrier family 39 (iron-regulated transporter), member 1) - Mus musculus (Mouse).
MTKARDQTHQEGCCGSLANYLTSAKFLLYLGHSLSTWGDRMWHFAVSVFLVELYGNSLLLTAVYGLVVAGSVLVLGAIIG
DWVDKNARLKVAQTSLVVQNVSVILCGIILMMVFLHKNELLTMYHGWVLTVCYILIITIANIANLASTATAITIQRDWIV
VVAGENRSRLADMNATIRRIDQLTNILAPMAVGQIMTFGSPVIGCGFISGWNLVSMCVEYFLLWKVYQKTPALAVKAALK
VEESELKQLTSPKDTEPKPLEGTHLMGEKDSNIRELECEQEPTCASQMAEPFRTFRDGWVSYYNQPVFLAGMGLAFLYMT
VLGFDCITTGYAYTQGLSGSILSILMGASAITGIMGTVAFTWLRRKCGLVRTGLFSGLAQLSCLILCVISVFMPGSPLDL
SVSPFEDIRSRFVNVEPVSPTTKIPETVFTTEMHMSNMSNVHEMSTKPIPIVSVSLLFAGVIAARIGLWSFDLTVTQLLQ
ENVIESERGIINGVQNSMNYLLDLLHFIMVILAPNPEAFGLLVLISVSFVAMGHLMYFRFAQKTLGNQIFVCGPDEKEVT
DENQPNTSVV'''

SEQ2 = '''>sp|P26678|PPLA_HUMAN CARDIAC PHOSPHOLAMBAN (PLB) - Homo sapiens (Human).
MEKVQYLTRSAIRRASTIEMPQQARQKLQNLFINFCLILICLLLICIIVMLL'''

SEQ3 = '''>sp|Q92508|PIEZ1_HUMAN Piezo-type mechanosensitive ion channel component 1 OS=Homo sapiens GN=PIEZO1 PE=1 SV=4
MEPHVLGAVLYWLLLPCALLAACLLRFSGLSLVYLLFLLLLPWFPGPTRCGLQGHTGRLLRALLGLSLLFLVAHLALQIC
LHIVPRLDQLLGPSCSRWETLSRHIGVTRLDLKDIPNAIRLVAPDLGILVVSSVCLGICGRLARNTRQSPHPRELDDDER
DVDASPTAGLQEAATLAPTRRSRLAARFRVTAHWLLVAAGRVLAVTLLALAGIAHPSALSSVYLLLFLALCTWWACHFPI
STRGFSRLCVAVGCFGAGHLICLYCYQMPLAQALLPPAGIWARVLGLKDFVGPTNCSSPHALVLNTGLDWPVYASPGVLL
LLCYATASLRKLRAYRPSGQRKEAAKGYEARELELAELDQWPQERESDQHVVPTAPDTEADNCIVHELTGQSSVLRRPVR
PKRAEPREASPLHSLGHLIMDQSYVCALIAMMVWSITYHSWLTFVLLLWACLIWTVRSRHQLAMLCSPCILLYGMTLCCL
RYVWAMDLRPELPTTLGPVSLRQLGLEHTRYPCLDLGAMLLYTLTFWLLLRQFVKEKLLKWAESPAALTEVTVADTEPTR
TQTLLQSLGELVKGVYAKYWIYVCAGMFIVVSFAGRLVVYKIVYMFLFLLCLTLFQVYYSLWRKLLKAFWWLVVAYTMLV
LIAVYTFQFQDFPAYWRNLTGFTDEQLGDLGLEQFSVSELFSSILVPGFFLLACILQLHYFHRPFMQLTDMEHVSLPGTR
LPRWAHRQDAVSGTPLLREEQQEHQQQQQEEEEEEEDSRDEGLGVATPHQATQVPEGAAKWGLVAERLLELAAGFSDVLS
RVQVFLRRLLELHVFKLVALYTVWVALKEVSVMNLLLVVLWAFALPYPRFRPMASCLSTVWTCVIIVCKMLYQLKVVNPQ
EYSSNCTEPFPNSTNLLPTEISQSLLYRGPVDPANWFGVRKGFPNLGYIQNHLQVLLLLVFEAIVYRRQEHYRRQHQLAP
LPAQAVFASGTRQQLDQDLLGCLKYFINFFFYKFGLEICFLMAVNVIGQRMNFLVTLHGCWLVAILTRRHRQAIARLWPN
YCLFLALFLLYQYLLCLGMPPALCIDYPWRWSRAVPMNSALIKWLYLPDFFRAPNSTNLISDFLLLLCASQQWQVFSAER
TEEWQRMAGVNTDRLEPLRGEPNPVPNFIHCRSYLDMLKVAVFRYLFWLVLVVVFVTGATRISIFGLGYLLACFYLLLFG
TALLQRDTRARLVLWDCLILYNVTVIISKNMLSLLACVFVEQMQTGFCWVIQLFSLVCTVKGYYDPKEMMDRDQDCLLPV
EEAGIIWDSVCFFFLLLQRRVFLSHYYLHVRADLQATALLASRGFALYNAANLKSIDFHRRIEEKSLAQLKRQMERIRAK
QEKHRQGRVDRSRPQDTLGPKDPGLEPGPDSPGGSSPPRRQWWRPWLDHATVIHSGDYFLFESDSEEEEEAVPEDPRPSA
QSAFQLAYQAWVTNAQAVLRRRQQEQEQARQEQAGQLPTGGGPSQEVEPAEGPEEAAAGRSHVVQRVLSTAQFLWMLGQA
LVDELTRWLQEFTRHHGTMSDVLRAERYLLTQELLQGGEVHRGVLDQLYTSQAEATLPGPTEAPNAPSTVSSGLGAEEPL
SSMTDDMGSPLSTGYHTRSGSEEAVTDPGEREAGASLYQGLMRTASELLLDRRLRIPELEEAELFAEGQGRALRLLRAVY
QCVAAHSELLCYFIIILNHMVTASAGSLVLPVLVFLWAMLSIPRPSKRFWMTAIVFTEIAVVVKYLFQFGFFPWNSHVVL
RRYENKPYFPPRILGLEKTDGYIKYDLVQLMALFFHRSQLLCYGLWDHEEDSPSKEHDKSGEEEQGAEEGPGVPAATTED
HIQVEARVGPTDGTPEPQVELRPRDTRRISLRFRRRKKEGPARKGAAAIEAEDREEEEGEEEKEAPTGREKRPSRSGGRV
RAAGRRLQGFCLSLAQGTYRPLRRFFHDILHTKYRAATDVYALMFLADVVDFIIIIFGFWAFGKHSAATDITSSLSDDQV
PEAFLVMLLIQFSTMVVDRALYLRKTVLGKLAFQVALVLAIHLWMFFILPAVTERMFNQNVVAQLWYFVKCIYFALSAYQ
IRCGYPTRILGNFLTKKYNHLNLFLFQGFRLVPFLVELRAVMDWVWTDTTLSLSSWMCVEDIYANIFIIKCSRETEKKYP
QPKGQKKKKIVKYGMGGLIILFLIAIIWFPLLFMSLVRSVVGVVNQPIDVTVTLKLGGYEPLFTMSAQQPSIIPFTAQAY
EELSRQFDPQPLAMQFISQYSPEDIVTAQIEGSSGALWRISPPSRAQMKRELYNGTADITLRFTWNFQRDLAKGGTVEYA
NEKHMLALAPNSTARRQLASLLEGTSDQSVVIPNLFPKYIRAPNGPEANPVKQLQPNEEADYLGVRIQLRREQGAGATGF
LEWWVIELQECRTDCNLLPMVIFSDKVSPPSLGFLAGYGIMGLYVSIVLVIGKFVRGFFSEISHSIMFEELPCVDRILKL
CQDIFLVRETRELELEEELYAKLIFLYRSPETMIKWTREKE'''

ALN1 = '''>P0C0Y8
MALLS--FER-KYRVPG--GTLVGGNLFDF-WVGPF--YVGFFGVATF--
-FFAALGIILI----AWSAVL------QGTWNPQLI-------SVYPPAL
EYGLG-GAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFA
FAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYN
PAHMIAISFFFTNALALALHGALVLSAANPEKGKEMR-------TPDHED
TFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWW
QWWVKLPWWANIPGGING'''
ALN2 = '''>P0C0Y9
MAEYQNIFSQVQVRGPADLGMTEDVNLANRSGVGPFSTLLGWFGNAQLGP
IYLGSLGVLSLFSGLMWFFTIGIWFWYQAGWNPAVFLRDLFFFSLEPPAP
EYGLSFAAPLKEGGLWLIASFFMFVAVWSWWGRTYLRAQALGMGKHTAWA
FLSAIWLWMVLGFIRPILMGSWSEAVPYGIFSHLDWTNNFSLVHGNLFYN
PFHGLSIAFLYGSALLFAMHGATILAVSRFGGERELEQIADRGTAAERAA
LFWRWTMGFNATMEGIHRWAIWMAVLVTLTGGIGILLSGTV-VDNWYVWG
QNHGMAP--------LN-'''

def simple_quod(*sequences, **kwargs):
	'''This example creates a standard hydropathy plot from several raw sequences and saves it with a specific filename if given

	*sequences: raw sequences
	**kwargs:
		window: window size for hydropathy averaging (default: 19)
		outfile: what to name the resulting figure'''

	outfile = 'test.png' if ('outfile' not in kwargs) else kwargs['outfile']
	window = 19 if ('window' not in kwargs) else kwargs['window']

	plot = quod.Plot()

	for i, seq in enumerate(sequences):
		plot.add(quod.What(seq, style=i, window=window))
	
	plot.render()

	plot.save(outfile)


#generates a standard QUOD plot
simple_quod(SEQ1, SEQ2, outfile='simple_quod.png')


#generates a QUOD plot with a window size of 9
simple_quod(SEQ1, SEQ2, window=9, outfile='noisy_quod.png')


def colored_quod(seq1, seq2, linecolor1='red', linecolor2='blue', tmscolor1='orange', tmscolor2='cyan', outfile='test.png'):
	'''This example creates a standard hydropathy plot for two raw sequences, making sure to specify colors manually

	seq1: sequence 1
	seq2: sequence 2
	linecolor1: line color for sequence 1
	linecolor2: line color for sequence 2
	tmscolor1: tms color for sequence 1 
	tmscolor2: tms color for sequence 2
	outfile: what to name the resulting figure'''

	plot = quod.Plot()

	plot.add(quod.What(seq1, tmscolor=tmscolor1, linecolor=linecolor1))
	plot.add(quod.What(seq2, tmscolor=tmscolor2, linecolor=linecolor2))

	plot.render()

	plot.save(outfile)

colored_quod(ALN1, ALN2, linecolor1='red', linecolor2='blue', tmscolor1='orange', tmscolor2='cyan', outfile='colored_quod.png')
colored_quod(ALN1, ALN2, linecolor1='red', linecolor2='blue', tmscolor1='red', tmscolor2='blue', outfile='colored_quod_gblast.png')
colored_quod(ALN1, ALN2, linecolor1='orange', linecolor2='cyan', tmscolor1='red', tmscolor2='blue', outfile='colored_quod_bizarre.png')


class Interval(object):
	'''This class provides an overview of the most commonly used QUOD plot features'''
	def __init__(self, left, right, y=None, scale=1, label=None, style=None, mode=None):
		'''left: left bound of interval
		right: right bound of interval
		y: y-level to place resulting object
		label: if applicable, what to label the resulting object
		style: color specifier (ints result in automatic rotating color assignment and strs result in explicit color assignment)
		mode: used for automatic feature selection'''
		self.left = left
		self.right = right
		self.y = y
		self.label = label
		self.style = style
		self.mode = mode
	def __iter__(self): return iter([self.left, self.right])

	def get_wall(self): 
		'''return an object corresponding to the bounds of an interval with arrows pointing in'''
		return quod.Wall([[self.left, self.right]], y=self.y, ylim=self.y, style=self.style)

	def get_tms(self): 
		'''return an object realized as transparent boxes'''
		h = quod.HMMTOP(gseq='', style=self.style, nohmmtop=True)
		h.spans = [[self.left, self.right]]
		return h

	def get_domain(self): 
		'''return an object consisting of some text and a strongly colored box'''
		y = -2 if self.y is None else self.y
		return quod.Region(spans=[[self.left, self.right]], yspan=[y-0.15, 0.15], label=self.label, style=self.style)

	def get_wedges(self): 
		'''return a list of Wedge objects to illustrate a deprecated QUOD feature'''
		return [quod.Wedge(x, y=self.y, scale=scale, style=self.style) for x in self]

	def get_auto(self):
		if self.mode == 'wall': return self.get_wall()
		elif self.mode == 'tms': return self.get_tms()
		elif self.mode == 'domain': return self.get_domain()
		elif self.mode == 'wedges': return self.get_wedges()
		else: return self.get_wall()

def walls_quod(*sequences, **kwargs):
	'''This example plots sequences with walls+wedges defining certain intervals

	*sequences: raw sequences
	**kwargs:
		intervals: a list of Interval objects to generate walls+wedges from (default: [])
		outfile: what to name the resulting figure (default: test.png)'''

	plot = quod.Plot()

	intervals = [] if 'intervals' not in kwargs else kwargs['intervals']
	outfile = 'test.png' if 'outfile' not in kwargs else kwargs['outfile']
	mode = 'wall' if 'mode' not in kwargs else kwargs['mode']

	for i in intervals: plot.add(i.get_auto())

	for n, seq in enumerate(sequences): plot.add(quod.What(seq, style=n))

	plot.render()

	plot.save(outfile)


#generates a plot with a basic wall covering 20-40
intervals = [Interval(20, 40, mode='wall')]
walls_quod(ALN1, ALN2, intervals=intervals, outfile='quod_basic_wall.png')


#generates a plot with an extra green TMS covering 20-40
intervals = [Interval(20, 40, mode='tms', style='green')]
walls_quod(ALN1, ALN2, intervals=intervals, outfile='quod_basic_tms.png')


#generates a plot with a purple domain/region marker covering 20-40
intervals = [Interval(20, 40, mode='domain', style='purple', label='Some domain')]
walls_quod(ALN1, ALN2, intervals=intervals, outfile='quod_basic_domain.png')


#generates a plot with all of the above, more or less
intervals = [Interval(20, 60, mode='wall')]
intervals.append(Interval(70, 110, mode='tms', style='red'))
intervals.append(Interval(120, 160, mode='domain', style='purple', label='Some domain'))
walls_quod(ALN1, ALN2, intervals=intervals, outfile='quod_all_above.png')


#generates a plot with above-y walls and under-y walls 
intervals = [Interval(20, 40, mode='wall', y=-2)]
intervals.append(Interval(120, 140, mode='wall', y=2))
walls_quod(ALN1, ALN2, intervals=intervals, outfile='quod_splitwall.png')


def entropy_quod(*sequences, **kwargs):
	'''This example plots the entropies of multiple sequences

	sequences: raw sequences
	**kwargs:
		outfile: where to put the resulting figure'''

	outfile = kwargs['outfile'] if 'outfile' in kwargs else kwargs['outfile']

	plot = quod.Plot()
	
	for i, seq in enumerate(sequences):
		plot.add(quod.What(seq, style=i, mode='entropy'))

	plot.render()

	plot.save(outfile)

#generates an entropy plot
entropy_quod(ALN1, ALN2, outfile='quod_entropy.png')
