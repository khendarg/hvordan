**Note**: In light of Microsoft's recent acquisition of Github, this project has migrated to https://www.gitlab.com/khendarg/hvordan . This repository will remain available on Github for another month.

# hvordan & Co.

A suite of niche tools of variable usefulness

## hvordan
HTML Visualization Of Reasonable, Decent Alignment Networks

This script generates HTML reports of individual Protocol2 results.

## Summary
This script generates HTML reports with hydropathy plots and representations of TCDB BLAST hits (mostly replicating the [TCDB BLAST tool](http://www.tcdb.org/progs/blast.php)) for _Protocol2_ results. 
This can be done in bulk on entire sets of _Protocol2_ results, on specific genes, on specific pairs of genes, or on specific ranges of GSAT Z-scores.

## Dependencies
The following programs need to be available in your path for this program to run properly:

1. **_blast+ 2.4.0 to 2.6.0_**  
Other versions of blast may require minor adaptations. 
Visit the
 [download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). 

2. **_TCDB protein database_**  
Given that blastdbcmd runs locally, the TCDB database must be available locally through the environment variable _$BLASTDB_. 
If possible, use [extractFamily.pl](https://github.com/SaierLaboratory/TCDBtools/blob/master/scripts/extractFamily.pl) to download the TCDB BLAST database:
```bash
extractFamily.pl -i all -f blast
```
Otherwise, you can manually download [all TCDB sequences](http://www.tcdb.org/public/tcdb) from the TCDB and run ```makeblastdb -in tcdb.fasta -out tcdb -dbtype prot```, but this script was not tested with such manually downloaded BLAST databases. 

3. **_Python 2.7+_**  
Visit the [official website](https://www.python.org/). 
This program was not tested with more recent versions of Python but was implemented with some forward compatibility.

4. **_Matplotlib 1.5.1-2.0.0+_**  
Visit the [official website](https://matplotlib.org/).
While not required for hvordan.py itself, Matplotlib is required for the graph plotting modules, and automating graph generation is the entire point of the script.

5. **_Biopython 1.70+_**
Visit the [official website](http://biopython.org/).
For now, index-finding for the ABCD plots will rely on pairwise2 on account of indices not being available on report.tbl.

## Instructions

1. Move or symlink hvordan.py into the same directory where quod.py and tcblast.py are stored.

2. (Strongly recommended) Keep only the first twelve columns of each famXpander results table and move them all into the same directory. 
```bash
mkdir ../famXpander_trimmed
for DIR in `ls`
    do cut -f1-12 $DIR/psiblast.tbl > ../famXpander_trimmed/$DIR.tbl
done
```
Failing that, symlink all of the table files into the same directory. Do note that the bottleneck in report generation is parsing these tables. **Avoid this if you value your time!**
```bash
mkdir ../famXpander_flat
for DIR in `ls`
	do ln -s $DIR/psiblast.tbl ../famXpander_flag/$DIR.tbl
done
```
At the moment, it is required that all table files be in the same directory.
3. (Optional) Set the environment variable $ENTREZ\_EMAIL to an address the NCBI can contact if you send too many requests. Otherwise, just use the ```-e``` argument every time. 
```bash
echo 'export ENTREZ_EMAIL=someone@example.com' >> ~/.profile
source ~/.profile
```
4. Make sure the TCDB BLAST database is installed either in the current directory or in $BLASTDB. If the database is not installed in either of these places, do this:
```bash
wget http://www.tcdb.org/public/tcdb -O tcdb.fasta
makeblastdb -in tcdb.fasta -out tcdb -dbtype prot
```
5. Run hvordan.py. (-h and --help provide more detailed documentation)
```bash
hvordan.py --p1d famXpander_trimmed --p2d 1.X.1_vs_2.Y.1/1.X.1_vs_2.Y.1 
```

## Command line options
The following options are available. 
You can also run the script without arguments (or with -h or --help) to display the options:

`--p1d`  directory containing _famXpander_ results (default: .)
`--p2d`  directory containing _Protocol2_ results (default: .)
`-o`     output directory (default: hvordan_out
`-f`     families to inspect, required if using --p2d on root _Protocol2_ directories
`-z`     minimum Z-score (default: 15)
`-Z`     maximum Z-score (default: None)
`-c`     force redownloads/redraws/regenerations (not fully implemented)
`-r`     resolution of plots in DPI (default: 100)
`-m`     maximum BLAST hits for the TCBLAST portion (default: 50)
`-e`     Email address the NCBI can contact (default: $ENTREZ_EMAIL if set)
`-i`     inspects only alignments containing at least one of these accessions
`-p`     inspects only these specific pairs of accessions

## tcblast
Making cheap knockoffs of popular tools since 2016!

Tcblast is not currently standalone.

## quod
Questionable Utility Of Doom

Makes average hydropathy graphs from sequences and sequence-containing files

## Summary
This script generates HTML reports with hydropathy plots and representations of TCDB BLAST hits (mostly replicating the [TCDB BLAST tool](http://www.tcdb.org/progs/blast.php)) for _Protocol2_ results. 
This script generates average hydropathy plots with arbitrary resolution in a variety of commonly used formats. 
This can be done on multiple sequences.
This tool mostly replicates [WHAT](http://biotools.tcdb.org/barwhat2.html).

## Dependencies
The following programs need to be available in your path for this program to run properly:

1. **_Python 2.7+_**  
Visit the [official website](https://www.python.org/). 
This program was not tested with more recent versions of Python but was implemented with some forward compatibility.

2. **_Matplotlib 1.5.1-2.0.0+_**  
Visit the [official website](https://matplotlib.org/).

## Command line options
The following options are available. 
You can also run the script without arguments (or with -h or --help) to display the options:

`-f` force all inputs to be interpreted as filenames (only useful for all-caps filenames with no extension)
`-s` force all inputs to be interpreted as sequences
`-v` verbose output
`-d` directory to store graphs in (recommended only with autogenerated filenames) (default: .)
`-o` filename of graph relative to the value of -d (default: (autogenerated))
`-q` disables automatic opening of graphs
`-a` viewer to be used for opening graphs
`-t` file format: {eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff} (default: png)
`-r` resolution in dpi (default: 80)
`-l` graph title
   ` (positional) sequences or filenames
`-m` manually add TMSs.

## Supported annotation features

### Manual TMS specification

TMS sets for different sequences are specified by space-separated comma-separated ranges. The following example specifies a TMS spanning residues 1-20 for the first/only sequence in seq1.fa:

```
./quod.py seq1.fa -m 1-20
```

This specifies TMSs for different sequences in seqcassete.faa:

```
./quod.py seqcassette.faa -m 1-20,30-50 3-23
```

This specifies TMSs for the second sequence but leaves the first to HMMTOP:

```
./quod.py seq1.fa seq2.fa -m skip 1-20,30-50,60-80
```

This specifies several green TMSs. Note the lack of quotes around the hexadecimal specification, which are unnecessary on a local build of bash 4.4.12(1)-release and which may be necessary on other bash builds:

```
./quod.py seq1.fa -m 1-20:green,30-50:#00ff00,60-80:darkgreen
```

### Vertical bars

Use space-separated x-values with `-b` to get vertical bars.

### Wedges/arrowheads

Use space-separated comma-separated `x,length` pairs with `-w` to get wedges/arrowheads. Negative lengths result in left-pointing wedges. Use `-W` to automatically draw bars for each wedge.
<!--General layout and various text copied from https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/famXpander.md-->

## haystack
A stack of needles may overflow

Compares two sequences based on CGAT (<- GSAT) Z-scores of TMSs (<- HMMTOP) using a linear gap cost Smith-Waterman implementation

## Summary
This script performs a crude Smith-Waterman alignment of two sequences based on their TMSs.

## Dependencies
The following programs need to be available in your path for this program to run properly:

1. **_Python 2.7+_**  
Visit the [official website](https://www.python.org/). 
This program was not tested with more recent versions of Python but was implemented with some forward compatibility.

## Command line options
The following options are available. 
You can also run the script without arguments (or with -h or --help) to display the options:

`fasta1` FASTA file of one record. Will not be shuffled
`fasta2` FASTA file of one record. Will be shuffled
`-s` number of shuffles to perform (default:2000}
`-v` verbose output
`-g` gap opening cost for needle (default:10.0)
`-e` gap extension cost for needle (default:0.5)
`-t` TMS gap cost for haystack (default:1.0)
