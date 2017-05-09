<!--General layout and various text copied from https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/famXpander.md-->

# Documentation for script: _hvordan.py_

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

## Instructions

1. Move or symlink hvordan.py into the same directory where quod.py and tcblast.py are stored.

2. (Strongly recommended) Keep only the first six columns of each famXpander results table and move them all into the same directory. 
```bash
mkdir ../famXpander_trimmed
for DIR in `ls`
    do cut -f1-6 $DIR/psiblast.tbl > ../famXpander_trimmed/$DIR.tbl
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

    --p1d  directory containing _famXpander_ results (default: .)
    --p2d  directory containing _Protocol2_ results (default: .)
    -o     output directory (default: hvordan_out)
	-f     families to inspect, required if using --p2d on root _Protocol2_ directories
	-z     minimum Z-score (default: 15)
    -Z     maximum Z-score (default: None)
    -c     force redownloads/redraws/regenerations (not fully implemented)
	-r     resolution of plots in DPI (default: 100)
    -e     Email address the NCBI can contact (default: $ENTREZ_EMAIL if set)
	-i     inspects only alignments containing at least one of these accessions
    -p     inspects only this specific pair of accessions
