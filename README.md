# hvordan & Co.

A suite of niche tools of variable usefulness

## hvordan
HTML Visualization Of Reasonable, Decent Alignment Networks

This script generates HTML reports of individual Protocol2 results.

### Dependencies:

 * matplotlib
 * NCBI blast
 * [curl](https://curl.haxx.se/)

### Instructions

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
## tcblast
Making cheap knockoffs of popular tools since 2016!

Tcblast is not currently standalone.

## quod
Questionable Utility Of Doom

Makes average hydropathy graphs from sequences and sequence-containing files
