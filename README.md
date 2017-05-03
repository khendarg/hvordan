# hvordan
HTML Visualization Of Reasonable, Decent Alignment Networks

This script generates HTML reports of individual Protocol2 results.

## Dependencies:

 * [quod](https://www.github.com/khendarg/quod)
 * * matplotlib
 * [tcblast](https://www.github.com/khendarg/tcblast)
 * * NCBI blast
 * [curl](https://curl.haxx.se/)

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
3. (Optional) Set the environment variable $ENTREZ\_EMAIl to an address the NCBI can contact if you send too many requests. Otherwise, just use the ```-e``` argument every time. 
```bash
echo 'export ENTREZ_EMAIL=someone@example.com' >> ~/.profile
source ~/.profile
```
4. Run hvordan.py. (-h and --help provide more detailed documentation)
```bash
hvordan.py --p1d famXpander_trimmed --p2d 1.X.1_vs_2.Y.1/1.X.1_vs_2.Y.1 
```
