<!--General layout and various text copied from https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/famXpander.md-->

# Documentation for script: _haystack.py_

## Summary
This script performs a crude, Smith-Waterman alignment of two sequences based on their TMSs.

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
