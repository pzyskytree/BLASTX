
README

# BLASTX

Enclosed in BLASTX folder are:

Programs (8):
BLAST.py
BWT.py
data_preprocess.py
dnaToProtein.py
extending.py
fm_index.py
singlescore.py
suffix_tree.py

FASTA files (3)
oneBacteraSeq.fasta
query_output.txt
output_result.txt

Misc files (2)
FM.txt
score_matrix

Folder (1)
Output*

*There are 9 different kinds of output possible from our implementation of BLAST, depending on what kind of settings were selected for data structure utilized and what kind of setting were selected for the extension step of BLAST. They have been put into this folder.

To make use of these programs, run the following commands:
python dnaToProtein.py
oneBacteriaSeq.fasta
python BLAST.py

After which, there will be a prompt for the user:
Please type your command:

The command will take two arguments, each with some options:
1. -l, -s, or -f
2. -d, -s, -g

The first are options on what data structures to use when implementing BLAST.py:
-l: look up table (online)
-s: suffix tree (offline)
-f: FM index (offline)

The second are options on how the program should operate when executing the extending part of the BLAST algorithm:
-d: extend with drop threshold, but no gap
-s: extend with score distance, but no gap
-g: extend with gap

Note: The input file must only have one query sequence in it, in FASTA format.

Sources for Code:

The implementation for a suffix tree (suffix_tree.py) came from: https://github.com/kvh/Python-Suffix-Tree

Some of dnaToProtein.py’s implementation came from/were modified from Ben Langmead:
http://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/FASTA.ipynb
https://piazza-resources.s3.amazonaws.com/j6dmgkudvu2bc/j7wfc3ivs851sd/hw1_sol.pdf?AWSAccessKeyId=AKIAIEDNRLJ4AZKBW6HA&Expires=1512800497&Signature=l4yLU8F9e5s4J%2FJMfxmcNkJzo%2F8%3D






