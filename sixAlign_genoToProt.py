#!/usr/bin/env python

"""
Get six alignments from one DNA string

Change DNA nucleotides to amino acid sequences
Altered version of Computational Genomics HW1 Q5.

Input: File of query sequences in FASTA formatted file. Query sequences
should be labeled with names.

Output: Query sequences with alignment label at end with string of proteins,
followed by a newline and protein sequence of query DNA sequence.

"""

import sys

# from comp-genomics-class notebooks --> FASTA
def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.iteritems():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa

def getComplement(dna):
    complement = ''
    for i in range(len(dna), 0, -1):
        if dna[i-1] == 'G':
            com = 'C'
        elif dna[i-1] == 'C':
            com = 'G'
        elif dna[i-1] == 'A':
            com = 'T'
        elif dna[i-1] == 'T':
            com = 'A'
        complement += com
    
    return complement

def getThreeMer(dna):
    kmer = []
    for i in range(len(dna)):
        if ((i + 2) < len(dna)):
            kmer.append(dna[i:i+3])
    return kmer

def threeAlign(dna):
    two = dna[1:]
    three = dna[2:]
    return [dna, two, three]

def sixAlignString(dna):
    allAlign = []
    complement = getComplement(dna)
    allAlign.extend(threeAlign(dna))
    allAlign.extend(threeAlign(complement))
    return allAlign

def sixAlignName(name, parseDict):
    seq = parseDict.get(name)
    frames = sixAlignString(seq)
    return frames

def allFrames(parseDict):
    allFrames = {}
    for name in parseDict:
        frames = sixAlignName(name, parseDict)
        count = 1
        for f in frames:
            allFrames[name + '.a' + str(count)] = f
            count += 1
    return allFrames

def printDictToFasta(fastaDic):
    fwrite = open("qstr", "w")
    sys.stdout = fwrite
    for key in fastaDic:
        value = fastaDic.get(key)
        if len(value) != 0:
            print '>', key
            i = 0
            while i < len(value):
                print value[i:i+70]
                i += 70
                # print ''


# start from HW1 Q5
geneticCode = { \
               'TTT':'F', 'CTT':'L', 'ATT':'I', 'GTT':'V',
               'TTC':'F', 'CTC':'L', 'ATC':'I', 'GTC':'V',
               'TTA':'L', 'CTA':'L', 'ATA':'I', 'GTA':'V',
               'TTG':'L', 'CTG':'L', 'ATG':'M', 'GTG':'V',
               'TCT':'S', 'CCT':'P', 'ACT':'T', 'GCT':'A',
               'TCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A',
               'TCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A',
               'TCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
               'TAT':'Y', 'CAT':'H', 'AAT':'N', 'GAT':'D',
               'TAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D',
               'TAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E',
               'TAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E',
               'TGT':'C', 'CGT':'R', 'AGT':'S', 'GGT':'G',
               'TGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
               'TGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G',
               'TGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G' }

def codonTranslate(s):
    s = s.upper()
    if s not in geneticCode:
        raise RuntimeError("No such codon: '%s'" % s)
    return geneticCode[s]

def stringTranslate(s):
    pstr = []
    for i in xrange(0, len(s)-2, 3):
        pc = codonTranslate(s[i:i+3])
        if pc == "Stop":
            return ''.join(pstr)
        else:
            pstr.append(pc)
    return ''.join(pstr)

# end from HW1 Q5

def convertToAminoAcid(fastaDic):
    proteinDict = {}
    for key in fastaDic:
        proteinDict[key] = stringTranslate(fastaDic.get(key))
    return proteinDict

# filename = sys.stdin.readline().rstrip()
filename = "ecoliQuerySeqs.fasta"
fastaFile = open(filename, 'r')
parseDict = parse_fasta(fastaFile)
fastaDic = allFrames(parseDict)
sixAlignProtein = convertToAminoAcid(fastaDic)
printDictToFasta(sixAlignProtein)