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

# READ IN FASTA FILES TO GET ALL DNA QUERY SEQUENCES
# ----------------------------------------------------------------------
# input: file name
# output: dict --> key (name of sequence), value (DNA sequence)
# ----------------------------------------------------------------------
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

      
# FIND SEQUENCE OF DNA AFTER START CODON AND BEFORE STOP CODON IN QUERY SEQ
# since only that sequence will be expressed as amino acids to create a protein
# ----------------------------------------------------------------------
def findStart(s):
    ''' Find the start codon in a given DNA sequence and output sequence 
    starting from there. '''
    found = False
    index = 0
    startIndex = 0
    while found == False:
        # find instance of start codon
        index = s.find('ATG', startIndex)
        # if start codon is not found
        if index < 0:
            return -1
        # check if start codon will be read in this particular alignment
        if (index % 3) == 0:
            found = True
            # cut off junk DNA before start codon
            s = s[index:]
            return s
        # set index of search to after incorrect start codon
        else:
            startIndex += index
            
def findEnd(s):
    ''' Find the stop codon in a given DNA sequence and output sequence before
    there. '''
    found = False
    index = 0
    startIndex = 0
    while found == False:
        # find instance of end codons
        index1 = s.find('TAA', startIndex)
        index2 = s.find('TAG', startIndex)
        index3 = s.find('TGA', startIndex)
        indexList = [index1, index2, index3]
        # if none of these stop codons are not found
        if max(indexList) == -1:
            return s
        # get the minimum index of stop codons
        index = min(i for i in indexList if i > 0)
        # check if stop codon will be read in this particular alignment
        if (index % 3) == 0:
            found = True
            # cut off junk DNA after stop codon
            s = s[:index + 3]
            return s
        # set index of search to after incorrect start codon
        else:
            startIndex += index
            
def startStopAll(dictionary):
    startStop = {}
    for key in dictionary:
        s = dictionary.get(key)
        startSeq = findStart(s)
        if startSeq != -1:
            endSeq = findEnd(startSeq)
            if endSeq != -1:
                startStop[key] = endSeq
    return startStop

def trimSeq(dictionary):
    ''' Remove any junk DNA and stop/stop codons from DNA sequence. '''
    updatedDict = {}    
    startStopDict = startStopAll(dictionary)    
    for k in startStopDict:
        seq = startStopDict.get(k)
        if seq[:3] == 'ATG':
            if (seq[-3:] == 'TAA' or seq[-3:] == 'TAG' or seq[-3:] == 'TGA'):
                updatedDict[k] = seq[3:-3]
            else:
                updatedDict[k] = seq[3:]
    return updatedDict       
            

# FIND ALL 6 FRAMES OF QUERY SEQUENCE
# ----------------------------------------------------------------------
def getComplement(dna):
    ''' Returns complement of a DNA string sequence. '''
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
    ''' Get three frames of a DNA sequence. '''
    two = dna[1:]
    three = dna[2:]
    return [dna, two, three]

def sixAlignString(dna):
    ''' Get all six frames of a DNA sequence. '''
    allAlign = []
    complement = getComplement(dna)
    allAlign.extend(threeAlign(dna))
    allAlign.extend(threeAlign(complement))
    return allAlign

def sixAlignName(name, parseDict):
    ''' Get six alignments for a sequence of a particular name. '''
    seq = parseDict.get(name)
    frames = sixAlignString(seq)
    return frames

def allFrames(parseDict):
    ''' Get all alignments of all key names. '''
    allFrames = {}
    for name in parseDict:
        frames = sixAlignName(name, parseDict)
        count = 1
        for f in frames:
            allFrames[name + '.a' + str(count)] = f
            count += 1
    return allFrames
    
# TRANSLATE DNA --> PROTEIN (SKIPPED DNA --> RNA)
# DNA --> RNA step was skipped to decrease unnecessary processes
# since exons and introns were not taken into consideration in this model
# ----------------------------------------------------------------------
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

# PRINT ALL OF THE ALIGNMENTS IN FASTA FORMAT
# ----------------------------------------------------------------------
def printDictToFasta(fastaDic):
    neatDict = {}
    for key in sorted(fastaDic.iterkeys()):
        value = fastaDic.get(key)
        if len(value) > 4:
            neatDict[key] = fastaDic.get(key)
            if len(value) > 4:
                print '>', key
                i = 0
            while i < len(value):
                print value[i:i+70]
                i += 70
    return neatDict

def action():
    #filename = sys.stdin.readline().rstrip()
    filename = 'oneBacteriaSeq.fasta'
    fastaFile = open(filename, 'r')
    parseDict = parse_fasta(fastaFile)
    fastaDic = allFrames(parseDict)
    trimDict = trimSeq(fastaDic)
    outputFileName = 'query_output.txt'
    sixAlignProtein = convertToAminoAcid(trimDict)
    sys.stdout = open('query_output.txt', "w")
    printDictToFasta(sixAlignProtein)
    return outputFileName

action()
    