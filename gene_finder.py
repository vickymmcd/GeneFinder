# -*- coding: utf-8 -*-
"""
Gene Finder finds amino acid sequences that are likely coded by the specified
DNA input

@author: Vicky McDermott

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        added doctest for input of empty string and random letter
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('')
    ''
    >>> get_complement('L')
    ''
    """
    if(nucleotide == 'A'):
        return 'T'
    elif(nucleotide == 'T'):
        return 'A'
    elif(nucleotide == 'C'):
        return 'G'
    elif(nucleotide == 'G'):
        return 'C'
    else:
        return ''


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        added doctest for input of empty string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("")
    ''
    """
    complement = ''
    i = len(dna)-1
    while i >= 0:
        complement = complement + get_complement(dna[i])
        i = i-1
    return complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        added doctest for input of empty string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("")
    ''
    """
    i = 3
    while i < len(dna):
        if dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA':
            stop_codon = i
            return dna[:stop_codon]
        i = i+3
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added doctest to ensure nested ORFs are not included
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGAAAATGAAATAG")
    ['ATGAAAATGAAA']
    """
    i = 0
    myORFS = []
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            myORF = rest_of_ORF(dna[i:])
            myORFS.append(rest_of_ORF(dna[i:]))
            i = i + len(myORF)
        else:
            i = i + 3
    return myORFS


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added doctest to ensure nested ORFs are not included
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs_both_strands("ATGCGAATGTAG")
    ['ATGCGAATG']
    """
    myOrfs = []
    i = 0
    for i in range(0, 3):
        myOrfs = myOrfs + find_all_ORFs_oneframe(dna[i:])
    return myOrfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added doctest to ensure empty string does not break it
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("")
    []
    """
    myOrfs = []
    myOrfs = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return myOrfs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

        dna: a DNA sequence
        returns: string containing longest ORF on both strands

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longOrfs = []
    longest = ''
    longOrfs = find_all_ORFs_both_strands(dna)
    for i in range(len(longOrfs)):
        if len(longOrfs[i]) > len(longest):
            longest = longOrfs[i]
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longestlen = 0
    for i in range(num_trials):
        dna = shuffle_string(dna)
        num = len(longest_ORF(dna))
        if num > longestlen:
            longestlen = num
    return longestlen


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i = 0
    aminos = ''
    while i < len(dna):
        if len(dna[i:i+3]) == 3:
            amino_acid = aa_table[dna[i:i+3]]
            aminos = aminos + amino_acid
        i = i + 3
    return aminos


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aminoseqs = []
    threshold = longest_ORF_noncoding(dna, 1500)
    myOrfs = find_all_ORFs_both_strands(dna)
    for i in range(len(myOrfs)):
        if len(myOrfs[i]) > threshold:
            aminoseqs.append(coding_strand_to_AA(myOrfs[i]))
    return aminoseqs


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
