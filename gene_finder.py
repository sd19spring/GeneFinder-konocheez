# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Jerry

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
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this

    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'C':
        return 'G'



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this

    #DNA = list(dna) #turns string into list
    DNA_pair = [] #creates new list for pairs
    i = 0
    for i in range(len(dna)):
        new_pair = get_complement(dna[i]) #fills in new list with pairs
        DNA_pair.insert(0, new_pair)
        i += 1
    """for base_pair in DNA_pair:
        n = (len(DNA_pair) - 1) #gets length of DNA_pair
        if n > 0:
            rev_DNA_pair = DNA_pair #creates new array for reverse DNA pairs
            rev_DNA_pair[0] = DNA_pair[n] #assigns first list item in reverse to final list item in DNA_pairs
            n = n - 1"""
    delimiter = ''
    s = delimiter.join(DNA_pair)
    return s

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this

    #SINCE IT ASS-U-MES IT BEGINS WITH START CODON,
    #I WON'T CODE TO START AT START CODON


    for i in range(0, len(dna), 3):
        if ((dna[i:i+3] == 'TAG') or (dna[i:i+3] == 'TGA') or (dna[i:i+3] == 'TAA')):
            return dna[:i]
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
    >>> find_all_ORFs_oneframe('ATGCATGAATGTAGATAGATGTGCCC')
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this

    all_ORFs_oneframe = []

    i = 0
    while i < len(dna):
        if (dna[i:i+3] == 'ATG'): #if the i:i+3 sequence is ATG, print rest of ORF
            all_ORFs_oneframe += [rest_of_ORF(dna[i:])]
            i += len(rest_of_ORF(dna[i:]))
        else:
            i += 3
    return all_ORFs_oneframe

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs('ATGCATGAATGTAG')
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this

    find_all_ORFs_notfunction = [] #list of indices at different starting i's
    for i in range(3):
        find_all_ORFs_notfunction.extend(find_all_ORFs_oneframe(dna[i:]))

    return find_all_ORFs_notfunction



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    #!/usr/bin/env python3
    import doctest
    doctest.testmod(verbose=True)
