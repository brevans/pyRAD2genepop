#!/usr/bin/env python
from collections import defaultdict as dd

def parse_pop_file(fn):
    '''
    read in a headerless file specifying individual -> population mapping
    return a dict thats keys are pop names, values are lists of individuals
    '''
    pops = dd(lambda : [])
    with open(fn, 'r') as fh:
        for l in fh:
            (indiv, pop) = l.split()
            pops[pop].append(indiv)
    return(pops)

class Locus(object):
    '''
    container and helper functions for parsing loci in an alleles file
    '''
    def __init__(self):

        self.inds = dd(lambda : [])
        #current number of haplotypes at this locus
        self.num_haplotypes = 0
        #haplotype -> hapnum mapping
        self.haps = {}

    def add_ind(self, l):
        (tmp_ind, hap) = l.split()
        if hap not in self.haps:
            self.num_haplotypes += 1
            self.haps[hap] = self.num_haplotypes
        ind = tmp_ind[1:-2]
        self.inds[ind].append(self.haps[hap])

class Alleles_File(object):
    '''
    collections of alleles, phased, for individuals
    -sample names are appended with r"_[01]" to specify haplotype 0 or 1
    -spaces separate sample IDs and sequence
    -loci are separeated by a line at the bottom that specifies the positions 
       of SNPs with either a * or -
    allele lines look like:
    >sample1_0   TCATACGTGACTGACTGACxxxxTCATACGTGACTGACTGAC
    >sample1_1   TCATAAGTGACTGACTGACxxxxTCATACGGGACTGACTGAC
    //                -                        *           |
    '''
    def __init__(self, fn):
        '''
        '''
        self.file_handle = open(fn, 'r')
        self.curr_loc = Locus()
        self.next_loc = self.get_next_locus()

    def get_next_locus(self):
        '''
        '''
        nl = Locus()
        try:
            currline = self.file_handle.next()
            while not currline.startswith("//"):
                nl.add_ind(currline)
                currline = self.file_handle.next()
        finally:
            if nl.num_haplotypes >0:
                return nl
            else:
                return None

    def __iter__(self):
        return self
    
    # Python 3 compatibility
    def __next__(self):
        return self.next()

    def next(self):
        '''
        '''
        self.curr_loc = self.next_loc
        if self.curr_loc != None:

            self.next_loc = self.get_next_locus()
            return self.curr_loc
        else:
            raise StopIteration()

if __name__ == '__main__':
    popfile = 'pops.txt'
    alleles_file = 'sample.alleles'

    num_poly = 0
    for i, locus in enumerate(Alleles_File(alleles_file)):

        if locus.num_haplotypes > 1:
            num_poly += 1
            print("locus#:", i+1, "numpoly:", num_poly, "individuals:", locus.inds)
    print("done.")
