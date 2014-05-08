#!/usr/bin/env python
from collections import defaultdict as dd
import argparse

def parse_pop_file(fh):
    '''
    read in a headerless file specifying individual -> population mapping
    return a dict thats keys are pop names, values are lists of individuals
    '''
    pops = dd(lambda : [])
    indivs = []
    for l in fh:
        (indiv, pop) = l.split()
        pops[pop].append(indiv)
        indivs.append(indiv)
    return(pops, indivs)

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
    def __init__(self, fh):
        '''
        initialize starting variables
        '''
        self.file_handle = fh
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

def get_genepop_matrix(alleles_file, individuals):
    gpop_field = '{:03d}' * 2
    loci_names = []
    hapmap = dd(lambda : [])
    for i, locus in enumerate(Alleles_File(alleles_file)):
        if locus.num_haplotypes > 1:
            loci_names.append(str(i+1))
            for ind in individuals:
                if ind in locus.inds:
                    hapmap[ind].append(gpop_field.format(locus.inds[ind][0],locus.inds[ind][1]))
                else:
                    hapmap[ind].append(gpop_field.format(0,0))
    return(loci_names, hapmap)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given a pyRAD alleles file and individual -> populations assignment, generate a genepop file.")
    parser.add_argument("-a","--alleles", help="Your alleles file.",
                        type=argparse.FileType('r'), required=True)
    parser.add_argument("-p","--pops", help="Populations assignment file. Must be white-space delimited, 2 columns. column 1: individual column 2: population",
                        type=argparse.FileType('r'), required=True)
    parser.add_argument("-o","--out", help="The name of the genepop file.",
                        type=argparse.FileType('w'), required=True)
    args = parser.parse_args()
    
    out = args.out
    
    populations, individuals = parse_pop_file(args.pops)

    loci_names, haplotypes = get_genepop_matrix(args.alleles, individuals)
    
    out.write("haplotypes\n")
    out.write('\n'.join(loci_names) + '\n')
    for pop in sorted(populations):
        out.write("POP\n")
        for ind in populations[pop]:
            out.write(ind + ' , ')
            out.write(' '.join(haplotypes[ind]) + '\n')
