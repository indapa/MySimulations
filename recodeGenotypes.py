#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
import collections
def main():
    usage = "usage: %prog [options] file.ped"
    parser = OptionParser(usage)
    parser.add_option("--pseudomap", type="string", dest="pseudomap", help="recode the file.ped genotypes to alleles in the pseudomap")

    (options, args)=parser.parse_args()
    pedfile=args[0]
    
    pedfh=open(pedfile,'r')
    mapfh=open(options.pseudomap,'r')

    #pseudomap contains the ref/alt bases based on making haplotypes from teh cosi haplotype simulation
    #and simulated reference genome
    
    Alleles=collections.namedtuple('Alleles', ['chrom','start','end','ref','alt'])
    loci=[]
    for line in mapfh:
        (chrom,start,end,reference,alternate)=line.strip().split('\t')
        loci.append( (chrom,start,end,reference,alternate) )

    LociNamedList=map(Alleles._make, loci)
    
    for line in pedfh:
        fields=line.strip().split('\t')
        genotypes=fields[6::]
    
        newgenotypes=[]
        pedinfo=fields[0:6]
        if len(genotypes) != len(LociNamedList):
            sys.stderr.write("unequal number of genotypes and loci!\n")
            print pedinfo
            sys.exit(1)

        for i in range(len(genotypes)):
            (a1,a2)=genotypes[i].split(' ')
            #print a1, a2, LociNamedList[i].ref, LociNamedList[i].alt
            
            if a1 == '1':
                a1=LociNamedList[i].alt
            else:
                a1=LociNamedList[i].ref

            if a2 == '1':
                a2=LociNamedList[i].alt
            else:
                a2=LociNamedList[i].ref
            newgenotype=' '.join([a1,a2])
            newgenotypes.append(newgenotype)
        newfields=pedinfo+newgenotypes
        print "\t".join(newfields)
            
                         


if __name__ == "__main__":
    main()
