#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
import bx.seq.twobit
import numpy as np
from common import *
import random



def JC_rates(refbase):
    if refbase=='A':
        return random.choice('CGT')
    elif refbase=='C':
        return random.choice('AGT')
    elif refbase == 'G':
        return random.choice('ACT')
    elif refbase == 'T':
        return random.choice('ACG')
    elif refbase == 'N':
        return 'N'
    else:
        sys.stderr.write(refbase + " not in DNA alphabet!\n")
        return 'N'


""" generate founder sequence haplotypes given an 2bit file of sequence based on cosi simulated haplotypes"""

def main():
    usage = "usage: %prog [options] file.2bit"
    parser = OptionParser(usage)
    parser.add_option("--cosihap", type="string", dest="hapfile", help="file with cosi simulated haplotypes")
    parser.add_option("--cosipos", type="string", dest="posfile", help="file with variant position information")

  
    (options, args)=parser.parse_args()
    tbf=args[0]
    
    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( tbf ) )
    except:
        sys.stderr.write("unable to open twobit file!\n")


    cosiposfh=open(options.posfile, 'r')

    altbases=[]
    altallelenum=[]
    positions=[]
    for line in cosiposfh:



       #refbase= get from tbf

        if 'SNP' in line:
            continue
        else:
            (snp, chrom, pos, allele1, freq1, allele2, freq2)=line.strip().split('\t')

            

            start=int(pos)-1
            end=int(pos)
            positions.append(start)
            try:
                sequence=twobit[chrom][start:end]
                sequence=sequence.upper()
            except:
                error="unable to fetch sequence from 2bit file!: " + chrom + " " + pos
                sys.stderr.write(error + "\n")
                exit(1)

            refAllele=sequence
            refbase=refAllele.upper()

            altbase=JC_rates(refbase)
            altbases.append(altbase)


            if freq1 < freq2:
                altAllele=allele1
            else:
                altAllele=allele2
            #print chrom, pos, refbase, altbase, altAllele, allele1, freq1, allele2, freq2
            altallelenum.append(altAllele)

    #print len(positions)

    happosfh=open(options.hapfile, 'r')

    for line in happosfh:
        (hapid, popid, genostr) = line.strip().split('\t')
        genotypes=genostr.split(' ')
        print len(genotypes)

        hapsequence=list( twobit[chrom][0:1000000] )
        for i in range(0,len(genotypes) ):
            print positions[i], altallelenum[i], altbases[i], genotypes[i]
            if altallelenum[i] == genotypes[i]:
                print "haplotype has alt allele at this position"
                hapsequence[ positions[i] ] = altbases[i]

        seqname="hap"+hapid
        writefasta("".join(hapsequence), seqname, "hap"+hapid+".fa")





                    
if __name__ == "__main__":
    main()
