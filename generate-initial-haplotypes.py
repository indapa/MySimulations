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

def sampleSNPFreq(N, sfs,L):
    """ snp freq are in proportion to the SFS 1/x where x is the sfs class 1/(1/N), 1/(2/N) ...  """
    sys.stderr.write("sampling SNP freq in proportion to neutral SFS ....\n")
    genotypes_list=[]
    freqProp= [float(x)/N for x  in range(1,N)]
    total=0
    ss=np.sum(np.array(sfs))
    #print len(sfs), len(freqProp), ss
    for i in range(len(sfs)):
        #print sfs[i], freqProp[i]
        j=i+1
        #print sfs[i], j
        Ei_allelefrq= np.random.uniform(0, freqProp[i], int(sfs[i])).round(decimals=3)
        print j, len(Ei_allelefrq)
        Ei_allelefrq_str=",".join(map(str, Ei_allelefrq))
        #print j, freqProp[i], Ei_allelefrq_str
        total+=len(Ei_allelefrq)
        #print len(Ei_allelefrq)
        #exp_val=np.mean(Ei_allelefrq)
        #print sfs[i],j
        for frq in Ei_allelefrq:
            # now to assign genotypes, we draw from a binomial in proportional to the sampled allele frq, which was sampled in proportion to neutral SFS
            genotypes=list( np.random.binomial(2, frq, size=N/2) )
            genotypes_list.append(genotypes)
            
    #print total
    return genotypes_list
""" given theta, N (i.e. number of chroms), and L (length of region), return the expected number of segregating sites using
    Watterson's theta For information on Watterson's estimate see http://en.wikipedia.org/wiki/Watterson_estimator   """
def expectedSegsites(theta, N):
    aK=0.
    for i in range(1, N):
        aK+=float(1./i)
    eS=int(theta  * aK)
    return eS

def returnNeutralSFS (thetaW, N):
    """
    return a netural spectrum given Watterson's estimate theta and #chromosomes
    For information on Watterson's estimate see http://en.wikipedia.org/wiki/Watterson_estimator
    """
    sys.stderr.write("calculating neutral SFS ...\n")
    neutralCounts=[]
    for i in range(1,N):
        neutralCounts.append( int ( float(thetaW)/float(i)) )
    return neutralCounts



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


""" generate founder haplotypes given an 2bit file of sequence """

def main():
    usage = "usage: %prog [options] file.2bit"
    parser = OptionParser(usage)
    parser.add_option("--theta", type="float", dest="theta", help="per-base theta value (default .001)", default=.001)
   
    parser.add_option("--chrombed", type="string", dest="chrombed", help="bed file with chromsizes", default="chromInfo.bed")
    parser.add_option("--chrom", type="string", dest="chrom", help="chromosome to simulate haplotypes for", default='chr22')
    parser.add_option("--N", type="int", dest="N", help="number of founder chromosomes to simulate", default=50)
    (options, args)=parser.parse_args()

    #if options.hap == None:
    #    sys.stderr.write("pleas provide --hap value!\n")
    #    exit(1)
    #bed file to wrtie record of mutation
    #mutationfh=open(options.hap+"_"+"mutations.bed", 'w')

    pedfilename='sim.ped'
    mapfilename='sim.map'

    pedfh=open(pedfilename, 'w')
    mapfh=open(mapfilename, 'w')

    tbf=args[0]
    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( tbf ) )
    except:
        sys.stderr.write("unable to open twobit file!\n")

    chromsizes={}#key chromName, value is length
   
    #first read the chromsizes in to dict
    chromfh=open(options.chrombed, 'r')
    for line in chromfh:
        (chrom, start, end)=line.strip().split('\t')
        end=int(end)
        chromsizes[chrom]=end
    
    for chr in chromsizes.keys():
        if chr != options.chrom: continue
        L=chromsizes[chr]
        thetaScaled=int( float(L)*options.theta )

        #print thetaScaled
        #now the number of mutations is found via Watterson's formula
        #segsites=expectedSegsites(thetaScaled,options.N )
        #print segsites
        
        #now get the  netural SFS
        sfs=returnNeutralSFS(thetaScaled, options.N)
        #number of seg sites is the count of each bin in the SFS
        segsites=np.sum(np.array(sfs))
        #print segsites

        #now conditioned upon the numebr of snps (segsites), place them randomly along the length of the chromosome
        # this results in a non-unique set of sites with the code below .... how to fix this???
        points=list(np.random.random_integers(1, high=L, size=segsites) )
        
      
        #now sample allele freq in proportion to the neutral SFS
        genotypes_list=sampleSNPFreq(options.N, sfs,L)

        
        #make a numpy 2d matrix
        genotype_matrix=[]
        for  genos in genotypes_list:
            
            #print chr, refAllele, altAllele, genos
            genotype_matrix.append(genos)
        #print chred), segsites, len(points),sfsSum

        #make a numpy 2-d array where nrows = #segsites ncols= N/2 of desired haplotypes
        genotype_nd_matrix=np.array(genotype_matrix)
        #print np.shape(genotype_nd_matrix)

        #list of tuples with [ (pos, ref, alt) ... ]
        site_and_alleles=[]

        #now assign ref and alt alleles according to Jukes-Cantor ....
        sys.stderr.write("assigning ref and alt alleles according Jukes-Cantor and writing out *.map file...\n")
        for i in range(len(points)):
            (start, end)= ( int(points[i])-1, int(points[i]) )
            assert( end > start ),"end greater than start!"

            try:
                sequence=twobit[chr][start:end]
                sequence=sequence.upper()
            except:
                sys.stderr.write("unable to fetch sequence from 2bit file!\n")
            refAllele=sequence
            refAllele=refAllele.upper()
            altAllele=JC_rates(refAllele)
            site_and_alleles.append ( (points[i], refAllele, altAllele) )
            mapstring="\t".join( [ chr, "snp"+str(i+1), '0', str(points[i]) ])
            mapfh.write(mapstring+"\n")

        #now iterate over the sample unphased genotypes by taking column slices for each individual
        # the columns represent the h_th individual's genotypes for all the segregating sites
        sys.stderr.write("generating fasta sequences of simulated haplotypes and writing genotypes to *.ped file ...\n")
        for h in range(0, options.N/2):
            sex=''
            j=h+1
            if h%2 == 0:
                sex='1'
            else:
                sex='2'
            maternal_haplotype=list(twobit[chr][0:L])
            maternal_name="_".join([chr, 'maternal', "sample"+str(j)])
            
            paternal_haplotype=list(twobit[chr][0:L])
            paternal_name="_".join([chr, 'paternal', "sample"+str(j)])

            
            #we take the h-th column, corresponding the genotypes of the h-th individual (this is the unphased haplotype)
            unphased_genotypes=genotype_nd_matrix[:,h]

            
            #now for ith segregating site get its genotype and correpsonding ref and alt allele
            genotypes=[]
            for i in range(len(unphased_genotypes)):
                #print unphased_genotypes[i], site_and_alleles[i]
                (pos, ref, alt) = site_and_alleles[i]
                maternal_genotype=''
                paternal_genotype=''
                # if the genotype is het, decide which haplotype to assign it to
                if unphased_genotypes[i] == 1:
                    if random.uniform(0,1) < .5:
                        maternal_haplotype[pos-1]=alt
                        maternal_genotype='1'
                        paternal_genotype='0'
                    else:
                        paternal_haplotype[pos-1]=alt
                        paternal_genotype='1'
                        maternal_genotype='0'
                elif unphased_genotypes[i] == 2:
                    maternal_haplotype[pos-1]=alt
                    paternal_haplotype[pos-1]=alt
                    paternal_genotype='1'
                    maternal_genotype='1'
                elif unphased_genotypes[i] == 0:
                    paternal_genotype='0'
                    maternal_genotype='0'
                    pass
                else:
                    sys.stderr.write("assigned genotype is greater than 2! " + unphased_genotype[i] + "\n")
                genotypes.append( " ".join([paternal_genotype, maternal_genotype]))
            genotype_string="\t".join(genotypes)

            pedinfo="\t".join([ str(j), str(j), '0', '0', sex, '1' ])
            pedstring="\t".join([pedinfo, genotype_string])
            pedfh.write(pedstring + '\n')

            writefasta("".join(maternal_haplotype), chr + " maternal", maternal_name+".fa")
            writefasta("".join(paternal_haplotype), chr + " paternal", paternal_name+".fa")





                    
if __name__ == "__main__":
    main()
