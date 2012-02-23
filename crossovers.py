#!/usr/bin/python2.6
import sys
import os
import string
import re
from optparse import OptionParser
import numpy as np
from Interval import Interval

""" generate crossover intervals that represent meiosis events to generate haploid gamete(s) """

def main():
    usage = "usage: %prog [options] chromInfo.bed"
    parser = OptionParser(usage)
    parser.add_option("--morgan", type="float", dest="morgan", default=100000000, help="physical length of 1 Morgan (default 100Mbp)")
    parser.add_option("--paternal", type="string", dest="paternal", help="fasta representing paternal haplotype")
    parser.add_option("--maternal", type="string", dest="maternal", help="fasta representing maternal haplotype")
    (options, args)=parser.parse_args()

    if options.maternal == None:
        sys.stderr.write("please provide maternal fasta file!\n")
        exit(1)
    if options.paternal == None:
        sys.stderr.write("please provide paternal fasta file!\n")
        exit(1)

    crossoverfh=open("crossover.log", 'w')

    chromSizes=[]
    chromfh=open(args[0])
    for line in chromfh:
        if '#' in line: continue
        fields=line.strip().split('\t')
        (chrom, start, end)=fields[0:3]
        chromSizes.append(int(end))

    lambdas=map(lambda x: x/options.morgan, chromSizes)
    #print lambdas
    for i in range(len(lambdas)):
        j=i+1
        chrom="chr"+str(j)
        lam=lambdas[i]
        #get the number of crossovers by drawing from poisson
        n=np.random.poisson(lam, 1)
        ncrossovers=n[0]
        #pick haplotype to pass
        flip=np.random.uniform(0,1)
        if flip <= .5: 
            haplotype="maternal"
        else:
            haplotype="paternal"
        #now get the crossover points
        if ncrossovers ==0:
            intervalstring="\t".join([haplotype, chrom, "+",'0', str(chromSizes[i]), str(ncrossovers),',',','])
            crossoverfh.write(chrom + " " + haplotype + " " + str(ncrossovers) + " " + str(chromSizes[i]) + "\n")
        else:
            chromsize=chromSizes[i]
            points=np.random.random_integers(1, chromsize, ncrossovers)
            points.sort()
            pointstrings=[]
            for p in points: pointstrings.append(str(p))
            outstring="\t".join(pointstrings)
            crossoverfh.write(chrom + " " + haplotype + " " + str(ncrossovers) + " " + ",".join(pointstrings) + "\n")
            if len(points) %2 == 0:
                #sys.stderr.write(haplotype + " " + chrom + " " + outstring +"\n")
                starts=list(points[1::2])
                ends=list(points[::2])
                ends.append(int(chromsize))
                starts.insert(0,0)
                startstrings=",".join(map(lambda x: str(x), starts))
                endstrings=",".join(map(lambda x: str(x), ends))
                intervalstring="\t".join([haplotype, chrom, "+",str(starts[0]), str(ends[-1]),str(ncrossovers), startstrings, endstrings])
               
            else:
                #sys.stderr.write(haplotype + " " + chrom + " " + outstring +"\n")
                starts=list(points[1::2])
                ends=list(points[0::2])
                starts.insert(0,0)
                startstrings=",".join(map(lambda x: str(x), starts))
                endstrings=",".join(map(lambda x: str(x), ends))
                #print starts
                #print ends
                intervalstring="\t".join([haplotype, chrom, "+",str(starts[0]), str(ends[-1]),str(ncrossovers), startstrings, endstrings])
            #print intervalstring
            intervalobj=Interval(intervalstring)
            for (name, start,end) in intervalobj.yieldStartEnds():
                print haplotype, chrom, start, end
        #print

if __name__ == "__main__":
    main()
