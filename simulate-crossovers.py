#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
import numpy as np
import bx.seq.twobit
from Interval import Interval
from common import *

""" given a bedstring and bx.seq.twobit object, extract the sequence from the intervals parsed from the bedstring """
def twoBitExtract(bedstring, twobit):
    (chr, start, end, name)=bedstring.split('\t')
    
    start=int(start)
    end=int(end)
 
    assert( end > start ),"end less  than start!"
    sequence=twobit[chr][start:end]
    return sequence.upper()


""" generate crossover intervals that represent meiosis events to generate haploid gamete(s)
    given the maternal and paternal 2bit file. """
def main():
    usage = "usage: %prog [options] chromInfo.bed"
    parser = OptionParser(usage)
    parser.add_option("--morgan", type="float", dest="morgan", default=100000000, help="physical length of 1 Morgan (default 100Mbp)")
    parser.add_option("--grandpaternaltwobit", type="string", dest="paternaltbf", help="2bit file representing (grand)paternal haplotype")
    parser.add_option("--grandmaternaltwobit", type="string", dest="maternaltbf", help="2bit file representing (grand)maternal haplotype")
    parser.add_option("--samplename", type="string", dest="samplename", help="sample name for which you are generating a gamete for", default="gamete")
    parser.add_option("--meiosis", type="int", dest="nmeiosis", help="meiotic event id (default 1) ", default=1)
    parser.add_option("--chromSize", type="int", dest="chromSize", help="length of chromosome", default=1000000)
    parser.add_option("--theta",  type="float", dest='theta', help='crossover rate', default=0)
    (options, args)=parser.parse_args()

    if options.maternaltbf == None:
        #pass
        sys.stderr.write("please provide maternal 2bit file!\n")
        exit(1)
    if options.paternaltbf == None:
        #pass
        sys.stderr.write("please provide paternal 2bit file!\n")
        exit(1)


    paternalname = str.split ( os.path.basename(options.paternaltbf) , '.')[0]
    maternalname = str.split ( os.path.basename(options.maternaltbf) , '.')[0]





    try:
        sys.stderr.write("opening paternal twobitfile...\n")
        twobit_paternal=bx.seq.twobit.TwoBitFile( open( options.paternaltbf ) )
    except:
        sys.stderr.write("unable to open paternal twobit file!\n")


    try:
        sys.stderr.write("opening maternal twobitfile...\n")
        twobit_maternal=bx.seq.twobit.TwoBitFile( open( options.maternaltbf ) )
    except:
        sys.stderr.write("unable to open maternal twobit file!\n")


    crossoverfh=open("sample."+options.samplename+".meiosis_"+str(options.nmeiosis)+".crossover.log", 'w')

    chromSizes=[]
    
    if options.theta == 0:
        lambdaval= float(options.chromSize)/float(options.morgan)
    else:
        lambdaval=options.theta


    sys.stderr.write("generating crossover intervals ...\n")
    
    gametesequence=[]

    chrom='1'
    
    crossoverfh.write("lambda value: " + str(lambdaval) + "\n")
    #get the number of crossovers by drawing from poisson
    n=np.random.poisson(lambdaval, 1)
    ncrossovers=n[0]
    #pick haplotype to pass
    flip=np.random.uniform(0,1)
    if flip <= .5: 
        haplotype="maternal"
    else:
        haplotype="paternal"
    intervalstring=''
    #now get the crossover points and record them in Interval format http://107.20.183.198/documentation.php#interval-format
    #if zero crossovers, then just past the entire chromosome
    if ncrossovers ==0:
        intervalstring="\t".join([haplotype, chrom, "+",'0', str(options.chromSize), str(ncrossovers),',',','])
        crossoverfh.write(chrom + " " + haplotype + " " + str(ncrossovers) + " " + str(options.chromSize) + "\n")
    #otherise randomly throw down crossover points according to uniform random
    else:
        chromsize=options.chromSize
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
    if haplotype=='maternal':
        inbtwn_name='paternal'
    else:
        inbtwn_name='maternal'
    print intervalstring
        
    

    #extract DNA sequence from the maternal/paternal 2bit files according to the crossover intervals
    for (bedstring) in intervalobj.toBedstring():
        crossoverfh.write(bedstring+"\n")
        if haplotype=='maternal':
            (chr, start, end, name)=bedstring.split('\t')
            chr=maternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_maternal)
        else:
            (chr, start, end, name)=bedstring.split('\t')
            chr=paternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_paternal)
    #print "=="
    for bedstring in intervalobj.yieldInBtwnBedstring(inbtwn_name):
        crossoverfh.write(bedstring+"\n")
        if inbtwn_name=='maternal':
            (chr, start, end, name)=bedstring.split('\t')
            chr=maternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_maternal)
        else:
            (chr, start, end, name)=bedstring.split('\t')
            chr=paternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_paternal)
    if intervalobj.getLastSubIntervalEnd() != -1 and intervalobj.getLastSubIntervalEnd() < chromsize:
        #print inbtwn_name, chrom,  intervalobj.getLastSubIntervalEnd(), chromsize
        bedstring="\t".join(map( lambda(x): str(x), [ chrom, intervalobj.getLastSubIntervalEnd(), chromsize, inbtwn_name ] ) )
        crossoverfh.write(bedstring+"\n")
        if inbtwn_name=='maternal':
            (chr, start, end, name)=bedstring.split('\t')
            chr=maternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_maternal)
        else:
            (chr, start, end, name)=bedstring.split('\t')
            chr=paternalname
            bedstring="\t".join([ chr, start,end,name])
            gametesequence += twoBitExtract(bedstring, twobit_paternal)
    #finally write out teh fasta file
    sys.stderr.write("writing gametic fasta file ... " + chrom + "\n")
    gametefastaname="_".join([chrom,"gamete", "meiosis"+str(options.nmeiosis), options.samplename])
    writefasta("".join(gametesequence), chrom + "." + options.samplename, gametefastaname+".fa")


if __name__ == "__main__":
    main()
