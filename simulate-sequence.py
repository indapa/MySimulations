#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
import numpy as np

""" simulate a sequence based on equal base proportions of nucleotides """

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--Afreq", type="float", dest="A", default=.25)
    parser.add_option("--Cfreq", type="float", dest="C", default=.25)
    parser.add_option("--Gfreq", type="float", dest="G", default=.25)
    parser.add_option("--Tfreq", type="float", dest="T", default=.25)
    parser.add_option("--seqname", type="string", dest="seqname", default="simref")
    parser.add_option("--length", type="int", dest="length", default=1000000)

    (options, args)=parser.parse_args()
    sequence=''


    print ">"+options.seqname

    for i in range(0, options.length):
        prob=np.random.uniform(0,1)
        base='N'
        if prob < 0.25:
            base='A'
        if prob >=.25 and prob < .50:
            base='C'
        if prob >=.50 and prob < .75:
            base='G'
        if prob >=.75:
            base='T'

        sequence += base

  

    i=0

    while i < options.length:
        subseq=sequence[i:i+60]
        print subseq
        i+=60

if __name__ == "__main__":
    main()
