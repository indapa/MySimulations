#!/usr/bin/env python
import sys, collections, argparse
import csv
""" a utility script to generate bamtools merge given a seqindex for pgmsnp test data """
def main():
    usage="usage: %prog [options] seqindex.txt"
    parser = argparse.ArgumentParser(description='generate bamtools merge script')
    parser.add_argument('-seqindex', metavar='seqindex', type=str,help='seqindex.txt')
    parser.add_argument('-scriptname', metavar='scriptname', type=str, help='name scriptname to write bamtools merge commands to',default='bam-merge.sh')
    parser.add_argument('-bamtoolspath', metavar='bamtoolspath', type=str, help='path to bamtools executable', default='/Users/indapa/software/bamtools/bin/bamtools')
    parser.add_argument('-mergedname', metavar='mergedname', type=str, help='name of merged BAM file')
    parser.add_argument('rm', metavar='sample', type=str, nargs='*',help='samplenames to remove from merge',default=[])
    args = parser.parse_args()

    seqindexfh=open(args.seqindex)
    

    #seqindexfh.readline()

    bamtoolsIn=[]
    record = collections.namedtuple('record', 'fid sample fullpath md5sum coverage')

    for line in seqindexfh:
        fields=line.strip().split('\t')
        b=record(*fields)
        if b.sample not in args.rm:
            bamtoolsIn.append(b.fullpath)


    bamtoolsIn=[ "-in " + rec for rec in bamtoolsIn]

            
    inputs=" ".join(bamtoolsIn)
    outfh=open(args.scriptname, 'w')
    commandline= args.bamtoolspath + " sort " +  inputs + " -out " + args.mergedname
    indexcommandline=args.bamtoolspath + " index " + "-in " + args.mergedname
    outfh.write(commandline+"\n")
    outfh.write(indexcommandline+"\n")

    

if __name__ == "__main__":
    main()
