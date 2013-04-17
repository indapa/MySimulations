#!/usr/bin/env python
import sys, os
from optparse import OptionParser
import bx.seq.twobit

""" construct a *.ped file with genotypes from haplotype 2bit files and 
    segregating sites *.map file generated/based on coalescent simulation from  cosi """

def main():

    usage="usage: %prog [options] reference.2bit"
    parser=OptionParser(usage)
    paternalTbfile='/Users/indapa/data/MySimulations/Simulation1/ChildOne/1_gamete_meiosis1_paternalgamete.2bit'
    parser.add_option("--hap1", type="string", dest="hapOne", help="2bit file haplotype one")
    parser.add_option("--hap2", type="string", dest="hapTwo", help="2bit file haplotype two")
    parser.add_option("--map", type="string", dest="mapfile", help="map file of polymorphic positions")
    parser.add_option("--sample", type="string", dest='samplename', help="samplename of individual")
    parser.add_option("--sex", type="string", dest='sex', help="sex of individual [1/2]")
    
    (options,args)=parser.parse_args()
    if not args:
        sys.stderr.write(usage+"\n")
        sys.exit(1)
    
    referencefile=args[0]
    
    hapOneTbf=bx.seq.twobit.TwoBitFile ( open(options.hapOne) )
    hapTwoTbf=bx.seq.twobit.TwoBitFile ( open(options.hapTwo) )
    referenceTbf=bx.seq.twobit.TwoBitFile( open(referencefile) )


    #get a list of segregating sites
    mapfh=open(options.mapfile,'r')
    chrom='simref.1'
    hapOneChrom=os.path.splitext(os.path.basename(options.hapOne))[0]
    hapTwoChrom=os.path.splitext(os.path.basename(options.hapTwo))[0]
    

    genotypes=[]
    pedstring="\t".join(['fam1', options.samplename, '0', '0', options.sex, '-9'])

    for line in mapfh:
        #print line.strip()
        (chr,start,end,ref,alt)=line.strip().split('\t')
        #print chr, start,end,ref,alt
        start=int(start)
        end=int(end)
        sequence=referenceTbf[chrom][start:end]
        hapOne_seq=hapOneTbf[hapOneChrom][start:end]
        hapTwo_seq=hapTwoTbf[hapTwoChrom][start:end]
        #print sequence, paternal_seq, maternal_seq
    
        hapone_genotype=''
        haptwo_genotype=''
    
        if hapOne_seq == sequence:
            hapone_genotype='2'
        else:
            hapone_genotype='1'
        
        if hapTwo_seq== sequence:
            haptwo_genotype='2'
        else:
            haptwo_genotype='1'
    
        diplotype=' '.join([hapone_genotype, haptwo_genotype])
    
        genotypes.append( diplotype )

    genotype_string=" ".join(genotypes)
    print "\t".join( [pedstring, genotype_string])

if __name__ == "__main__":
    main()