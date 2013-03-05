#now take the child haplotypes and make ped file with genotypes

import bx.seq.twobit

paternalTbfile='/Users/amit/data/MySimulations/Simulation1/ChildOne/1_gamete_meiosis1_paternalgamete.2bit'
maternalTbfile='/Users/amit/data/MySimulations/Simulation1/ChildOne/1_gamete_meiosis1_maternalgamete.2bit'
referencefile='/Users/amit/data/MySimulations/Simulation1/Reference/simref.1.2bit'

paternalTbf=bx.seq.twobit.TwoBitFile ( open(paternalTbfile) )
maternalTbf=bx.seq.twobit.TwoBitFile ( open(maternalTbfile) )
referenceTbf=bx.seq.twobit.TwoBitFile( open(referencefile) )

#get a list of segregating sites
mapfile='/Users/amit/data/MySimulations/Simulation1/out.pos-1.1.map'
mapfh=open(mapfile,'r')
chrom='simref.1'
maternal_chrom='1.maternalgamete'
paternal_chrom='1.paternalgamete'

genotypes=[]
ChildPedstring="\t".join(['fam1', 'ChildOne', 'FatherOne', 'MotherOne', '2', '-9'])

for line in mapfh:
    #print line.strip()
    (chr,start,end,ref,alt)=line.strip().split('\t')
    #print chr, start,end,ref,alt
    start=int(start)
    end=int(end)
    sequence=referenceTbf[chrom][start:end]
    paternal_seq=paternalTbf[paternal_chrom][start:end]
    maternal_seq=maternalTbf[maternal_chrom][start:end]
    #print sequence, paternal_seq, maternal_seq
    
    paternal_genotype=''
    maternal_genotype=''
    
    if paternal_seq == sequence:
        paternal_genotype='2'
    else:
        paternal_genotype='1'
        
    if maternal_seq== sequence:
        maternal_genotype='2'
    else:
        maternal_genotype='1'
    
    childgenotype=' '.join([paternal_genotype, maternal_genotype])
    
    genotypes.append( childgenotype )

genotype_string=" ".join(genotypes)
print "\t".join( [ChildPedstring, genotype_string])
