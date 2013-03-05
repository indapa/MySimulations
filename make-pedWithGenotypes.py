import itertools

file='/Users/amit/software/cosi_1.2/examples/bestfit/out.hap-1'
#http://stackoverflow.com/a/1657385/1735942
FatherPedstring="\t".join(['fam1', 'FatherOne', '0', '0', '1', '-9'])
MotherPedstring="\t".join(['fam1', 'MotherOne', '0', '0', '2', '-9'])
""" we read lines from the haplotype output file from cosi 2 at a time """
with open(file) as fh:
    for line1,line2 in itertools.izip_longest(*[fh]*2):
        (hapid1, popid1, genostr_one)=line1.strip().split('\t')
        (hapid2, popid2, genostr_two)=line2.strip().split('\t')

        if int(hapid1) == 0:
            pedstring=FatherPedstring
        else:
            pedstring=MotherPedstring
        if int(hapid2) > 3:
            break
         
        
        genotypes_one=genostr_one.split(' ')
        genotypes_two=genostr_two.split(' ')
        
        genotypes=[]
        for hap1,hap2 in itertools.izip(genotypes_one, genotypes_two):
            diplotype=' '.join([hap1,hap2])
            genotypes.append(diplotype)
        """ we collect the two haploytpes into genotypes """
        genotypes='\t'.join(genotypes)
        """ we print out the per-formatted file with genotypes: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped """
        print "\t".join( [pedstring, genotypes])
