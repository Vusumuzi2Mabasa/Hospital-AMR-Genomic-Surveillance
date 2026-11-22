import sys
sample=sys.argv[1]
with open(sample+'.mask.bed','w') as out:
    for line in open(sample+'.depth.tsv'):
        chrom,pos,depth=line.split()
        if int(depth)<5:out.write(f'{chrom}\t{int(pos)-1}\t{pos}\n')
