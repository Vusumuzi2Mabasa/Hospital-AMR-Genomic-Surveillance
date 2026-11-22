"""Deterministic paired-end FASTQ generator for short computational fixtures."""
import argparse,random,gzip
from prepare_wgs import fasta
p=argparse.ArgumentParser();p.add_argument('--fasta',required=True);p.add_argument('--sample',required=True);p.add_argument('--depth',type=int,required=True);a=p.parse_args()
r=random.Random(20261122+int(a.sample[-2:]));records=fasta(a.fasta);length=150;insert=320
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def errors(s):return ''.join(r.choice([b for b in 'ACGT' if b!=x]) if r.random()<1e-4 else x for x in s)
with gzip.open(a.sample+'_R1.fastq.gz','wt') as f1,gzip.open(a.sample+'_R2.fastq.gz','wt') as f2:
    count=0
    for contig,seq in records.items():
        for j in range(round(a.depth*len(seq)/(2*length))):
            start=r.randrange(0,len(seq)-insert+1);frag=seq[start:start+insert];count+=1
            for f,read,mate in ((f1,errors(frag[:length]),1),(f2,errors(rc(frag[-length:])),2)):
                f.write(f'@{a.sample}:{count}/{mate}\n{read}\n+\n'+('I'*length)+'\n')
print(a.sample,count,'pairs')
