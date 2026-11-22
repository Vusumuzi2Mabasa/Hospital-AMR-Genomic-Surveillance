"""Assembly QC, strict seven-locus allele calls and reference-projected consensus."""
import argparse,csv,json,re
from pathlib import Path

def fasta(path):
 records={};name=None
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];records[name]=''
  elif name is not None:records[name]+=line.strip().upper()
 return records
def write_fasta(path,name,seq):Path(path).write_text('>'+name+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
def reference(genome):
 records=fasta(genome);name=max(records,key=lambda k:len(records[k]));seq=records[name]
 assert len(seq)>4000000,'Expected a complete chromosome reference'
 write_fasta('reference.fasta','reference_chromosome',seq)
 Path('reference.json').write_text(json.dumps({'source_contig':name,'length':len(seq),'excluded_accessory_contigs':len(records)-1},indent=2)+'\n')
def analyse(sample,genome,sam,blast,profiles):
 records=fasta(genome);lengths=sorted(map(len,records.values()),reverse=True);total=sum(lengths);n50=next(lengths[i] for i in range(len(lengths)) if sum(lengths[:i+1])>=total/2)
 bases=''.join(records.values());qc={'sample':sample,'assembly_bases':total,'contigs':len(records),'contigs_ge500':sum(n>=500 for n in lengths),'N50':n50,'GC_fraction':(bases.count('G')+bases.count('C'))/len(bases),'ambiguous_fraction':sum(b not in 'ACGT' for b in bases)/len(bases)}
 # Exact full-length allele matches only. Multiple distinct alleles remain ambiguous.
 calls={l:set() for l in ('gapA','infB','mdh','pgi','phoE','rpoB','tonB')}
 for line in Path(blast).read_text().splitlines():
  q,s,ident,length,qlen,qstart,qend,sstart,send=line.split('\t');locus,allele=q.rsplit('_',1)
  if locus in calls and float(ident)==100 and int(length)==int(qlen) and int(qstart)==1 and int(qend)==int(qlen):calls[locus].add(allele)
 alleles={k:next(iter(v)) if len(v)==1 else ('ambiguous' if v else 'missing') for k,v in calls.items()}
 st=[]
 for row in csv.DictReader(Path(profiles).open(),delimiter='\t'):
  if all(row[k]==v for k,v in alleles.items()):st.append(row['ST'])
 qc.update({'sequence_type':';'.join(st) if st else 'unassigned',**{'allele_'+k:v for k,v in alleles.items()}})
 ref=next(iter(fasta('reference.fasta').values()));n=len(ref);out=bytearray(b'N'*n);seen=bytearray(n);indel=bytearray(n);mapped=0;aligned=0;different=0
 for line in Path(sam).open():
  if line.startswith('@'):continue
  f=line.rstrip().split('\t');flag=int(f[1])
  if flag&(4|256|2048) or int(f[4])<30 or f[2]!='reference_chromosome' or f[9]=='*':continue
  mapped+=1;r=int(f[3])-1;q=0;seq=f[9].encode()
  for count,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
   count=int(count)
   if op in ('M','=','X'):
    for offset in range(count):
     pos=r+offset;letter=seq[q+offset]
     if pos>=n:raise ValueError('Alignment exceeds reference')
     seen[pos]=min(255,seen[pos]+1)
     if seen[pos]==1 and letter in b'ACGT':out[pos]=letter
     else:out[pos]=ord('N')
     aligned+=1;different+=letter!=ord(ref[pos])
    r+=count;q+=count
   elif op in ('I','D','N'):
    start=max(0,r-5);end=min(n,r+(count if op in ('D','N') else 0)+5)
    indel[start:end]=b'\1'*(end-start)
    if op=='I':q+=count
    else:r+=count
   elif op=='S':q+=count
 for i in range(n):
  if indel[i] or seen[i]!=1:out[i]=ord('N')
 qc.update({'reference_bases':n,'mapped_primary_contigs':mapped,'callable_bases_before_shared_mask':sum(b in b'ACGT' for b in out),'multiple_alignment_bases':sum(k>1 for k in seen),'indel_masked_bases':sum(indel),'aligned_identity_approx':1-different/aligned if aligned else None})
 qc['callable_fraction_before_shared_mask']=qc['callable_bases_before_shared_mask']/n
 write_fasta(sample+'.consensus.fasta',sample,out.decode())
 Path(sample+'.qc.json').write_text(json.dumps(qc,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['reference','sample']);p.add_argument('--sample');p.add_argument('--genome',required=True);p.add_argument('--sam');p.add_argument('--blast');p.add_argument('--profiles');a=p.parse_args()
 if a.mode=='reference':reference(a.genome)
 else:analyse(a.sample,a.genome,a.sam,a.blast,a.profiles)
