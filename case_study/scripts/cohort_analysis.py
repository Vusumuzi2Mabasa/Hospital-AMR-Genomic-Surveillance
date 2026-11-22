"""Common callable chromosome, SNP distances and prespecified publication checks."""
import csv,glob,hashlib,json
from pathlib import Path
from assembly_analysis import fasta,write_fasta

def csvout(path,rows):
 with Path(path).open('w',newline='') as out:
  w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
cfg=json.loads(Path('study.json').read_text());ref=next(iter(fasta('reference.fasta').values()));n=len(ref);mask=bytearray(n)
for line in Path('self_repeats.tsv').read_text().splitlines():
 ident,length,qstart,qend,sstart,send=line.split('\t');qs,qe,ss,se=map(int,(qstart,qend,sstart,send))
 if float(ident)<95 or int(length)<100 or (qs==ss and qe==se):continue
 for start,end in ((qs,qe),(ss,se)):
  lo,hi=sorted((start,end));mask[lo-1:hi]=b'\1'*(hi-lo+1)
qcs=[json.loads(Path(f).read_text()) for f in glob.glob('*.qc.json')];qcs.sort(key=lambda r:cfg['samples'].index(r['sample']))
seqs={sid:next(iter(fasta(sid+'.consensus.fasta').values())) for sid in cfg['samples']}
assert set(seqs)=={r['sample'] for r in qcs}
passing=[r['sample'] for r in qcs if r['callable_fraction_before_shared_mask']>=cfg['minimum_callable_fraction']]
assert len(passing)>=4,'Insufficient QC-passing assemblies for phylogeny'
core=[i for i in range(n) if not mask[i] and all(seqs[s][i] in 'ACGT' for s in passing)]
assert len(core)>n*.8,'Shared core below prespecified 80% of reference chromosome'
snps=[i for i in core if len({seqs[s][i] for s in passing})>1]
assert snps,'No variable sites; phylogeny should not be inferred'
with Path('core_full.fasta').open('w') as out:
 for s in passing:
  seq=''.join(seqs[s][i] for i in core);out.write('>'+s+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
with Path('core_snps.fasta').open('w') as out:
 for s in passing:out.write('>'+s+'\n'+''.join(seqs[s][i] for i in snps)+'\n')
csvout('core_snp_sites.csv',[{'reference_position_1based':i+1,**{s:seqs[s][i] for s in passing}} for i in snps])
pairs=[]
for i,a in enumerate(passing):
 for b in passing[i+1:]:pairs.append({'sample_a':a,'sample_b':b,'common_core_bases':len(core),'snps':sum(seqs[a][p]!=seqs[b][p] for p in snps)})
csvout('core_pairwise_distances.csv',pairs)
amr=[]
for s in cfg['samples']:
 for r in csv.DictReader(Path(s+'.amr.tsv').open(),delimiter='\t'):
  symbol=r.get('Element symbol',r.get('Gene symbol'))
  if symbol is None:raise ValueError('Unrecognised AMRFinder schema; expected Element symbol or Gene symbol')
  amr.append({'sample':s,'normalised_gene_symbol':symbol,**r})
if amr:csvout('amrfinder_all_hits.csv',amr)
validation=[]
for r in qcs:
 s=r['sample'];hits=[a for a in amr if a['sample']==s and a.get('Type')=='AMR'];genes=sorted({a['normalised_gene_symbol'] for a in hits});kpc=[g for g in genes if g.startswith('blaKPC') or g.startswith('KPC')]
 r.update({'QC_pass':s in passing,'KPC_genes':';'.join(kpc),'AMR_genes':';'.join(genes),'AMR_hit_count':len(hits)})
 validation.append({'sample':s,'expected_ST':'258','observed_ST':r['sequence_type'],'ST_concordant':r['sequence_type']=='258','expected_family':'KPC','observed_KPC':';'.join(kpc) or 'not_detected','KPC_concordant':bool(kpc),'scope':'Publication-level lineage and KPC-family expectations; not a transmission validation'})
csvout('assembly_qc_typing.csv',qcs);csvout('publication_validation.csv',validation)
summary={'samples':len(qcs),'QC_passing':len(passing),'reference_chromosome_bases':n,'repeat_masked_bases':sum(mask),'common_callable_bases':len(core),'common_callable_fraction':len(core)/n,'variable_sites':len(snps),'pairwise_min':min(d['snps'] for d in pairs),'pairwise_max':max(d['snps'] for d in pairs),'ST_concordant':sum(v['ST_concordant'] for v in validation),'KPC_concordant':sum(v['KPC_concordant'] for v in validation),'published_full_cohort_variable_sites':41,'published_cohort_size':18,'exact_41_site_reproduction_claimed':False,'method':'Strict common reference-projected chromosome; MQ>=30; duplicates, indel flanks and self-BLAST repeats masked; no recombination inference','passed_computational_checks':True}
Path('case_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
