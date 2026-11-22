"""Independent truth audit of executed Linux sequencing tools."""
import glob,json,csv,hashlib
from pathlib import Path
from prepare_wgs import fasta,write_fasta
truth=json.loads(Path('wgs_truth.json').read_text());expected=fasta('wgs_true_alignment.fasta')
rows=[];alignment={}
for name in sorted(glob.glob('*.core.fasta')):
    sid=Path(name).name.split('.')[0];seq=next(iter(fasta(name).values()));ref=expected[sid]
    assert len(seq)==len(ref)==30000
    valid=[i for i,b in enumerate(seq) if b in 'ACGT'];wrong=sum(seq[i]!=ref[i] for i in valid)
    alignment[sid]=seq
    qc=json.loads(Path(sid+'.fastp.json').read_text());before=qc['summary']['before_filtering'];after=qc['summary']['after_filtering']
    rows.append({'sample':sid,'expected_depth':truth[sid]['depth'],'callable_sites':len(valid),'callable_fraction':len(valid)/len(seq),'wrong_callable_bases':wrong,'input_reads':before['total_reads'],'retained_reads':after['total_reads'],'q30_rate':after['q30_rate']})
assert len(rows)==12
high=[d for d in rows if d['expected_depth']==28];low=[d for d in rows if d['expected_depth']==2]
assert all(d['callable_fraction']>.9 for d in high),'High-depth consensus coverage failure'
assert all(d['wrong_callable_bases']<=2 for d in high),'Unexpected consensus errors'
assert all(d['callable_fraction']<.5 for d in low),'Low-depth negative controls should fail callability'
with open('wgs_qc.csv','w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write_fasta('observed_alignment.fasta',alignment)
# The tree includes QC-passing sequences only; it is not inferred from truth.
passing={d['sample'] for d in rows if d['callable_fraction']>=.9}
write_fasta('qc_passing_alignment.fasta',{sid:alignment[sid] for sid in alignment if sid in passing})
amr=[]
for path in sorted(glob.glob('*.blast.tsv')):
    sid=Path(path).name.split('.')[0];hits=set()
    for line in Path(path).read_text().splitlines():
        q,s,pident,length,qlen,slen,evalue,bitscore=line.split('\t')
        if float(pident)>=95 and int(length)/int(slen)>=.9:hits.add(s)
    want=set() if truth[sid]['gene']=='None' else {truth[sid]['gene']}
    amr.append({'sample':sid,'expected':';'.join(sorted(want)) or 'None','detected':';'.join(sorted(hits)) or 'None','correct':hits==want})
assert len(amr)==4 and all(d['correct'] for d in amr),'Assembly AMR positive/negative controls failed'
with open('assembly_amr_checks.csv','w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(amr[0]));w.writeheader();w.writerows(amr)
pairs=[]
ids=list(alignment)
for i,a in enumerate(ids):
    for b in ids[i+1:]:
        keep=[j for j in range(30000) if alignment[a][j] in 'ACGT' and alignment[b][j] in 'ACGT']
        pairs.append({'sample_a':a,'sample_b':b,'joint_callable':len(keep),'observed_snps':sum(alignment[a][j]!=alignment[b][j] for j in keep),'true_snps_on_callable_sites':sum(expected[a][j]!=expected[b][j] for j in keep),'true_full_core_snps':sum(x!=y for x,y in zip(expected[a],expected[b]))})
with open('observed_pair_distances.csv','w',newline='') as out:
    w=csv.DictWriter(out,fieldnames=list(pairs[0]));w.writeheader();w.writerows(pairs)
Path('wgs_validation.json').write_text(json.dumps({'passed':True,'samples':12,'high_depth':len(high),'low_depth_negative_controls':len(low),'assembly_controls':len(amr),'max_high_depth_wrong_bases':max(d['wrong_callable_bases'] for d in high),'scope':'Random 30kb core segments; four assemblies screened against two AMR references; not clinical genomes or a complete resistome'},indent=2)+'\n')
print('Read QC, mapping, consensus, low-depth and assembly AMR checks passed.')
