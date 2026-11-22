"""Build non-biological sequence fixtures and a two-gene public reference subset.

Random core segments are not bacterial genomes. Real public AMR sequences serve
only as computational positive controls; no species, plasmid or phenotype is inferred.
"""
import argparse,hashlib,json,random,csv
from pathlib import Path

def fasta(path):
    records={};header=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):header=line[1:];records[header]=''
        elif header is not None:records[header]+=line.strip().upper()
    return records

def write_fasta(path,records):
    with Path(path).open('w',newline='\n') as out:
        for header,seq in records.items():
            out.write('>'+header+'\n');out.write('\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--ncbi-cds');args=parser.parse_args()
    target=Path('data/reference/amr_targets.fasta')
    if args.ncbi_cds:
        selected={};meta=[]
        for header,seq in fasta(args.ncbi_cds).items():
            fields=header.split('|')
            for gene in ('blaKPC-2','blaNDM-1'):
                symbols=(gene,gene[3:])
                if any(s in fields for s in symbols) and gene not in selected:
                    selected[gene]=seq
                    meta.append({'gene':gene,'source_header':header,'length':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest()})
        assert set(selected)=={'blaKPC-2','blaNDM-1'},'Both exact reference symbols are required'
        write_fasta(target,selected)
        Path('provenance/amr_reference.json').write_text(json.dumps({'database':'NCBI AMRFinderPlus Reference Gene Database','release':'2026-08-07.1','source':'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/2026-08-07.1/AMR_CDS.fa','source_sha256':hashlib.sha256(Path(args.ncbi_cds).read_bytes()).hexdigest(),'records':meta,'purpose':'Two-target BLAST positive controls, not a full AMRFinderPlus analysis'},indent=2)+'\n')
    genes=fasta(target);assert len(genes)==2
    rng=random.Random(20261122);L=30000
    ancestor=''.join(rng.choices('ACGT',k=L))
    def mutate(seq,n):
        seq=list(seq)
        for j in rng.sample(range(len(seq)),n):seq[j]=rng.choice([b for b in 'ACGT' if b!=seq[j]])
        return ''.join(seq)
    founders=[mutate(ancestor,60),mutate(ancestor,60)]
    root=Path('data/wgs_truth');root.mkdir(parents=True,exist_ok=True)
    write_fasta('data/reference/synthetic_core.fasta',{'synthetic_core':ancestor})
    rows=[];truth={};core={}
    for i in range(1,13):
        sid=f'WGS{i:02d}'
        if i<=4:g=0;label='Cluster A'
        elif i<=8:g=1;label='Cluster B'
        elif i<=10:g=i-9;label='Close unlinked background'
        else:g=None;label='Distant background'
        seq=mutate(founders[g],i%3+1) if g is not None else mutate(ancestor,100)
        gene='blaKPC-2' if i in (1,2,3,9) else 'blaNDM-1' if i in (5,6,7,10) else ''
        records={'synthetic_core':seq}
        if gene:records['synthetic_accessory_control']=''.join(rng.choices('ACGT',k=1000))+genes[gene]+''.join(rng.choices('ACGT',k=1000))
        write_fasta(root/(sid+'.fasta'),records);core[sid]=seq;truth[sid]={'cluster':label,'gene':gene or 'None','depth':28 if i<=10 else 2,'assembly':i in (1,4,5,8)}
        rows.append({'sample':sid,'depth':truth[sid]['depth'],'assembly':int(truth[sid]['assembly'])})
    write_fasta('data/wgs_true_alignment.fasta',core)
    with Path('data/wgs_samples.tsv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=['sample','depth','assembly'],delimiter='\t');w.writeheader();w.writerows(rows)
    Path('data/wgs_truth.json').write_text(json.dumps(truth,indent=2)+'\n')
    print('Prepared 12 random-core fixtures, four assembly cases and two public AMR target references.')

if __name__=='__main__':main()
