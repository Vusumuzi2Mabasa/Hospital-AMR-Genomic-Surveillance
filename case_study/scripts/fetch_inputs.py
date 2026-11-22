"""Download accession-locked public assemblies and allele references; verify SHA256."""
import argparse,csv,hashlib,io,json,time,zipfile
from pathlib import Path
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'case_study'
def sha(data):return hashlib.sha256(data).hexdigest()
def get(url):
 for attempt in range(4):
  try:return urlopen(Request(url,headers={'User-Agent':'Hospital-AMR-public-reanalysis/1.0'}),timeout=120).read()
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--metadata-snapshot');parser.add_argument('--amr',action='store_true');args=parser.parse_args()
 cfg=json.loads((CASE/'config/study.json').read_text())
 cache=CASE/'cache';cache.mkdir(exist_ok=True);(cache/'genomes').mkdir(exist_ok=True)
 lockpath=CASE/'provenance/input_manifest.json'
 if lockpath.exists():lock=json.loads(lockpath.read_text())
 else:
  if not args.metadata_snapshot:raise SystemExit('Initial curation requires --metadata-snapshot; routine reproduction uses the committed manifest.')
  reports=json.loads(Path(args.metadata_snapshot).read_text())['reports'];records=[];chosen=[]
  for sample in cfg['samples']:
   candidates=[r for r in reports if r['accession'].startswith('GCA_') and r['organism'].get('infraspecific_names',{}).get('strain')==sample]
   assert len(candidates)==1,(sample,len(candidates))
   r=candidates[0];chosen.append(r);a=r['accession'];info=r['assembly_info'];bio=info.get('biosample',{})
   records.append({'sample':sample,'accession':a,'biosample':bio.get('accession',''),'bioproject':info.get('bioproject_accession',''),'release_date':info.get('release_date',''),'assembly_level':info['assembly_level'],'assembly_name':info['assembly_name'],'sequencing_technology':info.get('sequencing_tech',''),'archive_species':r['organism']['organism_name'],'ncbi_best_ani_species':r.get('average_nucleotide_identity',{}).get('best_ani_match',{}).get('organism_name',''),'url':f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{a}/download?include_annotation_type=GENOME_FASTA'})
  (CASE/'provenance/selected_ncbi_metadata.json').write_text(json.dumps(chosen,indent=2)+'\n')
  lock={'publication':cfg['publication_doi'],'assemblies':records,'references':[],'selection_frozen_before_analysis':True}
 def download_genome(row):
  path=cache/'genomes'/(row['sample']+'.fna')
  if not path.exists():
   data=get(row['url'])
   with zipfile.ZipFile(io.BytesIO(data)) as z:
    files=[n for n in z.namelist() if n.endswith('_genomic.fna')];assert len(files)==1,files
    path.write_bytes(z.read(files[0]))
  digest=sha(path.read_bytes())
  if row.get('sha256'):assert digest==row['sha256'],path
  row['sha256']=digest;return row
 with ThreadPoolExecutor(max_workers=2) as pool:lock['assemblies']=list(pool.map(download_genome,lock['assemblies']))
 names=['gapA.tfa','infB.tfa','mdh.tfa','pgi.tfa','phoE.tfa','rpoB.tfa','tonB.tfa','klebsiella.txt']
 old={r['name']:r for r in lock['references']};refs=[]
 (cache/'mlst').mkdir(exist_ok=True)
 for name in names:
  url=f"https://raw.githubusercontent.com/tseemann/mlst/{cfg['mlst_snapshot']}/db/pubmlst/klebsiella/{name}"
  p=cache/'mlst'/name
  if not p.exists():p.write_bytes(get(url))
  digest=sha(p.read_bytes())
  if name in old:assert digest==old[name]['sha256'],name
  refs.append({'name':name,'url':url,'sha256':digest})
 lock['references']=refs
 (cache/'mlst/alleles.fasta').write_bytes(b'\n'.join((cache/'mlst'/n).read_bytes() for n in names if n.endswith('.tfa')))
 lockpath.write_text(json.dumps(lock,indent=2)+'\n')
 with (CASE/'config/samples.tsv').open('w',newline='') as out:
  fields=['sample','accession','biosample','bioproject','release_date','assembly_level','sequencing_technology','archive_species']
  w=csv.DictWriter(out,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(lock['assemblies'])
 if args.amr:
  manifest=json.loads((CASE/'provenance/amrfinder_database_manifest.json').read_text());dest=cache/'amrfinder'/cfg['amrfinder_database'];dest.mkdir(parents=True,exist_ok=True)
  for item in manifest['files']:
   p=dest/item['name']
   if not p.exists():p.write_bytes(get(item['url']))
   assert sha(p.read_bytes())==item['sha256'],p
 print('Verified',len(lock['assemblies']),'public assemblies and',len(refs),'MLST source files.')
if __name__=='__main__':main()
