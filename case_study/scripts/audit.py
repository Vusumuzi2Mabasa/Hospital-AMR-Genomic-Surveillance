"""Independent release checks of raw AMR evidence, SNP arithmetic and preservation."""
import argparse,csv,hashlib,json,re
from pathlib import Path
from assembly_analysis import fasta
ROOT=Path(__file__).resolve().parents[2];CASE=ROOT/'case_study'
def read(name):return list(csv.DictReader((CASE/'results'/name).open(newline='')))
def true(v):return str(v).lower()=='true'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--ci',action='store_true');args=parser.parse_args();checks={}
 def check(name,value):
  checks[name]=bool(value)
  if not value:raise AssertionError(name)
 cfg=json.loads((CASE/'config/study.json').read_text());q=read('assembly_qc_typing.csv');v=read('publication_validation.csv');d=read('core_pairwise_distances.csv');s=json.loads((CASE/'results/case_summary.json').read_text());seqs=fasta(CASE/'results/core_snps.fasta')
 check('all_prespecified_samples_reported',set(r['sample'] for r in q)==set(cfg['samples']) and len(q)==12)
 check('all_pairwise_distances_recomputed',len(d)==66 and all(sum(a!=b for a,b in zip(seqs[r['sample_a']],seqs[r['sample_b']]))==int(r['snps']) for r in d))
 check('variant_site_count_recomputed',all(len(x)==s['variable_sites'] for x in seqs.values()) and len(read('core_snp_sites.csv'))==s['variable_sites'])
 check('common_denominator_consistent',all(int(r['common_core_bases'])==s['common_callable_bases'] for r in d))
 check('prespecified_QC_gate',all(true(r['QC_pass'])==(float(r['callable_fraction_before_shared_mask'])>=.9) for r in q))
 check('core_fraction_recomputed',abs(s['common_callable_bases']/s['reference_chromosome_bases']-s['common_callable_fraction'])<1e-12)
 raw_kpc={}
 for sample in cfg['samples']:
  raw=list(csv.DictReader((CASE/'results/per_sample'/(sample+'.amr.tsv')).open(),delimiter='\t'))
  assert raw and ('Element symbol' in raw[0] or 'Gene symbol' in raw[0]),'Unexpected AMRFinder schema'
  raw_kpc[sample]=sorted({r.get('Element symbol',r.get('Gene symbol','')) for r in raw if r.get('Element symbol',r.get('Gene symbol','')).startswith(('blaKPC','KPC'))})
 check('KPC_calls_match_original_AMRFinder_outputs',all(r['KPC_genes']==';'.join(raw_kpc[r['sample']]) for r in q))
 check('publication_concordance_recomputed',all(true(r['ST_concordant'])==(str(r['observed_ST'])=='258') and true(r['KPC_concordant'])==bool(raw_kpc[r['sample']]) for r in v))
 check('summary_concordance_matches_table',sum(true(r['ST_concordant']) for r in v)==s['ST_concordant'] and sum(true(r['KPC_concordant']) for r in v)==s['KPC_concordant'])
 tips=re.findall(r'KPNIH\d+',(CASE/'results/phylogeny/public.treefile').read_text());check('tree_tips_match_QC_samples',len(tips)==12 and set(tips)==set(seqs))
 check('four_case_figures_present',all((CASE/'figures'/(n+'.png')).exists() for n in ['public_phylogeny','public_snp_heatmap','public_amr_heatmap','public_assembly_qc']))
 if not args.ci:
  baseline=json.loads((CASE/'provenance/preserved_simulation_sha256.json').read_text(encoding='utf-8-sig'))
  changed=[r['path'] for r in baseline if digest(ROOT/r['path'])!=r['sha256']]
  check('original_simulation_files_unchanged',not changed)
  locked=json.loads((CASE/'provenance/input_manifest.json').read_text())
  check('downloaded_assemblies_match_locked_SHA256',all(digest(CASE/'cache/genomes'/(r['sample']+'.fna'))==r['sha256'] for r in locked['assemblies']))
  trace=list(csv.DictReader((CASE/'results/nextflow_trace.tsv').open(),delimiter='\t'))
  check('all_15_Nextflow_tasks_successful',len(trace)==15 and all(r['status'] in ('COMPLETED','CACHED') and r['exit']=='0' for r in trace))
 output=CASE/'provenance'/('CI_validation.json' if args.ci else 'case_validation.json')
 output.write_text(json.dumps({'passed':True,'checks':checks,'mode':'CI committed-output audit' if args.ci else 'Local complete case-study and preservation audit','check_count':len(checks),'publication_concordance_is_not_diagnostic_validation':True},indent=2)+'\n')
 print(len(checks),'public case-study release checks passed.')
if __name__=='__main__':main()
