import csv,importlib.util,json,os,tempfile,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('assembly_analysis',Path(__file__).with_name('assembly_analysis.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class ConsensusTests(unittest.TestCase):
 def run_fixture(self,cigar='32=',flag=0,mapq=60,duplicate=False):
  old=os.getcwd()
  with tempfile.TemporaryDirectory() as tmp:
   try:
    os.chdir(tmp);seq='ACGT'*8;m.write_fasta('reference.fasta','reference_chromosome',seq);m.write_fasta('sample.fna','contig',seq)
    row=f'q\t{flag}\treference_chromosome\t1\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\n';Path('sample.sam').write_text(row*(2 if duplicate else 1))
    loci=['gapA','infB','mdh','pgi','phoE','rpoB','tonB'];Path('alleles.tsv').write_text('\n'.join(f'{l}_1\tc\t100\t32\t32\t1\t32\t1\t32' for l in loci))
    Path('profiles.tsv').write_text('ST\t'+'\t'.join(loci)+'\n258\t'+'\t'.join(['1']*7)+'\n')
    m.analyse('X','sample.fna','sample.sam','alleles.tsv','profiles.tsv')
    return json.loads(Path('X.qc.json').read_text()),next(iter(m.fasta('X.consensus.fasta').values()))
   finally:os.chdir(old)
 def test_complete_exact_sequence_and_typing(self):
  q,s=self.run_fixture();self.assertEqual(s,'ACGT'*8);self.assertEqual(q['sequence_type'],'258');self.assertEqual(q['N50'],32)
 def test_low_map_quality_abstains(self):
  q,s=self.run_fixture(mapq=0);self.assertEqual(s,'N'*32);self.assertEqual(q['mapped_primary_contigs'],0)
 def test_secondary_alignment_excluded(self):
  q,s=self.run_fixture(flag=256);self.assertEqual(s,'N'*32)
 def test_overlapping_primary_contigs_are_masked(self):
  q,s=self.run_fixture(duplicate=True);self.assertEqual(s,'N'*32);self.assertEqual(q['multiple_alignment_bases'],32)
 def test_deletion_and_flanking_bases_are_masked(self):
  q,s=self.run_fixture(cigar='10=2D20=');self.assertEqual(s[5:17],'N'*12);self.assertGreater(q['indel_masked_bases'],0)
if __name__=='__main__':unittest.main()
