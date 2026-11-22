import unittest,tempfile,subprocess,sys,gzip,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class ReadsTest(unittest.TestCase):
    def test_truth_and_controls(self):
        import json
        t=json.loads((ROOT/'data/wgs_truth.json').read_text())
        self.assertEqual(len(t),12)
        self.assertEqual(sum(x['assembly'] for x in t.values()),4)
        self.assertEqual(sum(x['depth']==2 for x in t.values()),2)
    def test_fastq_structure_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as folder:
            cmd=[sys.executable,str(ROOT/'scripts/simulate_reads.py'),'--fasta',str(ROOT/'data/wgs_truth/WGS01.fasta'),'--sample','WGS01','--depth','2']
            subprocess.run(cmd,cwd=folder,check=True,capture_output=True)
            p=Path(folder)/'WGS01_R1.fastq.gz'
            first=gzip.decompress(p.read_bytes())
            subprocess.run(cmd,cwd=folder,check=True,capture_output=True)
            self.assertEqual(first,gzip.decompress(p.read_bytes()))
            lines=first.decode().splitlines();self.assertEqual(len(lines)%4,0)
            for i in range(0,len(lines),4):
                self.assertTrue(lines[i].startswith('@'));self.assertEqual(len(lines[i+1]),150);self.assertEqual(len(lines[i+3]),150)
if __name__=='__main__':unittest.main()
