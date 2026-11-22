nextflow.enable.dsl=2
params.outdir='case_study/results'
params.amr_db='case_study/cache/amrfinder/2026-08-07.1'
process REFERENCE_MASK {
 input:
 path genome
 path code
 output: path 'reference.fasta',emit: reference
 path 'self_repeats.tsv',emit: repeats
 path 'reference.json',emit: metadata
 script:
 """
 python3 ${code} reference --genome ${genome}
 blastn -query reference.fasta -subject reference.fasta -word_size 28 -dust no \
  -perc_identity 95 -outfmt '6 pident length qstart qend sstart send' > self_repeats.tsv
 """
}
process ASSEMBLY_CHARACTERISATION {
 tag { sample }
 publishDir "${params.outdir}/per_sample",mode:'copy',pattern:'*.json'
 publishDir "${params.outdir}/per_sample",mode:'copy',pattern:'*.amr.tsv'
 input:
 tuple val(sample),path(genome)
 path reference
 path alleles
 path profiles
 path amrdb
 path code
 output:
 path '*.consensus.fasta',emit: consensus
 path '*.qc.json',emit: qc
 path '*.amr.tsv',emit: amr
 script:
 """
 minimap2 -a --eqx -x asm5 -t 2 ${reference} ${genome} > ${sample}.sam
 blastn -query ${alleles} -subject ${genome} -perc_identity 100 \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send' > ${sample}.mlst.tsv
 amrfinder -n ${genome} -d ${amrdb} --threads 2 -o ${sample}.amr.tsv
 python3 ${code} sample --sample ${sample} --genome ${genome} \
  --sam ${sample}.sam --blast ${sample}.mlst.tsv --profiles ${profiles}
 """
}
process COMMON_CORE_AND_VALIDATION {
 publishDir params.outdir,mode:'copy',pattern:'*.csv'
 publishDir params.outdir,mode:'copy',pattern:'*.json'
 publishDir params.outdir,mode:'copy',pattern:'core_snps.fasta'
 input:
 path consensus
 path qc
 path amr
 path reference
 path repeats
 path config
 path cohort_code
 path assembly_code
 output:
 path 'core_full.fasta',emit: alignment
 path '*.csv'
 path '*.json'
 path 'core_snps.fasta'
 script:
 """
 python3 ${cohort_code}
 """
}
process PUBLIC_MAXIMUM_LIKELIHOOD_TREE {
 publishDir "${params.outdir}/phylogeny",mode:'copy',pattern:'*.treefile'
 publishDir "${params.outdir}/phylogeny",mode:'copy',pattern:'*.iqtree'
 publishDir "${params.outdir}/phylogeny",mode:'copy',pattern:'*.log'
 input: path alignment
 output: path 'public.*'
 script:
 """
 iqtree2 -s ${alignment} -m GTR+G -B 1000 -T 2 -keep-ident --seed 20261122 --prefix public
 """
}
workflow {
 samples=Channel.fromPath("${projectDir}/config/samples.tsv").splitCsv(header:true,sep:'\t').map{r->tuple(r.sample,file("${projectDir}/cache/genomes/${r.sample}.fna"))}
 assembly_code=Channel.value(file("${projectDir}/scripts/assembly_analysis.py"))
 REFERENCE_MASK(Channel.value(file("${projectDir}/cache/genomes/KPNIH1.fna")),assembly_code)
 ASSEMBLY_CHARACTERISATION(samples,REFERENCE_MASK.out.reference,Channel.value(file("${projectDir}/cache/mlst/alleles.fasta")),Channel.value(file("${projectDir}/cache/mlst/klebsiella.txt")),Channel.value(file(params.amr_db)),assembly_code)
 COMMON_CORE_AND_VALIDATION(ASSEMBLY_CHARACTERISATION.out.consensus.collect(),ASSEMBLY_CHARACTERISATION.out.qc.collect(),ASSEMBLY_CHARACTERISATION.out.amr.collect(),REFERENCE_MASK.out.reference,REFERENCE_MASK.out.repeats,Channel.value(file("${projectDir}/config/study.json")),Channel.value(file("${projectDir}/scripts/cohort_analysis.py")),assembly_code)
 PUBLIC_MAXIMUM_LIKELIHOOD_TREE(COMMON_CORE_AND_VALIDATION.out.alignment)
}
