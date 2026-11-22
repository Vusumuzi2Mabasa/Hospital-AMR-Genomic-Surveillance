nextflow.enable.dsl=2
params.outdir = 'results/wgs'

process READ_FIXTURES {
    tag { sample }
    input:
    tuple val(sample), val(depth), path(genome), val(assemble)
    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), val(assemble)
    script:
    """
    python3 ${projectDir}/scripts/simulate_reads.py --fasta ${genome} --sample ${sample} --depth ${depth}
    """
}
process FASTQ_QC {
    tag { sample }
    publishDir "${params.outdir}/fastp", mode:'copy', pattern:'*.json'
    input:
    tuple val(sample), path(r1), path(r2), val(assemble)
    output:
    tuple val(sample), path("${sample}.R1.fq.gz"), path("${sample}.R2.fq.gz"), val(assemble), emit: reads
    path "${sample}.fastp.json", emit: qc
    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample}.R1.fq.gz -O ${sample}.R2.fq.gz \
      --disable_adapter_trimming --length_required 100 --qualified_quality_phred 20 \
      --dont_eval_duplication --thread 2 --json ${sample}.fastp.json --html ${sample}.fastp.html
    """
}
process MAP_AND_CALL {
    tag { sample }
    publishDir "${params.outdir}/consensus", mode:'copy', pattern:'*.core.fasta'
    input:
    tuple val(sample), path(r1), path(r2), val(assemble)
    path reference
    output:
    path "${sample}.core.fasta", emit: core
    script:
    """
    cp ${reference} ref.fa
    bwa index ref.fa
    bwa mem -t 2 ref.fa ${r1} ${r2} | samtools sort -@ 1 -o ${sample}.bam
    samtools index ${sample}.bam
    samtools depth -aa -q 20 -Q 20 ${sample}.bam > ${sample}.depth.tsv
    python3 ${projectDir}/scripts/mask_low_depth.py ${sample}
    bcftools mpileup -Ou -f ref.fa -q 20 -Q 20 -a FORMAT/DP ${sample}.bam | \
      bcftools call -m --ploidy 1 -Oz -o raw.vcf.gz
    bcftools filter -i 'QUAL>=20 && FORMAT/DP>=5' raw.vcf.gz -Oz -o filtered.vcf.gz
    bcftools index filtered.vcf.gz
    bcftools consensus -f ref.fa -m ${sample}.mask.bed filtered.vcf.gz | \
      sed 's/^>.*/>${sample}/' > ${sample}.core.fasta
    """
}
process ASSEMBLE_AND_SCREEN {
    tag { sample }
    publishDir "${params.outdir}/assembly", mode:'copy', pattern:'*.blast.tsv'
    input:
    tuple val(sample), path(r1), path(r2), val(assemble)
    path amr
    output:
    path "${sample}.blast.tsv", emit: hits
    script:
    """
    spades.py --only-assembler -k 21,33,55 -1 ${r1} -2 ${r2} -t 2 -m 2 -o assembly
    blastn -query assembly/contigs.fasta -subject ${amr} \
      -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' > ${sample}.blast.tsv
    """
}
process AUDIT_SEQUENCES {
    publishDir params.outdir, mode:'copy'
    input:
    path cores
    path qc
    path hits
    path truth
    path expected
    output:
    path 'qc_passing_alignment.fasta', emit: alignment
    path '*.csv'
    path 'wgs_validation.json'
    script:
    """
    python3 ${projectDir}/scripts/evaluate_wgs.py
    """
}
process MAXIMUM_LIKELIHOOD_TREE {
    publishDir "${params.outdir}/phylogeny", mode:'copy'
    input:
    path alignment
    output:
    path 'observed.*'
    script:
    """
    iqtree2 -s ${alignment} -m GTR+G -B 1000 -T 2 --seed 20261122 --prefix observed
    """
}
workflow {
    samples=Channel.fromPath("${projectDir}/data/wgs_samples.tsv").splitCsv(header:true,sep:'\t').map{r->tuple(r.sample,r.depth.toInteger(),file("${projectDir}/data/wgs_truth/${r.sample}.fasta"),r.assembly=='1')}
    READ_FIXTURES(samples)
    FASTQ_QC(READ_FIXTURES.out)
    MAP_AND_CALL(FASTQ_QC.out.reads,Channel.value(file("${projectDir}/data/reference/synthetic_core.fasta")))
    ASSEMBLE_AND_SCREEN(FASTQ_QC.out.reads.filter{s,r1,r2,a->a},Channel.value(file("${projectDir}/data/reference/amr_targets.fasta")))
    AUDIT_SEQUENCES(MAP_AND_CALL.out.core.collect(),FASTQ_QC.out.qc.collect(),ASSEMBLE_AND_SCREEN.out.hits.collect(),Channel.value(file("${projectDir}/data/wgs_truth.json")),Channel.value(file("${projectDir}/data/wgs_true_alignment.fasta")))
    MAXIMUM_LIKELIHOOD_TREE(AUDIT_SEQUENCES.out.alignment)
}
