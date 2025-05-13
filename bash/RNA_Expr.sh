files=($(find trim_fq -maxdepth 1 -name "*.trim.fastq.gz" -type f -printf "%f\n"))

for file in ${files[@]}
do
	name=${file%%.*}
	echo ${name}
	bowtie2 -p 10 -x ~/genomes/hg19_bowtie2/hg19 -U trim_fq/${name}.trim.fastq.gz |samtools view -bS -o rna_bam/${name}.bam
	samtools sort rna_bam/${name}.bam -o rna_bam/${name}.sorted.bam
	samtools index rna_bam/${name}.sorted.bam
done

/work/tools/featureCounts/subread-2.0.1-Linux-x86_64/bin/featureCounts -T 15 -p -t exon -g gene_id -a /work/genomes/UCSC/hg19/gencode.v41lift37.annotation.gtf -o SRP295352_counts.txt rna_bam/*.sorted.bam