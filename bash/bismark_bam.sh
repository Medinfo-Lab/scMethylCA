files=($(find trim_fq -maxdepth 1 -type f -name "*.trim.fastq.gz" -printf "%f\n"))

for file in ${files[@]}
do
	name=${file%%.*}
	echo ${name}
	mkdir bam/${name}
	bismark --bowtie2 --parallel 10 --genome /home/zhoushitao/genomes/hg19_bismark --bam -o bam/${name} trim_fq/${name}.trim.fastq.gz 
done


files=($(find bam -maxdepth 1 -type d -printf "%f\n"))

for file in ${files[@]}
do
	echo ${file}
	mkdir deduplicate_bam/${file}
	deduplicate_bismark --bam bam/${file}/${file}.trim_bismark_bt2.bam --output_dir deduplicate_bam/${file}
done


files=($(find deduplicate_bam -maxdepth 1 -type d -printf "%f\n"))

for file in ${files[@]}
do
	echo ${file}
	mkdir methy/${file}
	bismark_methylation_extractor --genome_folder ~/genomes/hg19_bismark --bedGraph -s --comprehensive --buffer_size 10G --multicore 15 -o methy/${file} deduplicate_bam/${file}/${file}.trim_bismark_bt2.deduplicated.bam
done