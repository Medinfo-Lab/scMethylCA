files=($(find . -maxdepth 1 -type d -printf "%f\n"))

for file in ${files[@]}
do
	echo ${file}
	gunzip -c ${file}/${file}_1.trim_bismark_bt2_pe.deduplicated.bismark.cov.gz >${file}/${file}_R1.trim_bismark_bt2_pe.deduplicated.bismark.cov
	awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' ${file}/${file}_R1.trim_bismark_bt2_pe.deduplicated.bismark.cov >${file}/${file}_R1.trim_bismark_bt2_pe.deduplicated.bismark.chr1_22_XY.cov
done


files=($(find methy -maxdepth 1 -type d -printf "%f\n"))

for file in ${files[@]}
do
	echo ${file}
	mkdir cov/${file}
	coverage2cytosine --nome-seq --gc --genome_folder ~/genomes/hg19_bismark -o ${file} --dir cov/${file} methy/${file}/${file}.trim_bismark_bt2.deduplicated.bismark.chr1_22_XY.cov
done