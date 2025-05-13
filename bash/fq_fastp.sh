files=($(find fq -maxdepth 1 -type f -printf "%f\n"))

for file in ${files[@]}
do
	name=${file%%.*}
	echo ${name}
	fastp -i fq/${name}.fastq.gz -o trim_fq/${name}.trim.fastq.gz 
done
