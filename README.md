# Phoebe-zhennan-genome

./OrthoFinder/orthofinder -f proteinS/longestpep/ -S diamond -M msa -A mafft -t 8 -a 8





hisat2-build -p 20 reference.fa reference
hisat2 -x reference -p 20 --dta --summary-file SRRxxxxxxx_log -U SRRxxxxxxx.fastq.gz | samtools view -@ 20 -bh | samtools sort -@ 20 -o SRRxxxxxxx_sorted.bam

stringtie -o SRRxxxxxxx.gtf  -p 20 -A SRRxxxxxxx.tab -G reference.gff3 -B -e SRRxxxxxxx_sorted.bam
stringtie --merge -e -p 10 -G reference.gff3  -o stringtie_merged.gtf mergelist
stringtie -e -B -p 8 -G stringtie_merged.gtf -o SRRxxxxxxx_merged.gtf SRRxxxxxxx_sorted.bam

python /data/01/user119/tools/py/prepDE.py -i sample_list -l 150
