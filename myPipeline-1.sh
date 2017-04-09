#!/bin/sh
#PBS -N myPipeline-1
#PBS -r n
#PBS -e $PBS_JOBNAME.$PBS_JOBID.e
#PBS -o $PBS_JOBNAME.$PBS_JOBID.o
#PBS -m ae
#PBS -M fiyhigher@163.com
#PBS -l nodes=1:ppn=8,mem=30GB
#PBS -V

# add path of the scripts used in the pipeline
source /public/biodata/wei-lab-data/scripts/pipeline_bashrc.sh
source /public/users/liangminling/liangqx/bash_profile

# Have script report non-zero exit status if any command fails
set -e
set -o pipefail

##########################
SAMPLE='ERR011103'       
##########################

SRA='/public/users/liangminling/liangqx/ncbi/public/sra'
HOME='/public/users/liangminling/liangqx'

cd $HOME
mkdir ${SAMPLE}_genome

# split SRA paired-end 
cd $SRA
fastq-dump --split-files $SAMPLE.sra
mv ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq $HOME/${SAMPLE}_genome

cd $HOME/${SAMPLE}_genome
statfile=$SAMPLE.txt
echo -n "" > $statfile
echo "total_sequences" >> $statfile
wc -l ${SAMPLE}_1.fastq | awk '{print $1/4}' >> $statfile

#############################################################################################

# adapter, quality, length control
module load trim_galore
trim_galore --phred33 --length 50 --paired ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
echo "after_trimm_galore" >> $statfile
wc -l ${SAMPLE}_1_val_1.fq | awk '{print $1/4}' >> $statfile

if [ -f ${SAMPLE}_1_unpaired_1.fq ]; then rm ${SAMPLE}_1_unpaired_1.fq; fi;
if [ -f ${SAMPLE}_2_unpaired_2.fq ]; then rm ${SAMPLE}_2_unpaired_2.fq; fi;

rm ${SAMPLE}_1.fastq
rm ${SAMPLE}_2.fastq

#############################################################################################

# trim to same length (75bp), remove ambiguous reads, remove duplicated reads
prinseq-lite.pl -trim_to_len 75 -ns_max_p 10 -derep 12345 -fastq ${SAMPLE}_1_val_1.fq -fastq2 ${SAMPLE}_2_val_2.fq -out_bad null -out_good $SAMPLE.prinseq

# Get the number of sequences after removing  ambiguous reads
echo "after_prinseq" >> $statfile
wc -l $SAMPLE.prinseq_1.fastq | awk '{print $1/4}' >> $statfile

if [ -f $SAMPLE.prinseq_1_singletons.fastq ]; then rm $SAMPLE.prinseq_1_singletons.fastq; fi;
if [ -f $SAMPLE.prinseq_2_singletons.fastq ]; then rm $SAMPLE.prinseq_2_singletons.fastq; fi;

rm ${SAMPLE}_1_val_1.fq
rm ${SAMPLE}_2_val_2.fq

#############################################################################################

# remove human by Hisat2
/public/users/dengyuhua/software/hisat2-2.0.1-beta/hisat2 --dta-cufflinks --phred33 -p 8 -q -x /public/users/dengyuhua/hisat2.index/hg19/genome -1 $SAMPLE.prinseq_1.fastq -2 $SAMPLE.prinseq_2.fastq -S $SAMPLE.hisat.sam
module load samtools
# convert sam to bam
samtools view -bS $SAMPLE.hisat.sam > $SAMPLE.hisat.bam
# extract unmapped reads (non_host)
samtools view -b -f 12 $SAMPLE.hisat.bam > $SAMPLE.hisat.unmapped.bam
# convert bam to fastq
module load bedtools
bamToFastq -i $SAMPLE.hisat.unmapped.bam -fq $SAMPLE.hisat.unmapped_1.fastq -fq2 $SAMPLE.hisat.unmapped_2.fastq

# Get the number of paired sequences unmapped to the host
echo "host_unmapped_reads" >> $statfile
wc -l $SAMPLE.hisat.unmapped_1.fastq | awk '{print $1/4}' >> $statfile

rm $SAMPLE.prinseq_1.fastq
rm $SAMPLE.prinseq_2.fastq
rm $SAMPLE.hisat.bam
rm $SAMPLE.hisat.sam
rm $SAMPLE.hisat.unmapped.bam

#############################################################################################

# Further filter sequences for host genome contamination with BMTagger
bmtagger.sh -C /public/biodata/wei-lab-data/software/bmtools/bmtagger.conf -X -T /tmp -q1 -1 $SAMPLE.hisat.unmapped_1.fastq -2 $SAMPLE.hisat.unmapped_2.fastq -o $SAMPLE.bmt

# Get the number of paired reads passing bmtagger filtering
echo "after_bmtagger" >> $statfile
wc -l $SAMPLE.bmt_1.fastq | awk '{print $1/4}' >> $statfile

rm $SAMPLE.hisat.unmapped_1.fastq
rm $SAMPLE.hisat.unmapped_2.fastq

#############################################################################################

# convert fastq to fasta and filter low complexity reads
seqtk seq -A $SAMPLE.bmt_1.fastq > $SAMPLE.bmt_1.tmp.fasta
seqtk seq -A $SAMPLE.bmt_2.fastq > $SAMPLE.bmt_2.tmp.fasta
# run dustmasker to mask the low complexity sequences
dustmasker -in $SAMPLE.bmt_1.tmp.fasta -out $SAMPLE.bmt_1.dustout
dustmasker -in $SAMPLE.bmt_2.tmp.fasta -out $SAMPLE.bmt_2.dustout
filter_dustmasker_paired.py -1 $SAMPLE.bmt_1.dustout -2 $SAMPLE.bmt_2.dustout > $SAMPLE.dusted.keep_ids_1 2> $SAMPLE.dusted.keep_ids_2
extract_fullseq $SAMPLE.dusted.keep_ids_1 -keep -fastq -single $SAMPLE.bmt_1.fastq > $SAMPLE.dusted_1.fastq
extract_fullseq $SAMPLE.dusted.keep_ids_2 -keep -fastq -single $SAMPLE.bmt_2.fastq > $SAMPLE.dusted_2.fastq

# Get the number of paired reads after filtering for low complexity sequences
echo "after_low_complexity_filtering" >> $statfile
wc -l $SAMPLE.dusted_1.fastq | awk '{print $1/4}' >> $statfile

rm $SAMPLE.bmt_1.fastq
rm $SAMPLE.bmt_2.fastq
rm $SAMPLE.bmt_1.tmp.fasta
rm $SAMPLE.bmt_2.tmp.fasta
rm $SAMPLE.bmt_1.dustout
rm $SAMPLE.bmt_2.dustout
rm $SAMPLE.dusted.keep_ids_1
rm $SAMPLE.dusted.keep_ids_2

#############################################################################################

# align reads to the hmp reference
module load bwa
bwa mem -t 8 -k 20 /public/biodata/wei-lab-data/public_reference/bwa.index/hmp/hmp_genbank $SAMPLE.dusted_1.fastq $SAMPLE.dusted_2.fastq > $SAMPLE.mem.sam
module load samtools
samtools view -bS $SAMPLE.mem.sam > $SAMPLE.mem.bam
samtools view -b -F 12 -q 30 $SAMPLE.mem.bam > $SAMPLE.mem.mapped.bam
echo "mem_mapped_seqs" >> $statfile
samtools flagstat $SAMPLE.mem.mapped.bam | grep 'with itself and mate mapped' | cut -f 1 -d ' ' >> $statfile
samtools sort $SAMPLE.mem.mapped.bam $SAMPLE.mem.mapped.sorted
samtools index $SAMPLE.mem.mapped.sorted.bam
samtools idxstats $SAMPLE.mem.mapped.sorted.bam > $SAMPLE.idxstats.txt
python /public/biodata/wei-lab-data/scripts/idxstat_genome_summary.py -n /public/biodata/wei-lab-data/tax_data/nodes.dmp -m /public/biodata/wei-lab-data/tax_data/names.dmp -k a -i $SAMPLE.idxstats.txt > $SAMPLE.all.txt
python /public/biodata/wei-lab-data/scripts/idxstat_genome_summary.py -n /public/biodata/wei-lab-data/tax_data/nodes.dmp -m /public/biodata/wei-lab-data/tax_data/names.dmp -k b -i $SAMPLE.idxstats.txt > $SAMPLE.bacteria.txt
python /public/biodata/wei-lab-data/scripts/idxstat_genome_summary.py -n /public/biodata/wei-lab-data/tax_data/nodes.dmp -m /public/biodata/wei-lab-data/tax_data/names.dmp -k e -i $SAMPLE.idxstats.txt > $SAMPLE.eukaryota.txt
python /public/biodata/wei-lab-data/scripts/idxstat_genome_summary.py -n /public/biodata/wei-lab-data/tax_data/nodes.dmp -m /public/biodata/wei-lab-data/tax_data/names.dmp -k v -i $SAMPLE.idxstats.txt > $SAMPLE.virus.txt

rm $SAMPLE.mem.sam
rm $SAMPLE.mem.mapped.sorted.bam.bai
gzip $SAMPLE.dusted_1.fastq $SAMPLE.dusted_2.fastq

mv $HOME/*cluster.local* $HOME/standard_output
mv $statfile /public/users/liangminling/liangqx/statfile
mv $SAMPLE.all.txt /public/users/liangminling/liangqx/genome_summary/All
mv $SAMPLE.bacteria.txt /public/users/liangminling/liangqx/genome_summary/Bacteria
mv $SAMPLE.eukaryota.txt /public/users/liangminling/liangqx/genome_summary/Eukaryota
mv $SAMPLE.virus.txt /public/users/liangminling/liangqx/genome_summary/Virus