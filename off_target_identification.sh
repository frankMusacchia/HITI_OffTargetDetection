
###This script allows th accurate identification of off-target from CRISPr Analysis
#using BWA, PICARD, SAMTOOLS, AWK 
##
#Copyright Francesco Musacchia © 2021
#
#######CONFIGURATION PARAMETERS

#Software
bin_fold="/home/users/f.musacchia/bin/"
#Script that is able to launch jobs with PBS
exec_job_path="/data/bioinfo/fmusacchia/VarGeniusBeta/LIB/exec_job.pl"
#Use your java path
java_path="/opt/software/java/jdk1.8.0_31/jre/bin/java"
#Path to PICARD
picard_path="/home/users/musacchia/bin/picard.jar"
#Folder where the samples are given 
work_fold="/data/bioinfo/fmusacchia/Analisi/AuricchioLab/"
#temporary folder
tmp_dir="/home/users/f.musacchia/tmp"

#Foilder with fastq files
fastq_fold="/mnt/volume/input/AuricchioLab/AAV_CRISPR_201909/"

#Once the fastq uniq files are ready, this is the analysis part:
#mm10_genome_fa="/cineca/prod/biodata/igenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa"
mm10_genome_fa="/home/ngsworkspace/references/ucsc.mm10.fa"
aav_genome_fa="genome_aav.fasta"
#PARAMETERS
#alb or rho: in insert_f use lowercase, in suff use uppercase
insert_f="alb_insert.fasta"
suff="ALB"
ucsc_gene="ucsc_albumin.fasta"
#RHO
#insert_f="rho_insert.fasta"
#suff="RHO"
#ucsc_gene="ucsc_rhodopsine.fasta"

bwaprog="mem"
#bwaparams=" -w 1 -O 12 -B 8 -A 2 -d 200"
#Very stringent
bwaparams=" -O 60 -B 50 -E 10 -L 100 -c 2 -U 50 -d 200 -w 1"

####################### FASTQ SETTINGS
#Go to working folder, so that all log files are saved ther
cd /data/bioinfo/fmusacchia/Analisi/AuricchioLab/

##Tutti i fastq deveono essere in formato GZIP
#I first had to merge together the fastq files from different lanes with zcat
task="zcat"
for sample in AA_1 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8 AA_9;do
command="zcat $fastq_fold/${sample}/${sample}_L00*_R1* > $fastq_fold/${sample}/${sample}_R1.fq"
echo "Executing $command"
qsub  -r y -N ${sample}_${task}  \
 -v COMMAND="$command" -l walltime=24:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
command="zcat $fastq_fold/${sample}/${sample}_L00*_R2* > $fastq_fold/${sample}/${sample}_R2.fq"
echo "Executing $command"
qsub  -r y -N ${sample}_${task}  \
 -v COMMAND="$command" -l walltime=24:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;


#Then I used gzip to compress
task="gzip"
for sample in AA_1 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8 AA_9;do
command="gzip ${fastq_fold}/${sample}/${sample}_R1.fq"
echo "Executing $command"
qsub  -r y -N ${sample}_${task}  \
 -v COMMAND="$command" -l walltime=24:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
command="gzip ${fastq_fold}/${sample}/${sample}_R2.fq"
echo "Executing $command"

############ ############ ############ 
############ ANALYSIS 1: Obtain all the sequences containing the insert
############ ############ ############ 


#STEP1: Align the fastq against the insert (Sequence obtained from BenchLing)
task="bwamem"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/applications/bwa/0.7.15/intel--cs-xe-2015--binary/bin/bwa mem ${bwaparams} \
FRANCESCO/PROVE/AURICCHIO/${insert_f} \
${fastq_fold}/${sample}/${sample}_R1.fq.gz \
${fastq_fold}/${sample}/${sample}_R2.fq.gz > ${fastq_fold}/${sample}/${sample}_${suff}.sam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-o FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.log \
-e FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.err -v COMMAND="$command" -l walltime=24:00:00 \
-l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
echo "Executng $command"
done;
 
STEP 2 Filter only those regions overlapping the insert and not completely covered (AWK COMMAND)
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
grep 'insert' ${fastq_fold}/${sample}/${sample}_${suff}.sam | awk -F"\t" '$6 ~ /S/' > ${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_external.sam
done;

##STEP 3: Remove the insert and get a fasta. Use a perl script to get the fasta with the external mapped regions
#The fasta_utils.pl script additionally writes also a fasta file [SAMENAME].match.fa which contains the 
#part of the read that did not match and all the sequences will bring the read name to be reassociated later 
#The script produces external_externalOnly.fasta and external_externalOnly.fasta.match.fa
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
perl ${bin_fold}/fasta_utils.pl 8 ${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_external.sam 30
done;

##DA QUI IN CASO DI REANALISI POST Giugno2020 e da AGGIUSTARE - A maggio 2020 ho cambiato lo snip di codice precedente
#STEP2 ALIGN AGAINST THE GENE THE READS CLEANED OF THE INSERT
#Now I have to index the gene from UCSC
/cineca/prod/applications/bwa/0.7.15/intel--cs-xe-2015--binary/bin/bwa index FRANCESCO/PROVE/AURICCHIO/${ucsc_gene}


#The mapped fasta are now realigned against the GENE
task="bwamem"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/applications/bwa/0.7.15/intel--cs-xe-2015--binary/bin/bwa mem ${bwaparams} FRANCESCO/PROVE/AURICCHIO/${ucsc_gene} \
${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_externalOnly.fasta > ${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_externalOnly.sam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-o FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.log \
-e FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.err \
-v COMMAND="$command" -l walltime=2:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

Convert to BAM and sort 
task="sort"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/compilers/jre/1.8.0_73/none/bin/java -Xmx4g  -jar /cineca/prod/applications/picard/2.3.0/binary/bin/picard.jar SortSam \
TMP_DIR=/pico/scratch/userexternal/fmusacch  SORT_ORDER=coordinate \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_externalOnly.sam \
OUTPUT=${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_externalOnly_sort.bam"
qsub  -r y -N ${sample}_${task}_${suff}  -o FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.log \
-e FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.err -v COMMAND="$command" -l walltime=2:00:00 \
-l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

#INDEX
task="index"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/compilers/jre/1.8.0_73/none/bin/java -Xmx4g  -jar /cineca/prod/applications/picard/2.3.0/binary/bin/picard.jar BuildBamIndex  \
TMP_DIR=/pico/scratch/userexternal/fmusacch  \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_externalOnly_sort.bam"
qsub  -r y -N ${sample}_${task}_${suff}  -o FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.log \
-e FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.err -v COMMAND="$command" -l walltime=2:00:00 \
-l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

#RESULT: These BAM files obtained can be viewed in IGV to check how many times the insert has been detected in RHO and ALB
######################






############ ############ ############ 
#####################Analysis 2. Obtain all sequences containing the insert on the genome
#ALIGN AGAINST THE GENOME THE READS CLEANED OF THE INSERT
############ ############ ############ 


#We decided to remove short fasta sequences (<30nt) which could be aligned in different locus of the genome

#The mapped fasta are now realigned against the entire Mouse Genome
task="aln"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${bin_fold}/bwa-0.7.17/bwa mem ${bwaparams} ${mm10_genome_fa} \
${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_external_externalOnly.fasta > ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly.sam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} -v COMMAND="$command" -l walltime=2:00:00 -l select=1:ncpus=1:mem=20GB ${exec_job_path}
done;


#Count alignments on different chromosomes
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
grep '@' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly.sam > ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam
grep '\schr' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly.sam | grep 'MD:Z:[0-9]\+\s' |awk -F"\t" '$5 >30' | awk -F"\t" '$6 !~ /I|S|H/' >> ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  | cut -f3 | sort  | uniq -c > ${fastq_fold}/${sample}/${sample}_${suff}_WG_alnstats
done;

#Convert to BAM to see the alignments
#Convert to BAM and sort

task="sort"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${java_path} -Xmx4g  -jar ${picard_path} SortSam \
TMP_DIR=${tmp_dir} SORT_ORDER=coordinate \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam \
OUTPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch_sort.bam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} \
-v COMMAND="$command" -l walltime=2:00:00 -l select=1:ncpus=1:mem=20GB ${exec_job_path}
done;

#INDEX
task="index"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${java_path} -Xmx4g  -jar ${picard_path} BuildBamIndex  \
TMP_DIR=${tmp_dir}  \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch_sort.bam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} \
-v COMMAND="$command" -l nodes=1 ${exec_job_path}
done;


############ ############ ############
###################ANALYSIS 3  - Align the reads containing the CUT insert against the AAV
############ ############ ############ 


# ALIGN AGAINST THE GENE THE READS CLEANED OF THE INSERT
#Now I have to index the gene from UCSC
${bin_fold}/bwa-0.7.17/bwa index ${work_fold}/${aav_genome_fa}

#Getting a fasta file containing only the insert sequences corresponding with those fitered before
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
#Modify the fasta file to be grep-ed later: sequence and seq name on the same line 
sed -r ':a;$!{N;ba};s/((seq)[^\n]*)\n/\1 /g' ${fastq_fold}/${sample}/${sample}_${suff}_bothEndsMapped_external_externalOnly.fasta.match.fa > ${fastq_fold}/${sample}/temp.fa
#Get the sequence name for each exact match filtered
grep '^seq' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam | cut -f1 | cut -f2 -d'-'| sort > ${fastq_fold}/${sample}/filtseq.temp.txt
#Get the part of the insert for each sequence which was filtered
grep -f ${fastq_fold}/${sample}/filtseq.temp.txt ${fastq_fold}/${sample}/temp.fa | tr ' ' '\n' > ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.match.fa
done;

##For Double check, align also the part of reads that matched with the insert with the AAV genome
task="aln"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${bin_fold}/bwa-0.7.17/bwa mem ${bwaparams} ${work_fold}/${aav_genome_fa} \
${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.match.fa > ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_AAVmatch.sam"
qsub  -r y -N ${sample}_${task}_${suff} -v COMMAND="$command" -l nodes=1 ${exec_job_path}
done;

#Convert to BAM to see the alignments
#Convert to BAM and sort
task="sort"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${java_path} -Xmx4g  -jar ${picard_path} SortSam \
TMP_DIR=${tmp_dir}  SORT_ORDER=coordinate \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_AAVmatch.sam \
OUTPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_AAVmatch_sort.bam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} \
-v COMMAND="$command" -l nodes=1 ${exec_job_path}
done;

#INDEX
task="index"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${java_path} -Xmx4g  -jar ${picard_path} BuildBamIndex  \
TMP_DIR=${tmp_dir}  \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_AAVmatch_sort.bam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} \
-v COMMAND="$command" -l nodes=1 ${exec_job_path}
done;

####################




########### ############ ############
#ANALYSIS 4 -  ALIGN AGAINST THE GENE ALL THE READS TO FIND ALL THE READS WHICH HAVE ALIGNED AROUND THE "CUT"
########### ############ ############


#CUT regions identified using Benchling by Manel
#ALB
start=750
end=850
cut_reg="750-850"

#RHO
start="15"
end="165"
cut_reg="${start}_${end}"

#1. Alignment against the gene
task="bwamem"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="${bin_fold}/bwa-0.7.17/bwa mem ${bwaparams} \
FRANCESCO/PROVE/AURICCHIO/${ucsc_gene} \
${fastq_fold}/${sample}/${sample}_R1.fq.gz ${fastq_fold}/${sample}/${sample}_R2.fq.gz > ${fastq_fold}/${sample}/${sample}_${suff}_only.sam"
qsub -q ngs -r y -N ${sample}_${task}_${suff} \
-v COMMAND="$command" -l nodes=1 ${exec_job_path}
done;

#2. Get only sequences which aligned on the gene
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
grep 'mm10_dna' ${fastq_fold}/${sample}/${sample}_${suff}_only.sam > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.sam
done;

#3. use the perl script to get the fasta with the sequences with any type of alignment against the gene
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
perl FRANCESCO/PROVE/AURICCHIO/sam2fasta.pl \
${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.sam \
${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.fasta ALL
done;

#4. Align against the insert to get sequences with the insert or not
task="bwamem"
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/applications/bwa/0.7.15/intel--cs-xe-2015--binary/bin/bwa mem  ${bwaparams} \
FRANCESCO/PROVE/AURICCHIO/${insert_f} \
${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.fasta > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT.sam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-o FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.log -e FRANCESCO/PROVE/AURICCHIO/${sample}_${task}_${suff}.err \
-v COMMAND="$command" -l walltime=2:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

#5. Reads nel sito del taglio 
#a. Take only those reads which have aligned on the insert, get only those :
	#1. aligned on the first X nucleotides of the insert (for ALB I used $4<10 for RHO $4<20)
	#2-matching something external to the insert (only left), get only the read name, sort
#b. grep the extracted reads from the SAM file aligned against the gene and get the positions list
#c. get start and end (first and last element of previous out)
#d. get  from the SAM only reads starting after the start or before the end
#e. re-get all the reads from the SAM file
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
grep 'insert' ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT.sam | awk -F"\t" '$4<10' | awk -F"\t" '$6 ~ /S/' | awk -F"\t" '$6 !~ /S$/' | cut -f1  | sort | cut -f2 -d'-' > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT_cutreads.txt
grep -f ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT_cutreads.txt ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.sam  | cut -f4 | sort -k1,1n | uniq > ${fastq_fold}/${sample}/${sample}_${suff}_startend.txt
start=`head -n1 ${fastq_fold}/${sample}/${sample}_${suff}_startend.txt`
end=`tail -n1 ${fastq_fold}/${sample}/${sample}_${suff}_startend.txt`
awk -v s="$start" -v e="$end" -F"\t" '$4>s && $4<e ' ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.sam  | cut -f1 | sort | uniq  > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_readsin_${start}-${end}.txt
done;

#MANUAL COMMAND
#In some cases (no insert alignment) there is no start and end. HEnce to count the number of alignment reads in the CUT region I used the start and end from Benchlig
awk -v s="$start" -v e="$end" -F"\t" '$4>s && $4<e ' AA_11/AA_11_ALB_onlyALB.sam | wc -l


#Final STATS
#fter that to count the alignment on the insert, you need to grep 'insert' on the final SAM
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
grep 'insert' ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT.sam > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT_onlyINSERT.sam
done;

task="grep"
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
start=`head -n1 ${fastq_fold}/${sample}/${sample}_${suff}_startend.txt`
end=`tail -n1 ${fastq_fold}/${sample}/${sample}_${suff}_startend.txt`
command="grep -f ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_readsin_${start}-${end}.txt ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT_onlyINSERT.sam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-o ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_INSERT_onlyINSERT_${start}-${end}.sam \
-v COMMAND="$command" -l walltime=24:00:00 -l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;





####For Rhodopsine I have only checked the number of reads aligned in the region of the CUT which, from Benchling, 
#with Manel we saw is at the nucleotide 90 of the gene. 
#Hence after the alignment against the RHO gene I have only to convert to BAM and get "samtools view" a region 
#of 150 nt with 90 being the center: 15-165
#Starts from the STEP 2 of previous TASK
#cut_reg="15-165"
#3. Convert to BAM and sort
task="sort"
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/compilers/jre/1.8.0_73/none/bin/java -Xmx4g  -jar /cineca/prod/applications/picard/2.3.0/binary/bin/picard.jar SortSam \
TMP_DIR=/pico/scratch/userexternal/fmusacch  SORT_ORDER=coordinate \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_only${suff}.sam \
OUTPUT=${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_sort.bam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-v COMMAND="$command" -l walltime=2:00:00 \
-l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

#4: Index
task="index"
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
command="/cineca/prod/compilers/jre/1.8.0_73/none/bin/java -Xmx4g  -jar /cineca/prod/applications/picard/2.3.0/binary/bin/picard.jar BuildBamIndex  \
TMP_DIR=/pico/scratch/userexternal/fmusacch  \
INPUT=${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_sort.bam"
qsub  -r y -N ${sample}_${task}_${suff}  \
-v COMMAND="$command" -l walltime=2:00:00 \
-l select=1:ncpus=1:mem=20GB VarGeniusBeta/LIB/exec_job.pl
done;

#5. Now get only those reads belongin to the CUT region (which I saw in IGV)
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
samtools view ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_sort.bam "RHO_mm10_dna:15-165" > ${fastq_fold}/${sample}/${sample}_${suff}_only${suff}_${cut_reg}.sam
done;

#####################################
## STATISTICS

Count the number of reads inside the Fastq
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8 AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for f in ${fastq_fold}/${sample}/*fq.gz; do
echo $f;zcat $f | awk '{s++}END{printf "%.2f\n", s/4}';
done > ${fastq_fold}/${sample}/${sample}_reads_num.txt
done;

###############APRIL 2020
Manel asked me to count in ALB the total count of reads that map with the dsRED

Count only those regions overlapping the insert and not completely covered (AWK COMMAND)
suff="ALB"
fastq_fold="/gss/gss_work/DRES_TIGEM/ngsworkspace/input/AuricchioLab/AAV_CRISPR_201906/"
out_f="ALB_dsReads_mapped_all.txt"
echo "sample\toverlapped\tfull\n" > $out_f
for sample in AA_9 AA_10 AA_11 AA_12 AA_13 AA_14 AA_15 AA_16;do
for sample in AA_1 AA_2 AA_3 AA_4 AA_5 AA_6 AA_7 AA_8;do
echo $sample >> $out_f
grep 'insert' ${fastq_fold}/${sample}/${sample}_${suff}.sam | awk -F"\t" '$6 ~ /S/' | wc -l >> $out_f
grep 'insert' ${fastq_fold}/${sample}/${sample}_${suff}.sam | grep '151M' | wc -l >> $out_f
done;


#####################################
#On may 2020 Manel asked me to find for ALBumin the count of all the off-target regions
#scramble 11,12,15
#dsred 9,10,13,14,16

#1. Take only chromosome and position from the SAM file and repeat position column to generate a BED file
#2. Sort for mergeBed
#3. Apply mergeBed considering maximum distance among intervals 100nt and later cutting the frequencies minor than 100 as well

#First remove any existing file because we append
suff="ALB"
rm ${suff}_dsRed_mapOfftarget.bed 
rm ${suff}_scramble_mapOfftarget.bed 
rm ${fastq_fold}/${suff}_dsRed_WG_InsertAlnstats.txt
rm ${fastq_fold}/${suff}_scramble_WG_InsertAlnstats.txt
rm ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt

echo "sample sequencing allchr ALBsite"
for sample in AA_9 AA_10 AA_13 AA_14 AA_16;do
cut -f3,4 ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam | grep '^chr' |  awk -F"\t" ' { print $1"\t"$2"\t"$2 }' >> ${suff}_dsRed_mapOfftarget.bed
#Total reads across all samples
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  | cut -f3 | sort  | uniq -c >> ${fastq_fold}/${suff}_dsRed_WG_InsertAlnstats.txt

Count reads per sample without chr5 and only chr5
echo -n "${sample} gRNA " >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  |  wc -l >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  | grep 'chr5' | wc -l >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
done;

sort -k1,1 -k2,2n ${suff}_dsRed_mapOfftarget.bed > ${suff}_dsRed_mapOfftarget_sort.bed
${bin_fold}/bedtools2/bin/mergeBed -i ${suff}_dsRed_mapOfftarget_sort.bed -c 1 -o count -d 100 | awk -F"\t" '$4 >100  { print $1"\t"$2"\t"$3"\t"$4 }' | sort -k4,4nr > ${suff}_dsRed_mapOfftarget_sort_stats.txt.

#Scramble samples
for sample in AA_11 AA_12 AA_15;do
cut -f3,4 ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam | grep '^chr' |  awk -F"\t" ' { print $1"\t"$2"\t"$2 }' >> ${suff}_scramble_mapOfftarget.bed
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  | cut -f3 | sort  | uniq -c  >> ${fastq_fold}/${suff}_scramble_WG_InsertAlnstats.txt

#Count reads per sample without chr5 and only chr5
echo -n "${sample} scramble " >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  |  wc -l >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
awk -F"\t" '$3 ~ /chr/' ${fastq_fold}/${sample}/${sample}_${suff}_WG_bothEndsMapped_externalOnly_exactMatch.sam  | grep 'chr5' | wc -l >> ${fastq_fold}/${suff}_WG_InsertAlnCountPerSample.txt
done;
sort -k1,1 -k2,2n ${suff}_scramble_mapOfftarget.bed > ${suff}_scramble_mapOfftarget_sort.bed
${bin_fold}/bedtools2/bin/mergeBed -i ${suff}_scramble_mapOfftarget_sort.bed -c 1 -o count -d 100 | awk -F"\t" '$4 >100  { print $1"\t"$2"\t"$3"\t"$4 }' | sort -k4,4nr > ${suff}_scramble_mapOfftarget_sort_stats.txt

#Reorganize the files and count the number of reads per-chromosome for either dsRed and scramble
sed 's/ \+ //' ${fastq_fold}/${suff}_dsRed_WG_InsertAlnstats.txt | sort -k2,2 | awk '{a[$2]+=$1}END{for(i in a) print i,a[i]}' > ${fastq_fold}/${suff}_dsRed_WG_InsertAlnstats_all.txt
sed 's/ \+ //' ${fastq_fold}/${suff}_scramble_WG_InsertAlnstats.txt | sort -k2,2 | awk '{a[$2]+=$1}END{for(i in a) print i,a[i]}' > ${fastq_fold}/${suff}_scramble_WG_InsertAlnstats_all.txt






