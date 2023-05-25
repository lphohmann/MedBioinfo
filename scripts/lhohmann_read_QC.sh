#!/bin/bash
echo "script start: download and initial sequencing read quality control"
date

# main: quintessential commands

# 2. Programmatic access to NCBI SRA
# produce a simple text file usually with 4x2 lines (one DNA and one RNA sequencing run per patient), with one run accession number per line
sqlite3 -batch -noheader -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='lhohmann';" > ./lhohmann_run_accessions.txt
# mk directory
mkdir ../data/sra_fastq
# load module 
module load sra-tools
# split the reads in three files (two for paired end forward and reverse reads, one for the singleton reads if any) produce gzip-ed compressed files
cat lhohmann_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip \
--outdir ../data/sra_fastq/ --disable-multithreading --split-e
# check job details
#sacct --format=JobID,JobName%20,ReqCPUS,ReqMem,Timelimit,State,ExitCode,Start,elapsed,MaxRSS,NodeList,Account%15

# 3. Manipulating raw sequencing FASTQ files
# load module
module load seqkit
# stats for each fasta file
srun --cpus-per-task=2 --time=00:30:00 seqkit --threads 2 stats ../data/sra_fastq/*.gz
# comapre to metadata
sqlite3 -batch /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db "select * from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='lhohmann';"

# 4. Quality control
# make directory
mkdir ./fastqc
# load module
module load fastqc
# run fastqc
srun --cpus-per-task=2 --time=00:30:00 xargs -I{} -a lhohmann_run_accessions.txt fastqc --outdir ./fastqc/ \
--threads 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz

# (6. Adaptor and low quality read trimming)
# 7. Merging paired end reads
module load flash2
# create directory
mkdir ../data/merged_pairs
# merging all paired end FASTQ files
srun --cpus-per-task=2 --time=00:30:00 xargs -a lhohmann_run_accessions.txt -n 1 -I{} flash2 --threads=2 -z \
--output-directory=../data/merged_pairs/ --output-prefix={}.flash ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a lhohmann_flash2.log

# 8. Use read mapping to check for PhiX contamination
# download the PhiX genome and use bowtie2 map your FASTQ reads
# create ../data/reference_seqs subdirectory 
mkdir ../data/reference_seqs
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna
# load module
module load bowtie2
# create a bowtie2 indexed database from the reference sequences
mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB
# create directory in analyses/
mkdir ./bowtie
# aligns all your FASTQ merged reads against the index DB
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/lhohmann_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/lhohmann_bowtie_merged2PhiX.log

# which samples may contain traces of SARS-CoV-2 infection
# repeat the whole procedure to align merged reads to the reference SAR-CoV-2 genome (accession number NC_045512)
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/lhohmann_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/lhohmann_bowtie_merged2SC2.log

# 9. Combine quality control results into one unique report for all samples analysed
# load the multiqc module 
module load multiqc
srun multiqc --force --title "lhohmann sample sub-set" ../data/merged_pairs/ ./fastqc/ ./lhohmann_flash2.log ./bowtie/

# script ending
date
echo "script end."
