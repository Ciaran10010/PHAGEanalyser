#!/bin/bash

source /home/ciaran/anaconda3/bin/activate main

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo Please input the directory of the raw Illumina reads.
read rd

for phage_directory in "$rd"/*; do
	phage_folder=${phage_directory##*/}
	echo $phage_folder
	mkdir "$phage_directory"/fastqc
	mkdir "$phage_directory"/log
	touch "$phage_directory"/log/initial_fastqc.txt
	echo Performing FastQC analysis of initial reads
	fastqc -quiet -o "$phage_directory"/fastqc "$phage_directory"/*.fastq.gz &>> "$phage_directory"/log/initial_fastqc.txt
	
	read_1=$(ls "$phage_directory"/*1.fastq.gz)
	read_2=$(ls "$phage_directory"/*2.fastq.gz)
	
	touch "$phage_directory"/log/BBtools_initial_clean
	echo Cleaning reads with BBtools
	bbduk.sh in1="$read_1" in2="$read_2" out1="$phage_directory"/nophix_1.fastq.gz out2="$phage_directory"/nophix_2.fastq.gz ref="$script_dir"/BBMap_39.06/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=phiX_stats.txt &>> "$phage_directory"/log/BBtools_initial_clean
	bbduk.sh in1="$phage_directory"/nophix_1.fastq.gz in2="$phage_directory"/nophix_2.fastq.gz out1="$phage_directory"/clean_1.fastq.gz out2="$phage_directory"/clean_2.fastq.gz ref="$script_dir"/BBMap_39.06/bbmap/resources/adapters.fa ktrim=r hdist=1 tpe tbo minlen=100 qtrim=r1 trimq=28 &>> "$phage_directory"/log/BBtools_initial_clean
	
	touch "$phage_directory"/log/clean_fastqc.txt
	echo Performing FastQC analysis of the cleaned reads
	fastqc -quiet -o "$phage_directory"/fastqc "$phage_directory"/clean*.fastq.gz &>> "$phage_directory"/log/clean_fastqc.txt
	
	clean_read_1=$(ls "$phage_directory"/clean_1.fastq.gz)
	clean_read_2=$(ls "$phage_directory"/clean_2.fastq.gz)

	unzip $phage_directory/fastqc/clean_1_fastqc.zip -d $phage_directory/fastqc
	coverage=$( awk 'NR==10 {print $3}' $phage_directory/fastqc/clean_1_fastqc/fastqc_data.txt | awk -F '-' '{avg=($1+ $2)/2; printf (100*169000)/avg}')
	# echo $coverage
	rm -rf $phage_directory/fastqc/clean_1_fastqc

	seqtk sample -s 100 "$clean_read_1" "$coverage" > "$phage_directory"/sub_read_1.fastq
	seqtk sample -s 100 "$clean_read_2" "$coverage" > "$phage_directory"/sub_read_2.fastq
	
	gzip "$phage_directory"/sub_read_1.fastq
	gzip "$phage_directory"/sub_read_2.fastq
	
	source /home/ciaran/anaconda3/bin/activate spades
	
	touch "$phage_directory"/log/spades.txt
	echo Assembling the genome with SPAdes
	spades.py --only-assembler --isolate -1 "$phage_directory"/sub_read_1.fastq.gz -2 "$phage_directory"/sub_read_2.fastq.gz -o "$phage_directory"/spades_result &>>"$phage_directory"/log/spades.txt
	
	source /home/ciaran/anaconda3/bin/activate main
	
	gunzip "$phage_directory"/sub_read_1.fastq.gz
	gunzip "$phage_directory"/sub_read_2.fastq.gz
	
	touch "$phage_directory"/log/assembly_validation.txt
	echo Validating assembly
	bbmap.sh ref="$phage_directory"/spades_result/contigs.fasta in1="$phage_directory"/sub_read_1.fastq in2="$phage_directory"/sub_read_2.fastq covstats="$phage_directory"/contig_covstats.txt out="$phage_directory"/contig_mapped.sam &>> "$phage_directory"/log/assembly_validation.txt
	
	echo $coverage
	gzip "$phage_directory"/sub_read_1.fastq
	gzip "$phage_directory"/sub_read_2.fastq
	
	
	largest_contig=$(head -1 "$phage_directory"/spades_result/contigs.fasta)
	
	
	cat > "$phage_directory"/extract_contig.py <<EOF
from Bio.SeqIO.FastaIO import SimpleFastaParser 
import gzip 
import pandas as pd 
file = '$phage_directory/spades_result/contigs.fasta' 
handle = open(file) 
textfile = open("$phage_directory/phage_contig.fasta", "w") 
entries = list(SimpleFastaParser(handle)) 
for name, seq in entries: 	
	if "${largest_contig:1}" in name: 		
		print(f">{name}\n{seq}\n", file=textfile)
EOF


	python3 "$phage_directory"/extract_contig.py


	source /home/ciaran/anaconda3/bin/activate sam
	
	touch "$phage_directory"/log/sam_sort_index.txt
	echo Sorting and indexing the mapped reads for use with Pilon
	samtools view -bS -F4 "$phage_directory"/contig_mapped.sam | samtools sort - -o "$phage_directory"/contig_mapped_sorted.bam &>> "$phage_directory"/log/sam_sort_index.txt
	samtools index "$phage_directory"/contig_mapped_sorted.bam &>> "$phage_directory"/log/sam_sort_index.txt

	source /home/ciaran/anaconda3/bin/activate main
	
	touch "$phage_directory"/log/pilon.txt
	echo Polishing assembly with Pilon
	pilon --genome "$phage_directory"/phage_contig.fasta --frags "$phage_directory"/contig_mapped_sorted.bam --output "$phage_directory"/pilon_read --verbose --changes &>> "$phage_directory"/log/pilon.txt
	
	confirmed=$(awk '/^Confirmed/' $phage_directory/log/pilon.txt)
	awk '/^>/ {sub("_pilon$", "")} {print}' "$phage_directory"/pilon_read.fasta > "$phage_directory"/pilon_temp.fasta
	cat "$phage_directory"/pilon_temp.fasta > "$phage_directory"/pilon_read.fasta
	loop_number=1
	
	while true; do
		if [[ $confirmed =~ 100 ]] || [ "$loop_number" = 5 ]; then
			break
		else
			loop_number=$[loop_number+1]
			bbmap.sh ref="$phage_directory"/pilon_read.fasta in1="$phage_directory"/sub_read_1.fastq in2="$phage_directory"/sub_read_2.fastq covstats="$phage_directory"/contig_covstats.txt out="$phage_directory"/contig_mapped.sam &>> "$phage_directory"/log/assembly_validation.txt
			source /home/ciaran/anaconda3/bin/activate sam
	
			samtools view -bS -F4 "$phage_directory"/contig_mapped.sam | samtools sort - -o "$phage_directory"/contig_mapped_sorted.bam &>> "$phage_directory"/log/sam_sort_index.txt
			samtools index "$phage_directory"/contig_mapped_sorted.bam &>> "$phage_directory"/log/sam_sort_index.txt

			source /home/ciaran/anaconda3/bin/activate main
			
			touch "$phage_directory"/log/pilon$loop_number.txt
			echo Polishing assembly with Pilon
			pilon --genome "$phage_directory"/phage_contig.fasta --frags "$phage_directory"/contig_mapped_sorted.bam --output "$phage_directory"/pilon_read --verbose --changes &>> "$phage_directory"/log/pilon$loop_number.txt
			awk '/^>/ {sub("_pilon$", "")} {print}' "$phage_directory"/pilon_read.fasta > "$phage_directory"/pilon_temp.fasta
			cat "$phage_directory"/pilon_temp.fasta > "$phage_directory"/pilon_read.fasta
			confirmed=$(awk '/^Confirmed/' $phage_directory/log/pilon$loop_number.txt)
			echo $loop_number
		fi
	done
	
	 source /home/ciaran/anaconda3/bin/activate sam

	apc.pl "$phage_directory"/pilon_read.fasta

	mv permuted* "$phage_directory"

	source /home/ciaran/anaconda3/bin/activate main
	
	if [ -e $phage_directory/permuted.1.fa ];
	then
		phage_to_blast=$phage_directory/permuted.1.fa
	else
		phage_to_blast=$phage_directory/pilon_read.fasta
	fi
	
	echo 'BLASTing the genome (this will take some time)'
	blastn -query "$phage_to_blast" -remote -db nr -out "$phage_directory"/blast_result.txt -outfmt 6 -evalue 1e-30

	accession_number=$(awk -v line=3 'NR==3 {print $2}' "$phage_directory"/blast_result.txt)
	echo $accession_number
	datasets download virus genome accession "$accession_number" --include cds genome

	mv ncbi_dataset.zip "$phage_directory"

	unzip "$phage_directory"/ncbi_dataset.zip -d "$phage_directory"

	tail -n +2 "$phage_directory"/ncbi_dataset/data/cds.fna > "$phage_directory"/ncbi_dataset/data/processed_cds.fna

	first_gene_end=$(awk '/^>/ {print NR-1;exit}' "$phage_directory"/ncbi_dataset/data/processed_cds.fna)

	sed -n 1,"$first_gene_end"p "$phage_directory"/ncbi_dataset/data/processed_cds.fna > "$phage_directory"/first_gene.fasta

	start_of_phage=$(blastn -outfmt "7 sstart" -query  "$phage_directory"/first_gene.fasta -subject "$phage_directory"/permuted.1.fa | grep -v '^#')

	cat > "$phage_directory"/reorder.py <<EOF

from Bio import SeqIO

for contig_record in SeqIO.parse(open('$phage_to_blast'), 'fasta'):
    contig = str(contig_record.seq)

str1= contig[$start_of_phage:]
str2= contig[:$start_of_phage]
str_final='>'+contig_record.id + '\n' + str1 + str2
with open("$phage_directory/reordered_genome.fasta", "w") as out:
    out.write(str_final)

EOF

	python3 "$phage_directory"/reorder.py

	bbmap.sh ref="$phage_directory"/reordered_genome.fasta in1="$phage_directory"/clean_1.fastq in2="$phage_directory"/clean_2.fastq covstats="$phage_directory"/final_covstats.txt out="$phage_directory"/all_reads_contig_mapped.sam

	source /home/ciaran/anaconda3/bin/activate sam
	
	samtools view -bS -F4 "$phage_directory"/all_reads_contig_mapped.sam | samtools sort - -o "$phage_directory"/all_reads_contig_mapped_sorted.bam
	
	samtools index "$phage_directory"/all_reads_contig_mapped_sorted.bam
	
	samtools view -bS -f4 "$phage_directory"/all_reads_contig_mapped.sam | samtools sort - -o "$phage_directory"/non_mapped_sorted.bam
	
	samtools fastq "$phage_directory"/non_mapped_sorted.bam > "$phage_directory"/non_mapped_reads.fq
	
	source /home/ciaran/anaconda3/bin/activate main
	
	pilon --genome "$phage_directory"/reordered_genome.fasta --frags "$phage_directory"/all_reads_contig_mapped_sorted.bam --output "$phage_directory"/final_genome --verbose --changes
	
	sed "1s/.*/>contig/" "$phage_directory"/final_genome.fasta > "$phage_directory"/prokka_input.fasta
	
	prokka --outdir "$phage_directory"/prokka --prefix genome --kingdom viruses "$phage_directory"/prokka_input.fasta

	mkdir -p $script_dir/output/$phage_folder
	mv $phage_directory/* $script_dir/output/$phage_folder
	mv $script_dir/output/$phage_folder/1* $phage_directory
	mv $script_dir/output/$phage_folder/2* $phage_directory
	
	
done
