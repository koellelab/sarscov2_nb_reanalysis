# gets data and generates batch commands
bash get_data.sh
# submits batch commands to SGE job scheduler
SGE_Array -c "SARSCoV2_Code/batch_variants_wuhan_hu_1.sh" -b 15 --qsub_options="-q koelle@helens*"

# for each transmission pair we want to realign both the donor and the 
# recipient to the donor consensus sequence,
# thus, we need a file associating each sample wtih the desired reference
awk -F, '{print $1","$1"\n"$2","$1}' data/donor_recipient_files.csv > data/transmission_sample_realign_refs.csv
sort -u data/transmission_sample_realign_refs.csv -o data/transmission_sample_realign_refs.csv


mkdir data/seq_realign
# move and copy bam files
while IFS="" read -r p || [ -n "$p" ]
do
    sample=$(echo $p | awk -F, '{print $1}')
    ref=$(echo $p | awk -F, '{print $2}')
    dir="${sample}-${ref}"
    echo mkdir data/seq_realign/${dir}
    echo samtools sort -n data/seq/$sample/${sample}.bam -o data/seq_realign/${dir}/${sample}_name_sorted.bam
    echo samtools fastq -@ 8 data/seq_realign/${dir}/${sample}_name_sorted.bam -1 data/seq_realign/${dir}/${sample}_1.fastq -2 data/seq_realign/${dir}/${sample}_2.fastq -0 /dev/null -s /dev/null 
done < data/transmission_sample_realign_refs.csv > SARSCoV2_Code/generate_fastq.sh

# submits jobs to scheduler
SGE_Batch -c "bash SARSCoV2_Code/generate_fastq.sh" -r generate_fastq -q koelle@helens*



# gets reference sequences to generate reference fasta files
awk -F, '{print $1}' data/donor_recipient_files.csv > data/donor_samples.csv
sort -u data/donor_samples.csv -o data/donor_samples.csv
while IFS="" read -r p || [ -n "$p" ]
do 
	scp data/seq/${p}/${p}_ref.fasta  data/ref/${p}_ref.fasta    
done < data/donor_samples.csv

# indexes fasta
for i in data/ref/*.fasta; do bwa index $i; done

# generate commands to realign the transmission pairs
# to the donor consensus
while IFS="" read -r p || [ -n "$p" ]
do
    sample=$(echo $p | awk -F, '{print $1}')
    ref=$(echo $p | awk -F, '{print $2}')
    dir="${sample}-${ref}"
    echo bwa mem -k 17 -r 1.25 -M -t 4 data/ref/${ref}_ref.fasta data/seq_realign/${dir}/${sample}_1.fastq data/seq_realign/${dir}/${sample}_2.fastq \| samtools sort -o data/seq_realign/${dir}/${sample}_realign.bam \-
done < data/transmission_sample_realign_refs.csv > SARSCoV2_Code/realign.sh
SGE_Array -c SARSCoV2_Code/realign.sh -b 4 --qsub_options="-q koelle@helens*"

# FINALLY call variants
while IFS="" read -r p || [ -n "$p" ]
do
    sample=$(echo $p | awk -F, '{print $1}')
    ref=$(echo $p | awk -F, '{print $2}')
    dir="${sample}-${ref}"
    echo python3 SARSCoV2_Code/SarsVirSeq_VariantCaller.py -R -ref data/ref/${ref}_ref.fasta -dp data/seq_realign/${dir}/${sample}_realign.bam
done < data/transmission_sample_realign_refs.csv > SARSCoV2_Code/batch_variants_donor.sh
SGE_Array -c "SARSCoV2_Code/batch_variants_donor.sh" -b 12 --qsub_options="-q koelle@helens*"




