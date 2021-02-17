# get accession # for all data from 
# subjects involved in transmission pairs
# or subjects with longitudinal 
# data 
tail -n +2 data/abe2555_Data_file_S4.csv | \
 awk -F, '{print $2 "\n" $10}' | \
 sed '/^$/d' | sed 's/CeMM0//g' | sed 's/CeMM//g' | \
 awk '{print "CoV_"$0}' | \
 grep -f /dev/stdin data/filereport_read_run_PRJEB39849_tsv.txt | \
 awk -F'\t' '{print $3}' > data/accessions.csv

# gets subjects with lognitudinal data
tail -n +2 data/abe2555_Data_file_S5.csv | \
 awk -F, '{print $1}' | sed 's/CeMM0//g' | sed 's/CeMM//g' | \
 awk '{print "CoV_"$0}' |
 grep -f /dev/stdin data/filereport_read_run_PRJEB39849_tsv.txt | \
 awk -F'\t' '{print $3}' >> data/accessions.csv 

# unique accession numbers only
sort -u data/accessions.csv -o data/accessions.csv

# download uploaded BAM files from SRA
prefetch --option-file ./data/accessions.csv --type bam

for i in ./ERX*/*.bam; do 
    filename=$(echo $i | rev | cut -d/ -f1 | rev | cut -d. -f1)
    mkdir data/seq/${filename}
    scp $i data/seq/${filename}/${filename}.bam
done

# generate bash scripts
for i in data/seq/*/*_viral_rh.bam; do echo python3 SARSCoV2_Code/SarsVirSeq_VariantCaller.py -R -dp $i; done > SARSCoV2_Code/batch_variants_wuhan_hu_1.sh




