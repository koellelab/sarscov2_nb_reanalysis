# Calling variants
This directory includes all of the code necessary to replicate the bioinformatic analyses presented in Martin & Koelle 2021 which is a Technical Comment in response to Popa & Genger et al. Science Translational Medicine 2020. Here we do not include the raw viral reads which can be downloaded from NCBI Bioproject #PRJEB39849. 

## Main script
```
run_all.sh
```
This is a basic bash script which is responsible for rerunning the entire analysis. It does not take any arguments. This analysis was conducted on a computing system which relies on the SGE job scheduler and has custom utilities to interface with the job scheduler, hencer the non-standard "SGE_Batch"  and "SGE_Array" syntax. However, It should be relatively straightward to adapt this script to output commands for the job scheduler of your choice. Furthermore, this script was run by hand and therefore does not appropriately account for the fact the jobs submitted to the job scheudler must be completed before the next step can begin. 

## Download and call variants
```
get_data.sh
```
This is the first script to be called by `run_all.sh` and is responsible for downloading the data from the SRA. It performs the following steps:
1. Create a list of accessions numbers for all samples involved in transmission pairs or longitudinal sets. This relies on having Supplementary Data File S4 saved as CSV files as well as the ENA Read Report from Bioproject #RJEB39849 in the `data/` subdirectory. The output from this script `data/accessions.csv` is included in this repository. 
2. The SRA toolkit prefetch function is used to download the uploaded BAM files for each accession numbers. Assumes the `prefect` function is in your path
3. Downloaded files are moved to the `data/seq` directory in subdirectories which match the file name e.g. `data/seq/CoV_218_S63588_sorted_viral_rh/CoV_218_S63588_sorted_viral_rh.bam`.
4. Python3 commands to run each BAM file through the variant calling pipeline (`SARSCoV2_Code/SarsVirSeq_VariantCaller.py`) are generated and saved to the file `SARSCoV2_Code/batch_variants_wuhan_hu_1.sh`. This job is submitted by the run_all.sh script. 

```
SarsVirSeq_VariantCaller.py
```
This is the main variant calling pipeline and is almost entirely the same as the variant calling pipeline provided by Popa & Genger et al. at https://github.com/Bergthalerlab/SARSCoV2_Code. It is a Python3 script which relies on PyPiper as a pipeline manager (http://pypiper.databio.org/en/latest/). It takes two arguments: `-dp`: the path to the BAM file to analyze and `-ref`: the path to the reference FASTA to use, which by default is Wuhan-Hu-1 (saved in `data/SARS_CoV_2.fasta`). It requires a number of dependencies including Python3, PyPiper (v. 0.12.1), samtools/bcftools/tabix (v. 1.10), iVar (v. 1.3), bamUtil (v. 1.0.14), seqtk (v. 1.0-r72-dirty), and lofreq (v. 2.1.5). All outputs are placed in the same directory as the input BAM file. 

## Recall reletive to donor consensus
Once variants have been called and any variants present at >50% frequency used to generate a donor-specific reference sequence, the viral BAM files are sorted by name and used to generate FASTQ files with samtools. The `run_all.sh` script creates commands to do this which are stored in `SARSCoV2_Code/generate_fastq.sh`. These commands are submited as a batch job to SGE by the run_all.sh script. 

Once the viral fastqs have been generated, we move the donor consensus sequences to the data/ref subdirectory. Donors are identified based on the `data/donor_recipient_files.csv` file which is a parsed and slightly edited version of the data presented in Popa & Genger et al. Data File S4. All donors are listed in the file `data/donor_samples.csv`. These donor references are indexed by `bwa`. 

Next, we use `bwa mem` to realign reads from each transmission pair 