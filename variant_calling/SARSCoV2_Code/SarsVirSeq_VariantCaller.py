""" clean, align, variant call pipeline """
""" started on 24.03.2020 """
""" Alexandra Popa, Lukas Endler enriched from https://covid19.galaxyproject.org/genomics/4-variation/#workflows """
""" Edited by Michael Martin to work directly from the BAM files """
""" https://www.ncbi.nlm.nih.gov/bioproject/PRJEB39849 """
"""  it needs:  module load htslib/1.9
				module load samtools/1.9
				module unload gcc/6.1.0
				module load gcc/7.1.0
				module load fastqc/0.11.8
				module load python/3.7.2
				module load bwa/0.7.17
				module load bcftools/1.9
				module load vcftools/0.1.13
"""
""" 6. Get Coverages """
""" 7. iVar to correct for the primers """
""" 8. Realignment Viterbi """
""" 9. Mark indel qualities """
""" 10. Call the variants with LoFreq"""

from argparse import ArgumentParser
import os
import sys
import subprocess
import re
import pypiper as pypiper
#from envbash import load_envbash
#load_envbash('/scratch/lab_bergthaler/pipelines/virSeqVariantCaller/moduleLoad.sh')


########################
### Define software path ###
########################

#bamUtils for clipping overlaping edges
bamUTIL = 'bam'
#the sequences of the 98 set of primers used to amplify the viral sequences
PRIMERSBED = './data/nCoV-2019.scheme_modif.bed'



########################
### Argument Parsing ###
########################
parser = ArgumentParser(description='Pypiper arguments VirSeq.')
parser = pypiper.add_pypiper_args(parser)
parser.add_argument('-dp', '--data_path', dest='data_path', help='path to sequencing data file')
parser.add_argument('-ref', '--ref_path', dest='ref', help='path to sequencing data file', default="./data/ref/SARS_CoV_2.fasta")

args = parser.parse_args()


##################
### Initialize ###
##################
args.sample_name = args.data_path.split('/')[-1].split('.')[0]
outfolder = '/'.join(args.data_path.split('/')[0:-1])
refSARSCoV2 = args.ref


pm = pypiper.PipelineManager(name = 'VirSeq', outfolder = outfolder, args = args) # initialize pipeline manager instance

##################
### Index ###
##################
### We downloaded the _sorted_viral_rh.bam files from SRA 
### https://www.ncbi.nlm.nih.gov/bioproject/PRJEB39849
pathVirusReheaderBAM = args.data_path
indexVirusReheaderBAM = args.sample_name + '_sorted_viral_rh.bam.bai'
pathIndexVirusReheaderBAM = os.path.join(outfolder, indexVirusReheaderBAM)
cmd = "samtools index " + pathVirusReheaderBAM  
pm.timestamp(cmd)
pm.run(cmd, pathIndexVirusReheaderBAM)

################
### iVar for primer masking ###
################
pm.timestamp('### Trim primers: ')

virusTrimmBAM = args.sample_name + '_trim.bam'
pathVirusTrimmBAM = os.path.join(outfolder, virusTrimmBAM)
pathTmpMv = os.path.join(outfolder, virusTrimmBAM)


if os.path.isfile(pathVirusTrimmBAM) :
    pm.timestamp('### Trim file already exists exist!')
else:
	cmd = "ivar trim -b "+PRIMERSBED+" -s 4 -q 0 -m 20 -i "+pathVirusReheaderBAM+" -p " + pathVirusTrimmBAM
	pm.timestamp(cmd)
	pm.run(cmd, pathVirusTrimmBAM)

tmpBamTrim = args.sample_name + 'tmpTrim'
pathTrimTmp = os.path.join(outfolder, tmpBamTrim)
virusTrimmSortBAM = args.sample_name + '_sorted_viral_rh_trim.bam'
pathVirusTrimmSortBAM = os.path.join(outfolder, virusTrimmSortBAM)
cmd = " samtools view -Shb " + pathVirusTrimmBAM +  " | samtools sort -T " + pathTrimTmp + " > " + pathVirusTrimmSortBAM
pm.timestamp(cmd)
pm.run(cmd, pathVirusTrimmSortBAM)

indexVirusTrimmSortBAM = args.sample_name + '_sorted_viral_rh_trim.bam.bai'
pathIndexVirusTrimmSortBAM = os.path.join(outfolder, indexVirusTrimmSortBAM)
cmd = "samtools index " + pathVirusTrimmSortBAM 
pm.timestamp(cmd)
pm.run(cmd, pathIndexVirusTrimmSortBAM)


################
### Stats mapping Virus ###
################
pm.timestamp('### Stats mapping virus: ')
flagstatVirusBAM = args.sample_name + '_virus.flagstat'
pathFlagStatVirus = os.path.join(outfolder, flagstatVirusBAM)
cmd = "samtools flagstat " + pathVirusTrimmSortBAM  + " >  " + pathFlagStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathFlagStatVirus)

idxstatsVirusBAM = args.sample_name + '_virus.idxstats'
pathIdxStatVirus = os.path.join(outfolder, idxstatsVirusBAM)
cmd = "samtools idxstats " + pathVirusTrimmSortBAM + " >  " + pathIdxStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathIdxStatVirus)

fullStatsVirusBAM = args.sample_name + '_virus.stats'
pathFullStatVirus = os.path.join(outfolder, fullStatsVirusBAM)
cmd = "samtools stats " + pathVirusTrimmSortBAM + " >  " + pathFullStatVirus
pm.timestamp(cmd)
pm.run(cmd, pathFullStatVirus)

################
### Get Coverages ###
################
'''
pm.timestamp('### Get Coverages: ')
coverageFileBw = args.sample_name + '.bw'
pathCoverageBw = os.path.join(outfolder, coverageFileBw)
coverageFileBed = args.sample_name + '.bed'
pathCoverageBed = os.path.join(outfolder, coverageFileBed)


#this bamCoverage from deeptools creates a nice bwa coverage file
cmd = "bamCoverage --bam " + pathVirusTrimmSortBAM + " -o "  + pathCoverageBw + " --binSize 10"
pm.timestamp(cmd)
pm.run(cmd, pathCoverageBw)
cmd = "bamCoverage --bam " + pathVirusTrimmSortBAM + " -o "  + pathCoverageBed + " --outFileFormat bedgraph --binSize 10"
pm.timestamp(cmd)
pm.run(cmd, pathCoverageBed)
'''
################
### Clip overlap ###
################
# when the fragment size and the reads lengths are very close

pm.timestamp('### Clip overlaps: ')
clipViralBam = args.sample_name + '_clip_viral.bam'
pathClipVirus = os.path.join(outfolder, clipViralBam)
statsClip = args.sample_name + '_clip_viral_stats.txt'
pathStatsClip = os.path.join(outfolder, statsClip)
cmd = bamUTIL + " clipOverlap --in " + pathVirusTrimmSortBAM + " --out " + pathClipVirus + " --stats > " + pathStatsClip
pm.timestamp(cmd)
pm.run(cmd, pathClipVirus)

##index
clipIndexFile = args.sample_name + '_clip_viral.bam.bai'
pathIdxClip = os.path.join(outfolder, clipIndexFile)
cmd = "samtools index " + pathClipVirus
pm.timestamp(cmd)
pm.run(cmd, pathIdxClip)


################
### Create FASTA consensus from BAM ###
################
pm.timestamp('### FastaConsensus: ')

maskBed = args.sample_name + '_mask.bed'
pathMaskBed = os.path.join(outfolder, maskBed)
maskFasta = args.sample_name + '_ref_mask.fasta'
pathMaskFasta = os.path.join(outfolder, maskFasta )
fqConsensus = args.sample_name + '_cons.fq'
pathFQConsensus = os.path.join(outfolder, fqConsensus)

cmd= "samtools mpileup -d 0 -uf " + refSARSCoV2 + " " +  pathClipVirus + " | bcftools call -c | vcfutils.pl vcf2fq > " + pathFQConsensus
pm.timestamp(cmd)
pm.run(cmd, pathFQConsensus)

#https://www.biostars.org/p/367626/
fastaConsensus = args.sample_name + '_cons.fasta'
pathFastaConsensus = os.path.join(outfolder, fastaConsensus)

cmd= "seqtk seq -aQ64 -q20 -n N " + pathFQConsensus + " > " +  pathFastaConsensus
pm.timestamp(cmd)
pm.run(cmd, pathFastaConsensus)



################
### LoFreq - realignment ###
################

pm.timestamp('### LoFreq viterbi: ')
loFreqViterbiFile =args.sample_name + '_real_viterbi.bam'
pathViterbiFile = os.path.join(outfolder, loFreqViterbiFile)
tmpViterbi = "tmp_" + args.sample_name + '.bam'
pathTmpViterbi = os.path.join(outfolder, tmpViterbi)
cmd = "lofreq viterbi -f " +  refSARSCoV2 + " " + pathClipVirus + " | samtools sort -T " + pathTmpViterbi + " - > " + pathViterbiFile
pm.timestamp(cmd)
pm.run(cmd, pathViterbiFile)

##index
viterbiIndexFile = args.sample_name + '_real_viterbi.bam.bai'
pathIdxViterbi = os.path.join(outfolder, viterbiIndexFile)
cmd = "samtools index " + pathViterbiFile + " >  " + pathIdxViterbi
pm.timestamp(cmd)
pm.run(cmd, pathIdxViterbi)


################
### LoFreq - mark indel qualities ###
################

pm.timestamp('### LoFreq indelqual: ')
loFreqIndelQualFile = args.sample_name + '_real_viterbi_indelQual.bam'
pathIndelQualFile = os.path.join(outfolder, loFreqIndelQualFile)
cmd = "lofreq indelqual --dindel -f " +  refSARSCoV2 + " -o " + pathIndelQualFile + " " + pathViterbiFile
pm.timestamp(cmd)
pm.run(cmd, pathIndelQualFile)

##index
indelQualIndexFile = args.sample_name + '_real_viterbi_indelQual.bam.bai'
pathIdxIndelQual = os.path.join(outfolder, indelQualIndexFile)
cmd = "samtools index " + pathIndelQualFile 
pm.timestamp(cmd)
pm.run(cmd, pathIdxIndelQual)



################
### LoFreq - call variants ###
################
pm.timestamp('### LoFreq: ')
loFreqFile = args.sample_name + '_lofreq.vcf'
pathLoFreqFile = os.path.join(outfolder, loFreqFile)
loFreqFilter = args.sample_name + '_lofreq_filter.vcf'
pathLoFreqFilter = os.path.join(outfolder, loFreqFilter)

#- q Skip any base with baseQ smaller
#-Q Skip alternate bases with baseQ smaller
#-m Skip reads with mapping qualirty smaller
#-C Test only positions having at least this coverage
#-a P-Value cutoff / significance level
cmd = "lofreq call-parallel --pp-threads 4 -f " + refSARSCoV2 +  "  -C 75  --no-default-filter  --call-indels -o " + pathLoFreqFile + " " + pathIndelQualFile 
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqFile)

cmd = "lofreq filter -i " + pathLoFreqFile + " -v 75 -a 0.01 --no-defaults -Q 90 --print-all | bcftools filter -m + -s \"HRUN_gt_3\" -e \"HRUN > 4\" > " + pathLoFreqFilter
pm.timestamp(cmd)
pm.run(cmd, pathLoFreqFilter)

################
### NORMALIZE VCF ###
################
pm.timestamp('### Normalize: ')
normFile = args.sample_name + '_lofreq_filter_norm.vcf'
pathNormFile = os.path.join(outfolder, normFile)
 
# 1) left aligns and normalizes indels and deletes identical (pos, ref, alt) variants
# this removes any identical indels that result from normalization and is conssitent 
# with how these cases are handled in the original Popa et al. analysis (SARS_LowFreqScript_automated.sh)
# however there's is less explicit and caused by multiple applications of splitting and combining multiallelic sites
# 2) splits bialleleic sites
cmd = "bcftools norm -f " + refSARSCoV2 + " -d none " + pathLoFreqFilter + " | bcftools norm -m-any > " + pathNormFile
pm.run(cmd, pathNormFile)


################
### SAMPLE REFERENCE ###
################
pm.timestamp('### Sample reference: ')
consVarFile = args.sample_name + '_cons.vcf.gz'
pathConsVarFile = os.path.join(outfolder, consVarFile)
consVarIdxFile = args.sample_name + '_cons.vcf.gz.tbi'
pathConsVarIdxFile = os.path.join(outfolder,consVarIdxFile)
refFile = args.sample_name + '_ref.fasta'
pathRefFile = os.path.join(outfolder, refFile)

# generates a reference sequence from the variant calls for each sequence
# this is different from the consensus sequence above in that we do not allow
# for ambiguous nucleotide codes
# uses all variants at > 50% nucleotide frequency
cmd = "bcftools view -Oz -i 'AF > 0.5' " + pathNormFile + " > " + pathConsVarFile
pm.run(cmd, pathConsVarFile)

cmd = "tabix " + pathConsVarFile
pm.run(cmd, pathConsVarIdxFile)

cmd = "bcftools consensus -f " + refSARSCoV2 + " -p " + args.sample_name + " " + pathConsVarFile + " > " + pathRefFile
pm.run(cmd, pathRefFile)










pm.stop_pipeline()
