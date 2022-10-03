#!/bin/bash
################################################################################################
# Notes on required programs or dependencies (must be loaded prior to running script)
################################################################################################

## for loading virtual environment for python 2.7.12
module load gcc/6.2.0
module load python/2.7.12
# note the virtual environment below contains umi_tools (needs above modules to activate properly)
source /n/groups/walsh/indData/Tae/Virtual_Environments/O2/venv_python_2.7.12/bin/activate
# for loading modules needed for MIPgen_Pipeline
module load java/jdk-1.8u112
module load bwa/0.7.15
module load bedtools/2.26.0
module load picard/2.8.0
module load cutadapt/1.14
module load libpng/1.6.26
module load htslib/1.3.2
module load bcftools


export PATH="/n/groups/walsh/indData/Tae/Programs/O2/samtools-1.9/inst/bin:$PATH"
export PATH="/n/groups/walsh/indData/Tae/Programs/O2/gatk-4.1.2.0:$PATH"


################################################################################################
# Notes on changes
################################################################################################
## v8
# migrated to GATK 4
# apply base recalibration
# introduction of the use of functions! (helps better organize code)
# availability to perform both collapsed and uncollapsed versions of variant calling
# more flexible memory allocation
# changed bait pad to bait merge (no longer considering 100bp outside bait regions - just adds noise and artifacts)
# use of compressed variant called files
# record file size prior to deletion of intermediate files (good for debugging)
# try to retain recalibrated files (allow for removal of most intermediate files)
# deprecated 'split' and 'realign' options
# added CollectVariantCallingMetrics for each gvcf - useful for sample specific information before performing joint calling
# added XX:ParallelGCThreads limit to 1 to not use too many cores
# changed default threads to 1
# added NONE specification in backup_loc to turn off backup location
# added setting to immediately fail upon error
# added error messages for failures
# removed a bug usage() function would exit with a 0 (big bug)


## v7 add in additional fix/removal of unpaired reads

## v6
# filter out multimapping reads and fix mate information

## v5
# removed "R1/R2" and "fastq" assumption - can now handle "1/2" and "fq"

## v4
# only call variants on bait merged regions (using -L option - should lead to drastic increase in speed)
# introduce option to collapse, tag duplicates, and/or produce stats
# corrected number of threads for bamclipper
# introduced subsetting out off-targets from bam file (prior to bamclipper and collapsing or tagging)
# added in option to provide bait with padding (100bp) so that it won't be produced every time (and we won't miss regions that for some reason align just outside the bait region)
# separate out off-target

## v3
# remove need for bait/target files (produce from design file)

################################################################################################
# Fail immediately upon error
################################################################################################
set -e

################################################################################################
# Functions
################################################################################################
# perform base recalibration
base_recalibation(){
# input variables
# convert long options to short
while [ "$1" != "" ]
do
case "$1" in
--input_bam) # input vcf
shift; local input_bam=$1;;
--known_site_1) # known site file 1
shift; local known_site_1=$1;;
--known_site_2) # known site file 2
shift; local known_site_2=$1;;
--known_site_3) # known site file 3
shift; local known_site_3=$1;;
--output_prefix) # output directory
shift; local output_prefix=$1;;
--gref) # reference genome
shift; local gref=$1;;
--bait_merge) # intervals with padding
shift; local bait_merge=$1;;
--mem) # memory allocation
shift; local mem=$1;;
--threads) # number of threads allowed (for java in particular)
shift; local threads=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done

# extract read group info
local RG_ID=`samtools view -H $input_bam | grep "^@RG" | cut -d $'\t' -f 2 | cut -d ':' -f 2`
local RG_SM=`samtools view -H $input_bam | grep "^@RG" | cut -d $'\t' -f 3 | cut -d ':' -f 2`
local RG_PU=`samtools view $input_bam | head -1 | cut -d ':' -f 1`
echo "RG_ID: $RG_ID"
echo "RG_SM: $RG_SM"
echo "RG_PU: $RG_PU"

echo "reference: $gref"

gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
AddOrReplaceReadGroups \
-I ${input_bam} \
-O ${input_bam%.bam}.fixed.bam \
-RGPL illumina \
-RGID $RG_ID \
-RGSM $RG_SM \
-RGPU $RG_PU \
-RGLB $RG_ID

# use bedtools to select target regions only
samtools view \
-L $bait_merge \
-o ${input_bam%.bam}.fixed.bait_w_pad.bam \
-b ${input_bam%.bam}.fixed.bam

gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
BaseRecalibrator \
-R $gref \
-I ${input_bam%.bam}.fixed.bait_w_pad.bam \
-O ${output_prefix}.recal1.grp \
--known-sites $known_site_1 \
--known-sites $known_site_2 \
--known-sites $known_site_3 &&
echo "Base Recalibration first pass complete"

# Apply BaseRecalibration
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
ApplyBQSR \
-R $gref \
-I ${input_bam%.bam}.fixed.bait_w_pad.bam \
-bqsr ${output_prefix}.recal1.grp \
-O ${output_prefix}.recalibrated.bam &&
echo "Applied Base Recalibration"

# Test recalibration on adjusted bam to see improvement
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
BaseRecalibrator \
-R $gref \
-I ${output_prefix}.recalibrated.bam \
-O ${output_prefix}.recal2.grp \
--known-sites $known_site_1 \
--known-sites $known_site_2 \
--known-sites $known_site_3 &&
echo "Base Recalibration second pass complete"

# Assess the impact of the recalibration
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
AnalyzeCovariates \
-before ${output_prefix}.recal1.grp \
-after ${output_prefix}.recal2.grp \
-plots ${output_prefix}.recalQC.pdf &&
echo "Assessed Impact of Base Recalibration"

# remove intermediate files
rm ${input_bam%.bam}.fixed.bait_w_pad.bam
rm ${input_bam%.bam}.fixed.bam &&
echo "Removed intermediate files"
}

picard_collect_metrics(){
# input variables
# convert long options to short
while [ "$1" != "" ]
do
case "$1" in
--preclipped_bam) # preclipped input bam
shift; local preclipped_bam=$1;;
--postclipped_bam) # postclipped input bam
shift; local postclipped_bam=$1;;
--output_prefix) # output directory
shift; local output_prefix=$1;;
--gref) # reference genome
shift; local gref=$1;;
--bait_arms) # bait intervals with arms
shift; local bait_arms=$1;;
--bait) # bait intervals without arms
shift; local bait=$1;;
--target) # target intervals
shift; local target=$1;;
--mem) # memory usage
shift; local mem=$1;;
--threads) # number of threads allowed (for java in particular)
shift; local threads=$1;;
*)
(>&2 echo "error")
exit 1
;;
esac
shift
done

# Picard command (preclipped)
if [ -n "$bait_arms" ]
then
# Picard command (preclipped)
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
CollectHsMetrics \
-I ${preclipped_bam} \
-O ${output_prefix}.preclipped.picard \
-R $gref \
--BAIT_INTERVALS=$bait_arms \
--TARGET_INTERVALS=$target \
--PER_TARGET_COVERAGE=${output_prefix}.preclipped.PER_TARGET_COVERAGE.tab
else
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
CollectHsMetrics \
-I ${preclipped_bam} \
-O ${output_prefix}.preclipped.picard \
-R $gref \
--BAIT_INTERVALS=$bait \
--TARGET_INTERVALS=$target \
--PER_TARGET_COVERAGE=${output_prefix}.preclipped.PER_TARGET_COVERAGE.tab
fi

# Post clipping
gatk --java-options "-Xmx${mem} -XX:ParallelGCThreads=${threads}" \
CollectHsMetrics \
-I ${postclipped_bam} \
-O ${output_prefix}.clipped.picard \
-R $gref \
--BAIT_INTERVALS=$bait \
--TARGET_INTERVALS=$target \
--PER_TARGET_COVERAGE=${output_prefix}.clipped.PER_TARGET_COVERAGE.tab
}


# defaults
threads=1
################################################################################################
# parsing options
################################################################################################
usage() { echo "$0 Usage:" && grep ")\ #" $0; }

# convert long options to short
while [ "$1" != "" ]
do
	case "$1" in
		--loc) # file directory with fastqs
			shift; location=$1;;
		--file_name) # fastq name (without path)
			shift; input_fastq=$1;;
		--output_loc) # output directory
			shift; output_loc=$1;;
		--analysis_name) # analysis name
			shift; analysistype=$1;;
		--length) # length of barcode
			shift; length=$1;;
		--flex) # option not currently used
			shift; flex=$1;; # this option is not currently used
		--design_file) # MIP design file (without customization)
			shift; design=$1;;
		--bait) # bait file
			shift; bait=$1;;
		--bait_arms) # bait with arm sequence included (just slight improvement in preclipping metrics)
			shift; bait_arms=$1;;
		--bait_merge) # bait with arm sequence included
			shift; bait_merge=$1;;
		--target) # target file
			shift; target=$1;;
		--bedpe) # bedpe file (if not present will make from design)
			shift; bedpe=$1;;
		--call_variants) # perform variant calling
			call_variants="True";;
		--realign) # perform realignment
			realign="True";;
		--clean) # removes intermediate files
			clean="True";;
		--QC) # QC samples (Fastqc)
			QC="True";;
		--split) # split by chromosomes
			split="True";;
		--threads) # specify number of threads (default 16)
			shift; threads=$1;;
		--duplicates) # determine how to handle duplicates
			shift; duplicates=$1;;
		--recalibrate)
			recalibrate="True";;
		--mem) # memory allocated (allows for more flexible allocation for apps, based on max mem allowed)
			shift; mem=$1;;
		--backup_loc) # copy gvcf and stats to different (non-temporary) location - good for scratch drive use
			shift; backup_loc=$1;;
		*)
			usage
			(>&2 echo "error")
			exit 1
			;;
	esac
	shift
done
################################################################################################
# Print settings
################################################################################################ 
echo "Settings:
--loc location: $location
--file_name input_fastq: $input_fastq
--output_loc output_loc: $output_loc
--analysis_name analysistype: $analysistype
--length length: $length
--flex flex: $flex
--design_file design: $design
--bait bait: $bait
--bait_arms bait_arms: $bait_arms
--bait_merge bait_merge: $bait_merge
--target target: $target
--bedpe bedpe: $bedpe
--call_variants call_variants: $call_variants
--realign realign: $realign
--clean clean: $clean
--QC QC: $QC
--split split: $split
--threads threads: $threads
--duplicates duplicates: $duplicates
--recalibrate recalibrate: $recalibrate
--mem mem: $mem
--backup_loc backup_loc: $backup_loc"

################################################################################################
# settings for cluster
################################################################################################ 
# make output location with correct subfolder
# test to see if what type of naming reads follow (novagene vs macrogen)
grep -F ".fastq" <<< ${input_fastq} && fastq_test="True"
grep -F ".fq" <<< ${input_fastq} && fq_test="True"

if [ "$fastq_test" == "True" ]
then
# extract file name
input_fastq_name=${input_fastq%.fastq*}
# naming dependent of collapse status
if [ -n "$duplicates" ]
then
mkdir -p $output_loc/$analysistype/$input_fastq_name/$duplicates
output=$output_loc/$analysistype/$input_fastq_name/$duplicates
report=$output_loc/$analysistype/$input_fastq_name/$duplicates/$input_fastq_name
else
mkdir -p $output_loc/$analysistype/$input_fastq_name/no_collapse
output=$output_loc/$analysistype/$input_fastq_name/no_collapse
report=$output_loc/$analysistype/$input_fastq_name/no_collapse/$input_fastq_name
fi
elif [ "$fq_test" == "True" ]
then
# extract file name
input_fastq_name=${input_fastq%.fq*}
# naming dependent of collapse status
if [ -n "$duplicates" ]
then
mkdir -p $output_loc/$analysistype/${input_fastq_name}/$duplicates
output=$output_loc/$analysistype/${input_fastq_name}/$duplicates
report=$output_loc/$analysistype/${input_fastq_name}/$duplicates/${input_fastq_name}
else
mkdir -p $output_loc/$analysistype/${input_fastq_name}/no_collapse
output=$output_loc/$analysistype/${input_fastq_name}/no_collapse
report=$output_loc/$analysistype/${input_fastq_name}/no_collapse/${input_fastq_name}
fi
else
echo "can not detect fastq or fq suffix"
exit 1
fi

# test if using _R1. or _1. naming
grep -F "_R1." <<< ${input_fastq} && R1_test="True"
grep -F "_1." <<< ${input_fastq} && nonR1_test="True"

if [ "$R1_test" == "True" ]
then
# read 2 (assumes that its is equivalent to R1, but with just R2 instead)
input_fastq2=`echo $input_fastq | sed -e 's/_R1./_R2./g'`
ID=`echo $input_fastq2 | sed -e 's/_R2.*//g'`
SM=`echo $input_fastq2 | sed -e 's/_R2.*//g'`
elif [ "$nonR1_test" == "True" ]
then
# read 2 (assumes that its is equivalent to R1, but with just R2 instead)
input_fastq2=`echo $input_fastq | sed -e 's/_1./_2./g'`
ID=`echo $input_fastq2 | sed -e 's/_2.*//g'`
SM=`echo $input_fastq2 | sed -e 's/_2.*//g'`
else
echo "Can not detect R1 or 1 read number"
exit 1
fi

### Variables that depend on whether in O2 or Orchestra
if [[ "$HMS_CLUSTER" = o2 ]]
then
gref=/n/groups/walsh/ryan/genomes/g1k/human_g1k_v37_decoy.fasta
dict=/n/groups/walsh/ryan/genomes/g1k/human_g1k_v37_decoy.dict
mod_dir=/n/groups/walsh/indData/Tae/scripts
GATK=/n/groups/walsh/indData/Tae/scripts/GATK/GenomeAnalysisTK-3.7.jar
picard_tool=$PICARD/picard-2.8.0.jar
variant_calling="sbatch -p medium -t 1-18:00 -c 2 -N 1 --mem=20000 -o $report.out -e $report.err --wrap"
elif [[ "$HMS_CLUSTER" = orchestra ]]
then
gref=/groups/walsh/ryan/genomes/g1k/human_g1k_v37_decoy.fasta
mod_dir=/groups/walsh/indData/Tae/scripts
GATK=/groups/walsh/indData/Tae/scripts/GATK/GenomeAnalysisTK-3.7.jar
picard_tool=/opt/picard-2.7.1/bin/picard.jar
variant_calling="bsub -q mcore -W 42:00 -R \\\"rusage[mem=20000]\\\" -R \\\"select[scratch2]\\\" -R \\\"select[groups]\\\" -n 2 -o $report.1.out -e $report.1.err"
fi
echo "checking/generating necessary region files..."
# determine bedpe file usage (make or use one)
if [ -z "$bedpe" ]
then
# if bedpe is not provided, make bedpe from design file
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {if ($4 <= $8) print $3, $4, $5, $3, $8, $9, $1; else print $3, $8, $9, $3, $4, $5, $1;}' $design | sort -k1,1 -k2,2n > $output/${design##*/}.bedpe
bedpe=$output/${design##*/}.bedpe
fi

if [ -z "$bait" ]
then
# if bait is not provided, make bait from design file
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {print $3, $12, $13;}' $design | sort -k1,1 -k2,2n | uniq > $output/${design##*/}.bait.bed
# make interval file
java -Xmx2g -XX:ParallelGCThreads=${threads} -jar $picard_tool BedToIntervalList \
I=$output/${design##*/}.bait.bed \
O=$output/${design##*/}.bait.interval_list \
SD=$dict
bait=$output/${design##*/}.bait.interval_list
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {if ($4 <= $8) print $3, $4, $9; else print $3, $8, $5;}' $design | sort -k1,1 -k2,2n | uniq > $output/${design##*/}.bait_with_arms.bed
java -Xmx2g -XX:ParallelGCThreads=${threads} -jar $picard_tool BedToIntervalList \
I=$output/${design##*/}.bait_with_arms.bed \
O=$output/${design##*/}.bait_with_arms.interval_list \
SD=$dict
bait_arms=$output/${design##*/}.bait_with_arms.interval_list
elif [ -z "$bait_arms" ]
then
# if bait provided but not bait arms, then make from design file (these are bait files with arms retained) 
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {if ($4 <= $8) print $3, $4, $9; else print $3, $8, $5;}' $design | sort -k1,1 -k2,2n | uniq > $output/${design##*/}.bait_with_arms.bed
java -Xmx2g -XX:ParallelGCThreads=${threads} -jar $picard_tool BedToIntervalList \
I=$output/${design##*/}.bait_with_arms.bed \
O=$output/${design##*/}.bait_with_arms.interval_list \
SD=$dict
bait_arms=$output/${design##*/}.bait_with_arms.interval_lists
fi

if [ -z "$target" ]
then
# if target is not provided, make target from design file
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {print $3, $16, $17;}' $design | sort -k1,1 -k2,2n | uniq > $output/${design##*/}.target.bed
bedtools merge -i $output/${design##*/}.target.bed > $output/${design##*/}.target.merge.bed
java -Xmx2g -XX:ParallelGCThreads=${threads} -jar $picard_tool BedToIntervalList \
I=$output/${design##*/}.target.merge.bed \
O=$output/${design##*/}.target.interval_list \
SD=$dict
target=$output/${design##*/}.target.interval_list
fi

# Make bed file to call variants only at bait + surrounding region
if [ -z "$bait_merge" ]
then
awk -F '\t' 'BEGIN{OFS="\t";} NR == 1 {next;} {if ($4 <= $8) print $3, $4, $9; else print $3, $8, $5;}' $design | sort -k1,1 -k2,2n | uniq > $output/${design##*/}.bait.region.bed
bedtools merge -i $output/${design##*/}.bait.region.bed > $output/${design##*/}.bait.region.merge.bed
bait_merge=$output/${design##*/}.bait.region.merge.bed
fi
################################################################################################
# QC of sequencing data
################################################################################################
mkdir -p $output/stats
# Fastqc option
if [ "$QC" == "True" ]
then
	echo "generating QC metrics for fastq files..."
	if [[ $read1 == *.gz ]]
	then
		zcat $location/${input_fastq} | fastqc stdin --outdir=$output/stats/ &&
		zcat $location/${input_fastq2} | fastqc stdin --outdir=$output/stats/
	elif [[ $read1 == *.bz2 ]]
	then
		bzip2 -dkc $location/${input_fastq} | fastqc stdin --outdir=$output/stats/ &&
		bzip2 -dkc $location/${input_fastq2} | fastqc stdin --outdir=$output/stats/
	else
		fastqc $location/${input_fastq} --outdir=$output/stats/ &&
		fastqc $location/${input_fastq2} --outdir=$output/stats/
	fi	
fi

################################################################################################
# Pre-processing reads
################################################################################################
# Quality trimming and adapter removal
echo "Pergoming Quality trimming and adapter removal..."
cutadapt -q 15 -A XXXXX -o $output/$input_fastq_name.qual.fastq.gz -p $output/${input_fastq2%.fastq*}.qual.fastq.gz $location/${input_fastq} $location/${input_fastq2} > $output/${input_fastq_name}_quality_trim.txt
cutadapt -e 0.1 -m $(($length+10)) -g CATACGAGATCCGTAATCGGGAAGCTGAAG -a ACACTACCGTCGGATCGTGCGTGT -G GCTAAGGGCCTAACTGGCCGCTTCACTG -A CTTCAGCTTCCCGATTACGGATCTCGTATG --too-short-output $output/${input_fastq_name}.short.fastq.gz --too-short-paired-output $output/${input_fastq2%.fastq*}.short.fastq.gz -o $output/${input_fastq_name}.trim.fastq.gz -p $output/${input_fastq2%.fastq*}.trim.fastq.gz $output/${input_fastq_name}.qual.fastq.gz $output/${input_fastq2%.fastq*}.qual.fastq.gz > $output/${input_fastq_name}_adapter_trim.txt

# Use umitools to extract barcode
echo "using umitools to extract barcode..."
barcode=`printf 'N%.0s' $(seq 1 $length)`

# With regular MIPs the barcode is in the 5' end of Read 2 (so flip around in umitools so it can handle it correctly)
umi_tools extract --extract-method=string --stdin=$output/${input_fastq2%.fastq*}.trim.fastq.gz --read2-in=$output/${input_fastq_name}.trim.fastq.gz --bc-pattern=$barcode --stdout=$output/${input_fastq2%.fastq*}.umitools.fastq.gz --read2-out=$output/${input_fastq_name}.umitools.fastq.gz -L $output/${input_fastq2%.fastq*}.log

# remove files that are too short (short fragments)
echo "removing truncated reads..."
cutadapt -m $(($length+10)) --too-short-output $output/${input_fastq_name}.umitools.short.fastq.gz --too-short-paired-output $output/${input_fastq2%.fastq*}.umitools.short.fastq.gz -o $output/${input_fastq_name}.umitools.filtered.fastq.gz -p $output/${input_fastq2%.fastq*}.umitools.filtered.fastq.gz $output/${input_fastq_name}.umitools.fastq.gz $output/${input_fastq2%.fastq*}.umitools.fastq.gz > $output/${input_fastq_name}_too_short_trim.txt

# map reads (-M option for picard compatibility)
echo "making reads..."
bwa mem -t $threads -M -Y -R "@RG\tID:$ID\tSM:$SM\tPL:illumina" $gref $output/${input_fastq_name}.umitools.filtered.fastq.gz $output/${input_fastq2%.fastq*}.umitools.filtered.fastq.gz > $output/${input_fastq_name}.indexed.sam &&
samtools view -bS $output/${input_fastq_name}.indexed.sam | samtools sort -o $output/${input_fastq_name}.sort.bam &&
samtools index $output/${input_fastq_name}.sort.bam && rm $output/${input_fastq_name}.indexed.sam

# remove off-targets
echo "removing off-targets..."
bedtools intersect -a $output/${input_fastq_name}.sort.bam -b $bait_merge | samtools sort -o $output/${input_fastq_name}.ontarget.bam &&
samtools index $output/${input_fastq_name}.ontarget.bam

# keep off-targets for closer look
echo "retaining targets for further inspections..."
bedtools intersect -v -a $output/${input_fastq_name}.sort.bam -b $bait_merge | samtools sort -o $output/${input_fastq_name}.offtarget.bam &&
samtools index $output/${input_fastq_name}.offtarget.bam

# Clip extension and ligation arms
echo "clipping off extension and ligation arms..."
cd $output
sh /n/groups/walsh/ryan/scripts/bamclipper-master/bamclipper.sh -b $output/${input_fastq_name}.ontarget.bam -p $bedpe -n $threads -g /n/groups/walsh/indData/Tae/Programs/O2/parallel-20180122/src/parallel -s /n/app/samtools/1.3.1/bin/samtools -u 15 -d 15

# for unmapped reads post clipping, change MAPQ score to 0
echo "adjusting MAPQ score post clipping..."
samtools view -h $output/${input_fastq_name}.ontarget.primerclipped.bam | awk -F '\t' 'BEGIN{OFS="\t"} {if ($6 ~ /[*]/) sub(/.*/,0,$5); print $0;}' | samtools view -bS | samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.fix.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.fix.bam 

# filter out multimapping reads and fix mate information
echo "filtering oute multimapping reads and fixing mate information..."
java -XX:ParallelGCThreads=${threads} -jar $PICARD/picard-2.8.0.jar FixMateInformation \
I=$output/${input_fastq_name}.ontarget.primerclipped.fix.bam \
O=$output/${input_fastq_name}.ontarget.primerclipped.fixmate.bam \
ADD_MATE_CIGAR=true

# apply sample tools fixmate (for identifying newly unpaired reads)
samtools sort -n $output/${input_fastq_name}.ontarget.primerclipped.fixmate.bam -o $output/${input_fastq_name}.ontarget.primerclipped.fixmate.nsort.bam && \
samtools fixmate $output/${input_fastq_name}.ontarget.primerclipped.fixmate.nsort.bam $output/${input_fastq_name}.ontarget.primerclipped.fixmate.nsort.samfix.bam

# remove unmapped reads
echo "removing unmapped reads after fixes..."
samtools view -F 1804 -f 1 -b $output/${input_fastq_name}.ontarget.primerclipped.fixmate.nsort.samfix.bam | samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam

# determine how to handle duplicates (if at all)
if [ "$duplicates" == "collapse" ]
then
echo "collapsing reads..."
# collapse reads
umi_tools dedup --paired -I $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam -S $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam --extract-umi-method=read_id --method=directional
# sort and index
samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam && rm $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam
processed_bam_name=${processed_bam##*/}
elif [ "$duplicates" == "tag" ]
then
echo "tagging reads..."
# tag reads
umi_tools group --paired -I $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam --output-bam -S $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam --extract-umi-method=read_id --method=directional
# sort and index
samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam && rm $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam
processed_bam_name=${processed_bam##*/}
# note _stats is useful but more time and memory intensive
elif [ "$duplicates" == "collapse_stats" ]
then
echo "collapsing reads and generating stats (resource intensive)..."
# collapse reads
umi_tools dedup --paired -I $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam -S $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam --extract-umi-method=read_id --method=directional --output-stats=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.tsv
# sort and index
samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam && rm $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam
processed_bam_name=${processed_bam##*/}
# note _stats is useful but more time and memory intensive
elif [ "$duplicates" == "dedup_stats" ]
then
echo "tagging reads and generating stats (resource intensive)..."
# tag reads
umi_tools group --paired -I $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam --output-bam -S $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam --extract-umi-method=read_id --method=directional --group-out=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.tsv
# sort and index
samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam && rm $output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.bam
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.tagged.sort.bam
processed_bam_name=${processed_bam##*/}
elif [ "$duplicates" == "both" ]
then
echo "generating collapsed and uncollapsed reads..."
# tag reads
# collapse reads
umi_tools dedup --paired -I $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam -S $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam --extract-umi-method=read_id --method=directional
# sort and index
samtools sort -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam &&
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam && rm $output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.bam
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.collapsed.sort.bam
unprocessed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam
processed_bam_name=${processed_bam##*/}
unprocessed_bam_name=${unprocessed_bam##*/}
else
echo "leaving reads uncollapsed..."
# bams to be used down stream
pre_clip_bam=$output/${input_fastq_name}.sort.bam
processed_bam=$output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam
fi

################################################################################################
# Collect metrics
################################################################################################
# Collect Metrics
# Picard
# Post Barcode Extraction
if [ "$duplicates" == "both" ]
then
echo "generating metrics for both uncollapsed and collapsed reads..."
## note the specific use of uncollapsed and collapsed is only for when both are requested
# uncollapsed
picard_collect_metrics \
--preclipped_bam $pre_clip_bam \
--postclipped_bam $unprocessed_bam \
--output_prefix $output/stats/${unprocessed_bam_name%.bam} \
--gref $gref \
--bait_arms $bait_arms \
--bait $bait \
--target $target \
--mem $mem \
--threads $threads

# collapsed
picard_collect_metrics \
--preclipped_bam $pre_clip_bam \
--postclipped_bam $processed_bam \
--output_prefix $output/stats/${processed_bam_name%.bam} \
--gref $gref \
--bait_arms $bait_arms \
--bait $bait \
--target $target \
--mem $mem \
--threads $threads

# remove one of the preclipped files since not necessary (performed twice)
rm $output/stats/${unprocessed_bam_name%.bam}.preclipped.picard
rm $output/stats/${unprocessed_bam_name%.bam}.preclipped.PER_TARGET_COVERAGE.tab
else
echo "generating metrics for uncollapsed reads..."
picard_collect_metrics \
--preclipped_bam $pre_clip_bam \
--postclipped_bam $processed_bam \
--output_prefix $output/stats/${processed_bam_name%.bam} \
--gref $gref \
--bait_arms $bait_arms \
--bait $bait \
--target $target \
--mem $mem \
--threads $threads
fi


################################################################################################
# Perform Base recalibration
################################################################################################
if [ "$recalibrate" == "True" ]
then
# if both collapsed and uncollapsed
if [ "$duplicates" == "both" ]
then
echo "performing base recalibration for both uncollapsed and collapsed reads"
# base recalibration for unprocessed bam
base_recalibation \
--input_bam ${unprocessed_bam} \
--known_site_1 /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
--known_site_2 /n/groups/walsh/indData/Tae/References/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
--known_site_3 /n/groups/walsh/indData/Tae/References/GATK/1000G_phase1.indels.b37.vcf.gz \
--output_prefix ${unprocessed_bam%.bam} \
--gref $gref \
--bait_merge $bait_merge \
--mem $mem \
--threads $threads

# base recalibration for processed bam
base_recalibation \
--input_bam ${processed_bam} \
--known_site_1 /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
--known_site_2 /n/groups/walsh/indData/Tae/References/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
--known_site_3 /n/groups/walsh/indData/Tae/References/GATK/1000G_phase1.indels.b37.vcf.gz \
--output_prefix ${processed_bam%.bam} \
--gref $gref \
--bait_merge $bait_merge \
--mem $mem \
--threads $threads

# prepare input for variant calling
unprocessed_recalibrated_bam=${unprocessed_bam%.bam}.recalibrated.bam
processed_recalibrated_bam=${processed_bam%.bam}.recalibrated.bam
else
echo "performing base recalibration for collapsed reads..."
# base recalibration for processed bam
base_recalibation \
--input_bam ${processed_bam} \
--known_site_1 /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
--known_site_2 /n/groups/walsh/indData/Tae/References/GATK/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \
--known_site_3 /n/groups/walsh/indData/Tae/References/GATK/1000G_phase1.indels.b37.vcf.gz \
--output_prefix ${processed_bam%.bam} \
--gref $gref \
--bait_merge $bait_merge \
--mem $mem \
--threads $threads

# prepare input for variant calling
processed_recalibrated_bam=${processed_bam%.bam}.recalibrated.bam
fi
else
echo "skipping base recalibration"
# if not recalibrating, then pass on processed file
processed_recalibrated_bam=${processed_bam}
fi

################################################################################################
# Call Variants (actually running is optional (but will at least prepare jobs))
################################################################################################ 
# Variant calling
# Make independent jobs that can be submitted or used to rerun
if [ "$duplicates" == "both" ]
then
echo "preparing variant calling for both uncollapsed and collapsed reads..."
# uncollapsed
echo "gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
HaplotypeCaller \
-R $gref \
-I $unprocessed_recalibrated_bam \
-O ${unprocessed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-ERC GVCF \
--dont-use-soft-clipped-bases \
-G StandardAnnotation \
-G AS_StandardAnnotation \
-G StandardHCAnnotation \
-L $bait_merge && \
gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
CollectVariantCallingMetrics \
--DBSNP /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
-I ${unprocessed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-O ${unprocessed_recalibrated_bam%.bam} \
--GVCF_INPUT True \
-TI $target" > ${unprocessed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh

# collapsed
echo "gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
HaplotypeCaller \
-R $gref \
-I $processed_recalibrated_bam \
-O ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-ERC GVCF \
--dont-use-soft-clipped-bases \
-G StandardAnnotation \
-G AS_StandardAnnotation \
-G StandardHCAnnotation \
-L $bait_merge && \
gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
CollectVariantCallingMetrics \
--DBSNP /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
-I ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-O ${processed_recalibrated_bam%.bam} \
--GVCF_INPUT True \
-TI $target" > ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh
else
echo "preparing variant calling for collapsed reads..."
echo "gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
HaplotypeCaller \
-R $gref \
-I $processed_recalibrated_bam \
-O ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-ERC GVCF \
--dont-use-soft-clipped-bases \
-G StandardAnnotation \
-G AS_StandardAnnotation \
-G StandardHCAnnotation \
-L $bait_merge && \
gatk --java-options '-Xmx${mem} -XX:ParallelGCThreads=${threads}' \
CollectVariantCallingMetrics \
--DBSNP /n/groups/walsh/indData/Tae/References/GATK/dbsnp_138.b37.vcf.gz \
-I ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.gz \
-O ${processed_recalibrated_bam%.bam} \
--GVCF_INPUT True \
-TI $target" > ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh
fi

# Variant calling within script
if [ "$call_variants" == "True" ]
then
echo "perfoming variant calling..."
[ -f ${unprocessed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh ] && sh ${unprocessed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh 
[ -f ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh ] && sh ${processed_recalibrated_bam%.bam}.snps.indels.g.vcf.sh
fi

################################################################################################
# Perform local realignment (DEPRECATED)
################################################################################################

# local realignment for downstream use of mutect or similar programs (very)
if [ "$realign" == "True" ]
then
echo "perfoming local realignment..."
java -jar $GATK -R $gref -T RealignerTargetCreator -I $processed_recalibrated_bam -known /n/scratch2/rd140/gnomad/gnomad.genomes.indels.vcf -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam.intervals &&
java -jar $GATK -R $gref -T IndelRealigner -I $processed_recalibrated_bam -known /n/scratch2/rd140/gnomad/gnomad.genomes.indels.vcf -targetIntervals $output/${input_fastq_name}.ontarget.primerclipped.cleaned.bam.intervals -o $output/${input_fastq_name}.ontarget.primerclipped.cleaned.realign.bam &&

samtools view -H $processed_bam > $output/${input_fastq_name}.header.txt
samtools reheader $output/${input_fastq_name}.header.txt $output/${input_fastq_name}.ontarget.primerclipped.cleaned.realign.bam > $output/${input_fastq_name}.ontarget.primerclipped.cleaned.realign.bam
samtools index $output/${input_fastq_name}.ontarget.primerclipped.cleaned.realign.bam
fi

################################################################################################
# Remove all intermediate files (optional)
################################################################################################
# regardless of cleaning, record file sizes of all files at this point
echo "generating record of file sizes..."
find $output -maxdepth 1 -mindepth 1 -type f -exec ls -l {} + > $output/stats/${processed_bam_name%.bam}.file_size_record.txt

# retain different files, depending on whether the files was recalibrated or not
if [ "$clean" == "True" ]
then
if [ "$recalibrate" == "True" ]
then
echo "deleting intermediate files (keeping recalibrated files)..."
find $output -maxdepth 1 -mindepth 1 -type f ! -name "*.recalibrated.*" ! -name "*.pdf" ! -name "*.out" ! -name "*.err" ! -name "*.log" ! -name "*offtarget*" -delete
else
echo "deleting intermediate files (keeping primerclipped files)..."
find $output -maxdepth 1 -mindepth 1 -type f ! -name "*.primerclipped.*" ! -name "*.pdf" ! -name "*.out" ! -name "*.err" ! -name "*.log" ! -name "*offtarget*" -delete
fi
fi

################################################################################################
# Backup final output and stats (optional, but recommended if using scratch drive)
################################################################################################
if [ -n "$backup_loc" ] && [ "$backup_loc" != "NONE" ]
then
echo "generating backup copy..."
mkdir -p ${backup_loc}/${input_fastq_name}
# save stats and record keeping
find $output -mindepth 1 -type f ! -name "*.bam" ! -name "*.sam" ! -name "*.fastq*" ! -name "*.fq*" | xargs cp -t ${backup_loc}/${input_fastq_name}/
# save final preprocessing output
find $output -mindepth 1 -type f \( -name "${unprocessed_recalibrated_bam##*/}" -o -name "${processed_recalibrated_bam##*/}" \) | xargs cp -t ${backup_loc}/${input_fastq_name}/
fi
################################################################################################
# Split by Chromosomes (optional) - note that this is not necessary given -L option in GATK, unless file is very large (DEPRECATED)
################################################################################################
if [ "$split" == "True" ]
then
echo "WARNING: USING DEPRECATED OPTION - LIKELY INCOMPATIBLE WITH PIPELINE"
exit
chrom_list=`cat $output/${input_fastq_name}.ontarget.primerclipped.cleaned.snps.indels.g.vcf | grep -v "^#" | cut -d $'\t' -f 1 | grep -v "^[MGNh]" | sort | uniq`
# for autosomal and sex chromosomes (separate into different folders)
for file in `find $project_folder -maxdepth 1 -type f -name "*g.vcf"`
do
filename=${file##*/}
for chrom in $chrom_list
do
split_output=$output/split/$chrom
mkdir -p $split_output
# add header
grep "^#" $file > $split_output/${filename%.g.vcf*}.chrom_$chrom.g.vcf
# add chromosome specific info
awk -F '\t' -v chrom=$chrom '$1 == chrom {print $0;}' $file >> $split_output/${filename%.g.vcf*}.chrom_$chrom.g.vcf &&
echo "chromosome $chrom of $file separated"
done
# for decoy, MT sequences
output_N=$output/split/N
mkdir -p $output_N
grep "^#" $file > $output_N/${filename%.g.vcf*}.chrom_N.g.vcf
grep -v "^#" $file | grep "^[MGNh]" >> $output_N/${filename%.g.vcf*}.chrom_N.g.vcf && 
echo "chromosome N of $file separated"
done
fi
