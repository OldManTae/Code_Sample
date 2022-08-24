#!/bin/bash
#######################################################
# This script is to analyze tag seq data

# outline of what the script does
# 1. use known genomic features to extract barcodes
# 2. cluster barcodes so to remove singletons due to sequencing or PCR/RT error

# V4/10nt version (07/30/21)
# changed the script to only process 10nt barcodes. This should only be used for caMPRA using the synthesized barcodes
# removed scatter and gather since these are now only 10nt barcodes. Hopefully this will help with accuracy while still being relatively fast

# V3 (11/14/20)
# changed script to be snakemake compatible
# umitools group and count options

# V2 (06/19/20)
# changed so all array versions of the functions internally call for ${SLURM_ARRAY_TASK_ID} (so no need to call it in the sbatch wrap command)
#######################################################

# remove primer sites (SfiI and Luciferase)(and heterospacer is present)
function trim_primer_and_spacers(){
# take inputs
while [ "$1" != "" ]
do
case "$1" in
--read1) # read1 fastq file
shift; local read1=$1;;
--read2) # read2 fastq file
shift; local read2=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done

# variables to name output files
local filename=${read1##*/}
local output_file_prefix=${filename%_R1*}

# make output directory if it doesn't already exist
mkdir -p $output_dir

## Process with umi_tools # CGACGCTCTTCCGATCT SfiI site # GCGATCGC partial AsiSI site
umi_tools extract \
--extract-method=regex \
--stdin=${read1} \
--read2-in=${read2} \
--bc-pattern='(?P<discard_1>.{0,200})(?P<discard_2>(CGACGCTCTTCCGATCT){s<=1})(?P<umi_1>.{0,10})(?P<discard_3>(GCGATCGC){s<=1})' \
--bc-pattern2='(?P<discard_1>.{0,200})(?P<discard_2>(CGACGCTCTTCCGATCT){s<=1})(?P<umi_1>.{0,10})(?P<discard_3>(GCGATCGC){s<=1})' \
--stdout=$output_dir/${output_file_prefix}_R1.umi_tools.fastq.gz \
--read2-out=$output_dir/${output_file_prefix}_R2.umi_tools.fastq.gz \
-L $output_dir/${output_file_prefix}.umi_tools.log \
--either-read
}


# function to pull out and visualize barcode prior to mapping?
function extract_raw_barcodes(){
# take inputs
while [ "$1" != "" ]
do
case "$1" in
--read1) # read1 fastq file
shift; local read1=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done


# variables to name output files
local filename=${read1##*/}
local output_file_prefix=${filename%_R1*}

# pull out the barcodes
zcat $read1 | grep "^@" | cut -d ' ' -f 1 | awk 'BEGIN{OFS="\t"}{split($1, a, "_", seps)} {print ">"$1"\n"a[2]}' | gzip > $output_dir/${output_file_prefix}.barcodes.all.txt.gz
zcat $read1 | grep "^@" | cut -d ' ' -f 1 | awk 'BEGIN{OFS="\t"}{split($1, a, "_", seps)} length(a[2]) == 10 {print ">"$1"\n"a[2]}' | gzip > $output_dir/${output_file_prefix}.barcodes.10N.txt.gz
}

# make to custom reference
function map_reads() {
# set defaults
local threads=1

# take inputs
while [ "$1" != "" ]
do
case "$1" in
--read1) # read1 fastq file
shift; local read1=$1;;
--read2) # read2 fastq file
shift; local read2=$1;;
--threads) # number of computational threads
shift; local threads=$1;;
--ref_fasta) # reference fasta
shift; local ref_fasta=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done


# variables to name output files
local filename=${read1##*/}
local output_file_prefix=${filename%_R1*}

# make output directory if it doesn't exist
mkdir -p $output_dir

# map to custom reference
bwa mem -t $threads $ref_fasta $read1 $read2 | samtools sort -o $output_dir/${output_file_prefix}.bam -
samtools index $output_dir/${output_file_prefix}.bam
}

# filter for 10N barcodes (works with bam since they are easier to keep track of)
# fastq records consist of multiple lines, which makes parsing more difficult
# also adds in a gene tag (combining segregating and tagging, since it is possible to do it at the same time)
function segregate_barcodes_tag(){
# take inputs
while [ "$1" != "" ]
do
case "$1" in
--input_bam) # input bam file
shift; local input_bam=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
--gene) # gene name
shift; local gene=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done

# variables to name output files
local filename=${input_bam##*/}
local output_file_prefix=${filename%.bam*}

# make output directory if it doesn't exist
mkdir -p $output_dir

# filter for 10N barcodes
samtools view -H $input_bam > $output_dir/${output_file_prefix}.10N.sam
samtools view $input_bam | awk -F "\t" -v gene=$gene 'BEGIN{OFS="\t"} {split($1, a, "_", seps)} length(a[2]) == 10 {print $0 "\tBX:Z:" gene}' >> $output_dir/${output_file_prefix}.10N.sam

# convert sam file into bam file
samtools sort -o $output_dir/${output_file_prefix}.10N.bam $output_dir/${output_file_prefix}.10N.sam
samtools index $output_dir/${output_file_prefix}.10N.bam
rm $output_dir/${output_file_prefix}.10N.sam
}

# tag same barcode (using different methods)
function cluster_barcodes(){
# default value
local dist=1
# take inputs
while [ "$1" != "" ]
do
case "$1" in
--input_bam) # input bam file
shift; local input_bam=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
--dist) # hamming distance threshold
shift; local dist=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done

# variables to name output files
local filename=${input_bam##*/}
local output_file_prefix=${filename%.bam*}

# make output directory if it doesn't exist
mkdir -p $output_dir

## group reads by different modes
# unique - reads group share the exact same UMI
# percentile - reads groups share the exact same UMI. UMIs with counts < 1% of the median counts for UMIs at the same position are ignored
# cluster - identify clusters of connected UMIs (based on hamming distance threshold). Each network is a read group
# adjacency - cluster UMIs as above. For each cluster, select the node (UMI) with the highest counts. Visit all nodes one edge edge away. If all nodes have been visited, stop. Otherwise, repeat with remaining nodes until all nodes have been visited. Each step defines a read group
# directional - identify clustes of connected UMIS (based on hammig distance threshold) and umi A counts >= (2* umi B counts) -1. Each network is a read group

# from published data, methods percentile - directional all appear to work fine (just unique performs terribly)
umi_tools group \
--per-gene \
--gene-tag=BX \
--paired \
--edit-distance-threshold=$dist \
--extract-umi-method=read_id \
--method=directional \
-I $input_bam \
--group-out=$output_dir/${output_file_prefix}.grouped.dist_$dist.tsv

echo "complete!"
}

# OBSOLETE
################ alternative grouping
# note that this does not have the intended effect as it doesn't group by barcode in the final output (too collapsed)
function cluster_count_barcodes(){
# default value
local dist=1
# take inputs
while [ "$1" != "" ]
do
case "$1" in
--input_bam) # input bam file
shift; local input_bam=$1;;
--output_dir) # output directory
shift; local output_dir=$1;;
--dist) # hamming distance threshold
shift; local dist=$1;;
*)
usage
(>&2 echo "error")
exit 1
;;
esac
shift
done

# variables to name output files
local filename=${input_bam##*/}
local output_file_prefix=${filename%.bam*}

# make output directory if it doesn't exist
mkdir -p $output_dir

## group reads by different modes
# unique - reads group share the exact same UMI
# percentile - reads groups share the exact same UMI. UMIs with counts < 1% of the median counts for UMIs at the same position are ignored
# cluster - identify clusters of connected UMIs (based on hamming distance threshold). Each network is a read group
# adjacency - cluster UMIs as above. For each cluster, select the node (UMI) with the highest counts. Visit all nodes one edge edge away. If all nodes have been visited, stop. Otherwise, repeat with remaining nodes until all nodes have been visited. Each step defines a read group
# directional - identify clustes of connected UMIS (based on hammig distance threshold) and umi A counts >= (2* umi B counts) -1. Each network is a read group

# from published data, methods percentile - directional all appear to work fine (just unique performs terribly)
umi_tools count \
--per-gene \
--gene-tag=BX \
--paired \
--edit-distance-threshold=$dist \
--extract-umi-method=read_id \
--method=directional \
-I $input_bam \
-S $output_dir/${output_file_prefix}.count.dist_$dist.tsv

echo "complete!"
}








