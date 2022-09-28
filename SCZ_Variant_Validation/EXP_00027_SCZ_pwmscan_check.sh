# use pwmscan to check all fasta regions for SCZ MPRA experiment for TFBS
conda activate JASPAR-UCSC-tracks

mkdir -p EXP_00027_SCZ_MPRA_Check
cd /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/SCZ_variant_validation

/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/profiles/get-profiles.py
ln -s /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/bin_old/pwmscan/bin /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks
/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/scan-sequence.py \
"/n/scratch3/users/t/ts267/EXP_00027_SCZ_MPRA_Check/SCZ_variant_validation.fasta" \
./profiles/ \
--output-dir ./tracks/SCZ_variant_validation/ --threads 4 --latest --taxon vertebrates

# it worked! now concatenate all the files
for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/SCZ_variant_validation -type f ! -name "*SCZ_variant_validation_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/SCZ_variant_validation/SCZ_variant_validation_all_motif.tsv.gz
done








bedtools getfasta \
-fi /n/groups/shared_databases/genomes/hg19.fa \
-bed "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/UC-Targets-hg19.100way.phastCons.Human_Brain_Vista.0.57_cut_off.50bp_merge_sort_merged.named.chr.bed" \
> "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NEx_VE_PAR_hg19_chr.fa" 

bedtools getfasta \
-fi /n/groups/shared_databases/genomes/hg19.fa \
-bed "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NE-Targets_all_chrom-Conserved_Neural_Enhancers_target_sort_merged.named.chr.bed" \
> "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NE_hg19_chr.fa" 


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NEx_VE_PAR_hg19_chr.fa" \
./profiles/ \
--output-dir ./tracks/NEx_VE_PAR_hg19/ --threads 4 --latest --taxon vertebrates



./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NE_hg19_chr.fa" \
./profiles/ \
--output-dir ./tracks/NE_hg19/ --threads 4 --latest --taxon vertebrates



for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_VE_PAR_hg19 -type f ! -name "*NEx_VE_PAR_hg19_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_VE_PAR_hg19/NEx_VE_PAR_hg19_all_motif.tsv.gz
done

for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NE_hg19 -type f ! -name "*NE_hg9_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NE_hg19/NE_hg19_all_motif.tsv.gz
done



# for background sequences

/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/NEx_PAR_VE_BiasAway_k_f8r3qmer_20220327.fa
/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/NE_kmer_shuffle_BiasAway_k_ag7l3ph1_20220327.fa
/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/HAR_kmer_shuffle_BiasAway_k_7w1a5sh4_20220327.fa


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/HAR_kmer_shuffle_BiasAway_k_7w1a5sh4_20220327.fa" \
./profiles/ \
--output-dir ./tracks/HAR_hg19_background/ --threads 4 --latest --taxon vertebrates


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/NE_kmer_shuffle_BiasAway_k_ag7l3ph1_20220327.fa" \
./profiles/ \
--output-dir ./tracks/NE_hg19_background/ --threads 4 --latest --taxon vertebrates


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/NEx_PAR_VE_BiasAway_k_f8r3qmer_20220327.fa" \
./profiles/ \
--output-dir ./tracks/NEx_PAR_VE_hg19_background/ --threads 4 --latest --taxon vertebrates




for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_hg19_background -type f ! -name "*HAR_hg19_background_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_hg19_background/HAR_hg19_background_all_motif.tsv.gz
done

for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NE_hg19_background -type f ! -name "*NE_hg19_background_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NE_hg19_background/NE_hg19_background_all_motif.tsv.gz
done

for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_PAR_VE_hg19_background -type f ! -name "*NEx_PAR_VE_hg19_background_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_PAR_VE_hg19_background/NEx_PAR_VE_hg19_background_all_motif.tsv.gz
done


# get separate sequences specifically for NEx and VE
# so it looks like NEx_VE_PAR migth acutlly just be VE
# need to make one for NEx only
bedtools getfasta \
-fi /n/groups/shared_databases/genomes/hg19.fa \
-bed "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NEx_hg19.chr.txt" \
> "/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NEx_hg19_chr.fa" 

./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Design_Files/Merged_BED_Files/NEx_hg19_chr.fa" \
./profiles/ \
--output-dir ./tracks/NEx_hg19/ --threads 4 --latest --taxon vertebrates


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/background_sequences/NEx_BiasAway_k_vos8mlhl_20220327.fa" \
./profiles/ \
--output-dir ./tracks/NEx_hg19_background/ --threads 4 --latest --taxon vertebrates


for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_hg19 -type f ! -name "*NEx_hg19_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_hg19/NEx_hg19_all_motif.tsv.gz
done

for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_hg19_background -type f ! -name "*NEx_hg19_background_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/NEx_hg19_background/NEx_hg19_background_all_motif.tsv.gz
done




##################
# human chimp variants
/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.bed
/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_INDEL.bed


# add extra bases for the motif scan
awk -F '\t' 'BEGIN{OFS="\t"} {print $1, $2-10, $3+10, $4}' /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.bed > /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.10bp_padding.bed 

bedtools getfasta \
-fi /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/genomes/hg38/hg38.fa \
-bed "/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.10bp_padding.bed"  \
> "/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.10bp_padding.fa"


./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.10bp_padding.fa" \
./profiles/ \
--output-dir ./tracks/HAR_panTro6_hg38_SNV/ --threads 4 --latest --taxon vertebrates


for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_panTro6_hg38_SNV -type f ! -name "*HAR_panTro6_hg38_SNV_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_panTro6_hg38_SNV/HAR_panTro6_hg38_SNV_all_motif.tsv.gz
done


awk -F '\t' 'BEGIN{OFS="\t"} {print $1, $2-5, $3+5, $4}' /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.bed > /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.5bp_padding.bed 

bedtools getfasta \
-fi /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/genomes/hg38/hg38.fa \
-bed "/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.5bp_padding.bed"  \
> "/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.5bp_padding.fa"

./scan-sequence.py \
"/n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/test/HAR_panTro6_hg38_SNV.5bp_padding.fa" \
./profiles/ \
--output-dir ./tracks/HAR_panTro6_hg38_SNV_5bp/ --threads 4 --latest --taxon vertebrates


for file in `find /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_panTro6_hg38_SNV_5bp -type f ! -name "*HHAR_panTro6_hg38_SNV_5bp_all_motif.tsv.gz" `
do
cat $file >> /n/groups/walsh/indData/Tae/Programs/O2/JASPAR-UCSC-tracks/tracks/HAR_panTro6_hg38_SNV_5bp/HAR_panTro6_hg38_SNV_5bp_all_motif.tsv.gz
done