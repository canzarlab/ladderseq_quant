# Input:
# $1 : kallisto index
# $2 : Path to directory with reads

# Merge reads
python interleave_fqs.py $2R1_band1.fq $2R1_band2.fq $2R1_band3.fq $2R1_band4.fq $2R1_band5.fq $2R1_band6.fq $2R1_band7.fq | sed '/^$/d' > $2R1_merged.fq
python interleave_fqs.py $2R2_band1.fq $2R2_band2.fq $2R2_band3.fq $2R2_band4.fq $2R2_band5.fq $2R2_band6.fq $2R2_band7.fq | sed '/^$/d' > $2R2_merged.fq

# Generate pseudobam file from merged reads
kallisto-ls pseudo --pseudobam -i $1 -o $2 $2R1_merged.fq $2R2_merged.fq

# Generate list of uniquely mapping reads
samtools view $2pseudoalignments.bam | awk 'NR%2 == 0 && $12=="NH:i:1" {print $1,"\t",$3,"\t",$12}' > $2uniquelyMappingReads.txt

# Estimating betas
Rscript EstimateBetas.R $2uniquelyMappingReads.txt $2
