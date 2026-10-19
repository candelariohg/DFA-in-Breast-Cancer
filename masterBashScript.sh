#!/bin/bash


######################################################3
# Generating the list of coding genes
# #####################################################

#for N in {1..22} X; do
#for N in 22; do
#  awk -v N="$N" '
#    BEGIN {FS=OFS="\t"}
#    $6 == N && $4 == "protein_coding" {count[$2]++; names[NR]=$2}
#    END {
#      for (i=1; i<=NR; i++) {
#        if (count[names[i]] == 1)
#          print names[i]
#      }
#    }' biomart.tsv > "chr${N}CodingGenes.txt"
#done



##############################################################
# Generating the expression files 
# #############################################################

#Define subtypes
#subtypes=("Basal" "Healthy" "Her2" "LumB" "LumA" "Normal")

# Loop through chromosomes 1 to 22 and X
#for chr in {1..22} X; do
#for chr in 22; do
#    gene_file="chr${chr}CodingGenes.txt"

    # Read gene names into an array
#    mapfile -t genes < "$gene_file"

    # Loop through each subtype
#    for subtype in "${subtypes[@]}"; do
#        subtype_file="${subtype}_PBCMC_ARSYM.txt"
#        output_file="chr${chr}${subtype}Expression.tsv"

        # Extract matching rows where gene name is in first column
#        awk -v OFS="\t" -v genes="${genes[*]}" '
#            BEGIN {
#                split(genes, garr)
#                for (i in garr) gmap[garr[i]] = 1
#            }
#            gmap[$1] { print $0 }
#        ' "$subtype_file" > "$output_file"
#    done
#done



##########################################################
## Running the Mutual Information program
###########################################################

#subtypes=("Basal" "Her2" "LumA" "LumB" "Normal" "Healthy")

# Chromosomes from 1 to 22 and X
#for chr in {1..22} X; do
#for chr in 22; do
#    for subtype in "${subtypes[@]}"; do
#        filename="chr${chr}${subtypeExpression}.tsv"
#        Rscript mI.R "$filename"
#    done
#done


##########################################################
## Running the distance program 
###########################################################

#subtypes=("Basal" "Her2" "LumA" "LumB" "Normal" "Healthy")

# Chromosomes from 1 to 22 and X
#for chr in {1..22} X; do
#    for subtype in "${subtypes[@]}"; do
#        filename="chr${chr}${subtypeExpressionMI}.tsv"
#        Rscript mIandGeneDistance.R "$filename"
#    done
#done



##########################################################
## Running the program for the figures  
###########################################################

#subtypes=("Basal" "Her2" "LumA" "LumB" "Normal" "Healthy")

# Chromosomes from 1 to 22 and X
#for chr in {1..22} X; do
#    for subtype in "${subtypes[@]}"; do
#        filename="chr${chr}${subtypeExpressionMI_pairs}.tsv"
#        Rscript mIvsDistance-DistancevsMI.R "$filename"
#    done
#done


# The previous part with corrections:
#
##########################################################
## Running the program for the figures  
###########################################################

subtypes=("Basal" "Her2" "LumA" "LumB" "Normal" "Healthy")

# Chromosomes from 1 to 22 and X
for chr in {1..22} X; do
    for subtype in "${subtypes[@]}"; do
        filename="chr${chr}${subtype}.tsv"
        Rscript mIvsDistance-DistancevsMI.R "$filename"
    done
done





##########################################################
## Running the program for the other figures  
###########################################################

subtypes=("Basal" "Her2" "LumA" "LumB" "Normal" "Healthy")

# Chromosomes from 1 to 22 and X
for chr in {1..22} X; do
    for subtype in "${subtypes[@]}"; do
        filename="chr${chr}${subtypeExpressionMI_pairs}.tsv"
        Rscript mIvsOrdinal-DistanceVsOrdinal.R "$filename"
    done
done



