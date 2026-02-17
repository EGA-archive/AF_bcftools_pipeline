#!/bin/bash -ue
set -euo pipefail

INPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC-sampleQC.vcf.gz"
OUTPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC-sampleQC-ploidy_fixed.vcf.gz"

{
  echo ""
  echo "=== FIX_PLOIDY: Creating gender.txt from metadata ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Create gender.txt file: sample<TAB>sex
# Convert M/F or 1/2 to M/F format for bcftools +fixploidy
awk -F, 'NR>1 {
  sample = $1
  sex_code = $2

  # Remove whitespace
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", sex_code)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", sample)

  # Convert to M/F format
  sex = ""
  if (sex_code == "M" || sex_code == "1") sex = "M"
  else if (sex_code == "F" || sex_code == "2") sex = "F"

  # Output: sample<TAB>sex
  if (sex != "") print sample "\t" sex
}' "meta.csv" > gender.txt

echo "✓ Gender file created" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Show sample count
n_samples=$(wc -l < gender.txt)
n_males=$(awk '$2=="M"' gender.txt | wc -l)
n_females=$(awk '$2=="F"' gender.txt | wc -l)

{
  echo "  Total samples in metadata file: $n_samples"
  echo "  Males: $n_males"
  echo "  Females: $n_females"
  echo "  Only samples that passed the QC will be used for the groupings "
  echo ""
  echo "Example gender assignments:"
  head -n 5 gender.txt
  echo ""
  echo "=== Fixing ploidy for X/Y chromosomes ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Fix ploidy 
# If GRCh37, ploidy_cmd points to the provided rules for GRCh37
# If GRCh38, ploidy_cmd points to the provided rules for GRCh38 
bcftools +fixploidy "${INPUT_VCF}" -Oz -o "$OUTPUT_VCF" -- -s gender.txt -p ploidy_grch37.txt

# Index the output
tabix -p vcf "$OUTPUT_VCF"

{
  echo "✓ Ploidy correction complete"
  echo "  Output: $OUTPUT_VCF"
  echo "  Males will have haploid genotypes (0 or 1) on X and Y"
  echo "  Females will have diploid genotypes (0/0, 0/1, 1/1) on X"
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
