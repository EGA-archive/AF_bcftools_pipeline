#!/bin/bash -ue
set -euo pipefail

INPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC-sampleQC-ploidy_fixed.vcf.gz"
OUTPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC-sampleQC-ploidy_fixed-AF_recalc.vcf.gz"

{
  echo ""
  echo "=== AF recalculating: Creating groups.txt from metadata ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Build groups file for bcftools +fill-tags:
#   <sample> <TAB> <group1,group2,...>
# 
# For each sample, we add:
#   - SEX group (MALE or FEMALE)
#   - ANCESTRY group (EUR, AFR, etc.)
#   - ANCESTRY_SEX combination (EUR_MALE, EUR_FEMALE, etc.)
#
# SEX in metadata: 1=MALE, 2=FEMALE

awk -F, 'NR>1 {
  sample = $1
  sex_code = $2
  ancestry = $3

  # Remove any whitespace
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", sex_code)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", ancestry)

  # Determine sex label - handle both M/F and 1/2 formats
  sex = ""
  if (sex_code == "M" || sex_code == "1") sex = "M"
  else if (sex_code == "F" || sex_code == "2") sex = "F"

  # Build comma-separated group list
  groups = ""

  # Add sex group
  if (sex != "") {
    groups = sex
  }

  # Add ancestry group
  if (length(ancestry) > 0) {
    if (groups != "") groups = groups ","
    groups = groups ancestry
  }

  # Add ancestry_sex combination
  if (length(ancestry) > 0 && sex != "") {
    if (groups != "") groups = groups ","
    groups = groups ancestry "_" sex
  }

  # Output: sample<TAB>groups
  print sample "\t" groups
}' "meta.csv" > groups.txt

echo "✓ Groups file created with ancestry, sex, and ancestry+sex combinations" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Show a few example lines for verification
{
  echo "Example group assignments:"
  head -n 5 groups.txt
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

{
  echo ""
  echo "=== Adding stratified allele frequencies to VCF (dropping all FORMAT/GT columns) ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Recalculate AF/AC/AN (and other tags) across:
#   - All samples (default)
#   - Per-group: MALE, FEMALE, EUR, AFR, EUR_MALE, EUR_FEMALE, etc.
# Then drop genotypes (-G) and strip any FORMAT header remnants (-x FORMAT).

bcftools +fill-tags "${INPUT_VCF}" -Ou -- -S groups.txt       | bcftools view -G -Ou       | bcftools annotate -x FORMAT       -Oz -o "$OUTPUT_VCF"

# Index the INFO-only VCF
tabix -p vcf "$OUTPUT_VCF"

{
  echo "✓ AF annotation complete with stratified groups"
  echo "  Output: $OUTPUT_VCF"
  echo ""
  echo "AF tags created for:"
  echo "  - Overall: AF, AC, AN"
  echo "  - By sex: AF_MALE, AF_FEMALE, AC_MALE, AC_FEMALE, AN_MALE, AN_FEMALE"
  echo "  - By ancestry: AF_EUR, AF_AFR, etc."
  echo "  - By ancestry+sex: AF_EUR_MALE, AF_EUR_FEMALE, etc."
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
