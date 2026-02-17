#!/bin/bash -ue
set -euo pipefail

VCF_IN="subset-1000_AFs-females_split-multiallelic-GTmasked.vcf.gz"
VCF_OUT="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC.vcf.gz"
before_count=$(bcftools index -n "$VCF_IN") # original number of variants
OP="|" # Combine conditions with OR: EXCLUDE variant if any condition is true

# Candidate site-level rules (only enabled if the INFO tag is present)
# QD  : Qual/Depth below threshold
# DP  : Site depth below threshold (requires INFO/DP to be present)
# MQ  : Mapping quality below threshold
# FS  : Phred-scaled strand bias above threshold
# ReadPosRankSum : Read position bias less than threshold
declare -A site_conditions=(
  [QD]='INFO/QD < 2.0'
  [DP]='INFO/DP < 10'
  [MQ]='INFO/MQ < 40'
  [FS]='INFO/FS > 60'
  [ReadPosRankSum]='INFO/ReadPosRankSum < -8.0'
)

has_info() { bcftools view -h "$1" | grep -q "^##INFO=<ID=$2,"; }

{
  echo ""
  echo "=== VARIANT_QC on: $VCF_IN ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

expr_parts=("QUAL < 30") # Always include QUAL threshold (QUAL is a core VCF field, not in INFO)
echo "✓ QUAL — adding: QUAL < 30" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Add INFO-based rules only if that tag exists
for tag in "${!site_conditions[@]}"; do
  if has_info "$VCF_IN" "$tag"; then
    echo "✓ $tag (INFO) found — adding: ${site_conditions[$tag]}" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
    expr_parts+=("${site_conditions[$tag]}")
  else
    echo "x $tag (INFO) not found — skipping" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
  fi
done

expr_str="${expr_parts[0]}"
for cond in "${expr_parts[@]:1}"; do expr_str+=" $OP $cond"; done

{
  echo ""
  echo "Final filter expression:"
  echo "$expr_str"
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# HARD FILTER: remove variants where expr is TRUE
bcftools view -e "$expr_str" "$VCF_IN" -Oz -o "$VCF_OUT"
tabix -p vcf "$VCF_OUT"
after_count=$(bcftools index -n "$VCF_OUT")

{
  echo ""
  echo "✓ Filtering complete. Output: $VCF_OUT"
  echo "Removed: $(( before_count - after_count ))"
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
