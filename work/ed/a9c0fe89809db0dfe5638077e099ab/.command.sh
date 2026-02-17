#!/bin/bash -ue
set -euo pipefail

VCF_IN="subset-1000_AFs-females_split-multiallelic.vcf.gz"
VCF_OUT="subset-1000_AFs-females_split-multiallelic-GTmasked.vcf.gz"
OP="|"

# Define candidate per‑genotype rules (only applied if the FORMAT tag exists)
#  - GQ: genotype quality below threshold
#  - DP: per‑genotype depth below threshold
#  - AD: allele balance (ALT fraction) below threshold when total AD>0
#        ALT fraction = AD[1] / (AD[0] + AD[1])

declare -A gt_conditions=(
  [GQ]='FMT/GQ < 20'
  [DP]='FMT/DP < 10'
  [AD]='GT="het" & (FMT/AD[*:0]+FMT/AD[*:1])>0 & ((FMT/AD[*:1]/(FMT/AD[*:0]+FMT/AD[*:1])) < 0.2 | (FMT/AD[*:1]/(FMT/AD[*:0]+FMT/AD[*:1])) > 0.8)'
)

{
  echo ""
  echo "=== GENOTYPE_QC on: $VCF_IN ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Build the final mask expression only with tags present in the header
gt_expr_parts=()
for tag in "${!gt_conditions[@]}"; do
  if bcftools view -h "$VCF_IN" | grep -q "^##FORMAT=<ID=${tag},"; then
    echo "✓ ${tag} (FORMAT) found — adding rule" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
    gt_expr_parts+=("${gt_conditions[$tag]}")
  else
    echo "x ${tag} (FORMAT) not found — skipping" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
  fi
done

# If no rules apply (none of the tags are present), pass-through the file unchanged

if (( ${#gt_expr_parts[@]} == 0 )); then
  {
    echo "x No FORMAT-based rules available; no masking performed."
  } >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
  cp -f "$VCF_IN" "$VCF_OUT"
  cp -f "${VCF_IN}.tbi" "${VCF_OUT}.tbi" 2>/dev/null || tabix -p vcf "$VCF_OUT"
  exit 0
fi

# Combine rules with OR (mask if ANY condition is true)
gt_expr="${gt_expr_parts[0]}"
for cond in "${gt_expr_parts[@]:1}"; do gt_expr+=" $OP $cond"; done

# Mask failing genotypes to ./.
# -t q : expression applies per-sample/per-genotype (FORMAT context)
# -n . : set GT to missing when -e expr evaluates to true
# -i   : apply the action to genotypes where the expression is TRUE -> mask variants that fail the QC 

{
  echo "Final per-genotype mask expression:"
  echo "$gt_expr"
  echo ""
  echo "Running bcftools +setGT (mask to ./.): 
  bcftools +setGT "$VCF_IN" -Oz -o "$VCF_OUT" -- -t q -n . -i "$gt_expr""
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

bcftools +setGT "$VCF_IN" -Oz -o "$VCF_OUT" -- -t q -n . -i "$gt_expr" 2>> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
tabix -p vcf "$VCF_OUT"

{
  echo "✓ Genotype masking complete."
  echo "Output: $VCF_OUT"
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
