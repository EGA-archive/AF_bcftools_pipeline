#!/bin/bash -ue
set -euo pipefail
{
  echo ""
  echo "=== Splitting multiallelic variants for subset-1000_AFs-females.gonl-EGAD00001000743-concat.vcf.gz ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# bcftools norm -m -any: split multi-allelic records into multiple lines (one ALT allele per line)
bcftools norm -m -any "subset-1000_AFs-females.gonl-EGAD00001000743-concat.vcf.gz" -Oz -o "subset-1000_AFs-females_split-multiallelic.vcf.gz"
tabix -p vcf "subset-1000_AFs-females_split-multiallelic.vcf.gz"

{
  echo "✓ Split + index done"
  echo "Output: subset-1000_AFs-females_split-multiallelic.vcf.gz"
  echo "Variants now: $(bcftools index -n "subset-1000_AFs-females_split-multiallelic.vcf.gz")"
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
