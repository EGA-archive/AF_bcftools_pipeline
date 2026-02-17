#!/bin/bash -ue
# Fail fast + propagate errors in pipelines
set -euo pipefail

INPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC.vcf.gz"
OUTPUT_VCF="subset-1000_AFs-females_split-multiallelic-GTmasked-variantQC-sampleQC.vcf.gz"

# Workspace for intermediate per-sample metrics
TMP="qc_tmp"; mkdir -p "$TMP"

{
  echo ""
  echo "=== SAMPLE_QC on: $INPUT_VCF ==="
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# -----------------------------
# Choose thresholds by seq type
# -----------------------------
if [[ "WGS" == "WGS" ]]; then
  COV_THRESHOLD=15
  HET_HOM_THRESHOLD=3.3
  CALL_RATE_THRESHOLD=0.95
  SINGLETONS_THRESHOLD=100000
  CONTAM_THRESHOLD=0.05
  echo "✓ Seq type: WGS thresholds applied" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
else
  COV_THRESHOLD=10
  HET_HOM_THRESHOLD=10
  CALL_RATE_THRESHOLD=0.95
  SINGLETONS_THRESHOLD=5000
  CONTAM_THRESHOLD=0.00015
  echo "✓ Seq type: WES thresholds applied" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
fi

# Generate bcftools stats for the input VCF
echo "Generating bcftools stats..." >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
bcftools query -l "$INPUT_VCF" > list-of-samples.txt
bcftools stats -S list-of-samples.txt "$INPUT_VCF" > bcftools-stats.txt

# --- detect FORMAT presence in the VCF header ---
DP_FOUND=$(bcftools view -h "$INPUT_VCF" | grep -c '^##FORMAT=<ID=DP,' || true)
GT_FOUND=$(bcftools view -h "$INPUT_VCF" | grep -c '^##FORMAT=<ID=GT,' || true)


if [[ "$DP_FOUND" -gt 0 ]]; then
  echo "✓ DP format found - Performing coverage filtering" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
else
  echo "x DP format NOT found - Coverage filtering NOT performed" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
fi
if [[ "$GT_FOUND" -gt 0 ]]; then
  echo "✓ GT format found - Performing call rate, het/hom ratio and singleton filtering" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
else
  echo "x GT format NOT found - Call rate, het/hom ratio and singleton filtering NOT performed" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
fi

# 1) Extract SN + PSC sections from bcftools stats
awk '
  printing && /^#/ { printing=0 }
  /^# SN/ { printing=1 }
  /^# PSC[[:space:]]+\[2\]id[[:space:]]+\[3\]sample/ { printing=1 }
  printing
' bcftools-stats.txt > bcftools-stats_min.txt

# 2) Read number of samples and number of records from SN section
read nsamples_orig nvariants < <(
  awk '
    $3=="number" && $4=="of" && $5=="samples:" { ns=$6 }
    $3=="number" && $4=="of" && $5=="records:" { nr=$6 }
    END { print ns, nr }
  ' bcftools-stats_min.txt
)

# 3) Compute total genotypes
ngenotypes=$(( nsamples_orig * nvariants ))

# 4) Add rHetHom [15] and CallRate [16] to PSC section
awk -v OFS="\t" -v ngenotypes="$ngenotypes" '
  # Extend PSC header
  /^# PSC[[:space:]]+\[2\]id[[:space:]]+\[3\]sample[[:space:]]+\[4\]nRefHom/ {
    in_psc=1
    print $0, "[15] rHetHom", "[16] CallRate"
    next
  }

  # Any other header stops PSC mode
  /^#/ { in_psc=0; print; next }

  # PSC rows: compute new metrics
  in_psc && $1=="PSC" {
    r = ($5+0)==0 ? "NA" : sprintf("%.6f", $6/$5)        # [15]
    cr = sprintf("%.6f", (1 - ($14/ngenotypes)))         # [16], using [14]=nMissing
    print $0, r, cr
    next
  }

  { print }
' bcftools-stats_min.txt > sample-qc-stats.txt

# 5) Contamination check 
CHARR_TSV="sceVCF-results.tsv"
: > "$CHARR_TSV"  # create empty file by default

if bcftools view -h "$INPUT_VCF" | grep -q "^##FORMAT=<ID=AD,"; then
  if [[ -n "./sceVCF" ]]; then
    SCE_CMD=""
    if [[ -x "./sceVCF" ]]; then
      SCE_CMD="./sceVCF"
    elif [[ -x "./sceVCF/sceVCF" ]]; then
      SCE_CMD="./sceVCF/sceVCF"
    fi

    if [[ -n "$SCE_CMD" ]]; then
      echo "✓ sceVCF found ($SCE_CMD) — running contamination check" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
      "$SCE_CMD" -o "$CHARR_TSV" "$INPUT_VCF"
      awk -F'	' -v OFS='	' '
      FNR==NR {
        # Read sceVC results into array (skip meta)
        if ($1 ~ /^#/ || $1 ~ /^##/) next
        charr[$1] = $10
        next
      }
      # When we hit the #PSC header, append [17] CHARR
      /^# PSC/ {
        print $0, "[17] CHARR"
        next
      }
      # For PSC data lines, append CHARR value
      /^PSC/ {
        val = ($3 in charr ? charr[$3] : "NA")
        print $0, val
        next
      }
      # Print all other lines unchanged
      { print }
    ' charr_full.tsv sample-qc-stats.txt > sample-qc-stats.txt.tmp && mv sample-qc-stats.txt.tmp sample-qc-stats.txt

    else
      echo "x sceVCF not found or not executable at: ./sceVCF — skipping" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
    fi
  else
    echo "x Contamination check not running (sceVCF_path empty)" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
  fi
else
  echo "x AD (FORMAT) not found — skipping contamination check (required for sceVCF)" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
fi


# --- now filter samples from the PSC section according to your rules ---
# Rules:
# - If DP present:        [10] average depth < COV_THRESHOLD           -> fail "low_cov"
# - If GT present:        [15] rHetHom > HET_HOM_THRESHOLD             -> fail "high_rHetHom"
#                         [16] CallRate < CALL_RATE_THRESHOLD          -> fail "low_callrate"
# - Always check:         [11] nSingletons > SINGLETONS_THRESHOLD      -> fail "high_singletons"

awk -v OFS="\t" \
    -v dp_found="$DP_FOUND" \
    -v gt_found="$GT_FOUND" \
    -v cov_thr="$COV_THRESHOLD" \
    -v het_hom_thr="$HET_HOM_THRESHOLD" \
    -v cr_thr="$CALL_RATE_THRESHOLD" \
    -v sing_thr="$SINGLETONS_THRESHOLD" \
    -v contam_thr="$CONTAM_THRESHOLD" '
  BEGIN {
    print "sample","reason(s)","avg_depth","rHetHom","call_rate","nSingletons", "CHARR" > "sample-qc-fails.tsv"
  }

  # Detect PSC header and whether [17] CHARR exists
  /^# PSC[[:space:]]+\[2\]id[[:space:]]+\[3\]sample/ {
    in_psc = 1
    # Check if header line includes CHARR (robust to extra spacing)
    if ($0 ~ /\[17\][[:space:]]+CHARR/) has_charr = 1
    next
  }

  # Any other header ends PSC block
  /^#/ { in_psc = 0; next }


  in_psc && $1=="PSC" {
    sample = $3
    avgd   = $10+0
    nsing  = $11+0
    rHH    = ($15=="NA" ? "NA" : $15+0)
    cr     = $16+0
    chval = ($17=="" ? "nan" : $17+0)

    fail=0
    reasons=""


    if (dp_found>0 && avgd < cov_thr) {
      fail=1; reasons = reasons (reasons?";":"") "low_cov"
    }
    if (gt_found>0 && rHH!="NA" && rHH > het_hom_thr) {
      fail=1; reasons = reasons (reasons?";":"") "high_rHetHom"
    }
    if (gt_found>0 && cr < cr_thr) {
      fail=1; reasons = reasons (reasons?";":"") "low_callrate"
    }
    if (gt_found>0 && nsing > sing_thr) {
      fail=1; reasons = reasons (reasons?";":"") "high_singletons"
    }

    # CHARR if present (col 17), treat NA/-nan as not valid
    charr_raw = (has_charr && NF>=17 ? $17 : "NA")
    charr_valid = (charr_raw!="NA" && charr_raw!="-nan")
    if (charr_valid) charr = charr_raw+0

    # --- Contamination check via CHARR ---
    # Only apply if CHARR column exists and is numeric
    if (has_charr && charr_valid && charr > contam_thr) {
      fail=1; reasons = reasons (reasons?";":"") "high_contam"
    }


    if (fail) {
      print sample, reasons, avgd, rHH, cr, nsing, chval >> "sample-qc-fails.tsv"
      bad[sample]=1
    } else {
      good[sample]=1
    }
    next
  }

  END {
    # emit PASS list
    for (s in good) if (!(s in bad)) print s > "sample-keep.txt"
  }
' sample-qc-stats.txt

echo "QC filtering done." >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
echo "- Failing samples: sample-qc-fails.tsv" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
echo "- Passing samples: sample-keep.txt" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"

# Filter VCF to keep only passing samples
if [[ -s sample-keep.txt ]]; then
    echo "Filtering VCF to keep passing samples..." >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
    bcftools view -S sample-keep.txt -Oz -o "$OUTPUT_VCF" "$INPUT_VCF"
    bcftools index -t "$OUTPUT_VCF"
    echo "✓ Filtered VCF created: $OUTPUT_VCF" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
else
    echo "WARNING: No samples passed QC filters!" >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
    # Create empty VCF with just header
    bcftools view -h "$INPUT_VCF" | bgzip > "$OUTPUT_VCF"
    bcftools index -t "$OUTPUT_VCF"
fi

# add final stats
failed_samples=$(wc -l < sample-qc-fails.tsv) # tsv contains header
nsamples_qc=$((nsamples_orig - failed_samples + 1)) # add one for the 1 substracted because of the header
ngenotypes=$(( nsamples_qc * nvariants ))
{
echo ""
echo " === STATISTICS AFTER QC === "
echo "Variant Number: $nvariants"
echo "Sample Number: $nsamples_qc"
echo "Genotype Number: $ngenotypes" 
} >> "/home/mireia//GitHub/bcftools-pipeline/AF_bcftools_pipeline/output/subset-1000_AFs-females.log"
# Clean up temporary files
rm -rf "$TMP"
