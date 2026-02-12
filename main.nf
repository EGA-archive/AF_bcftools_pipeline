// Ensure the output dir exists (Groovy)
file(params.output).mkdirs()

// =====================
// Process: INDEX_VCF
// Purpose: Ensure each input VCF has a .tbi index and write a header + quick stats to a shared log file.
// Inputs:
//   - vcf        : the VCF (.vcf.gz) to index
//   - has_index  : boolean flag saying whether an index already exists
//   - vcf_idx    : path to the existing index (not used explicitly here; see note below)
//   - log_file   : path (as a value) where messages are appended/written
// Outputs:
//   - vcf        : the same input VCF path
//   - ${vcf}.tbi : expected tabix index filename
//   - log_file   : same value passed through for downstream appends
// Notes:
//   - When has_index=true this script assumes the index is already named "${vcf}.tbi" (or otherwise present in the workdir).


process INDEX_VCF {
    tag "$vcf"

    input:
    tuple path(vcf), val(has_index), path(vcf_idx), val(log_file)

    output:
    tuple path(vcf), path("${vcf}.tbi"), val(log_file)

    debug true

    script:
    if (has_index) {
        """
        touch "${log_file}"
        nvariant="\$(bcftools index -n ${vcf})"
        nsample="\$(bcftools query -l ${vcf} | wc -l)"
        ngenotypes="\$(( nsample * nvariant ))"

        {
            echo " === EGA BCFTools Pipeline ==="
            echo "Developed by: Mireia Marin Ginestar (mireia.marin@crg.eu)"
            echo "version 3.0.0"
            echo ""
            echo "✓ Index already exists for ${vcf}"
            echo ""
            echo " === ORIGINAL STATISTICS === "
            echo "Variant Number: \$nvariant"
            echo "Sample Number: \$nsample"
            echo "Genotype Number: \$ngenotypes"
            echo ""
        } > "${log_file}"
        """
    } else {
        """
        set -euo pipefail
        touch "${log_file}"
        echo "✗ No index found for ${vcf} — creating..." >> "${log_file}"
        tabix -p vcf "${vcf}"
        nvariant="\$(bcftools index -n ${vcf})"
        nsample="\$(bcftools query -l ${vcf} | wc -l)"
        ngenotypes="\$(( nsample * nvariant ))"
        {
          echo " === EGA BCFTools Pipeline ==="
          echo "Developed by: Mireia Marin Ginestar (mireia.marin@crg.eu)"
          echo "version 3.0.0"
          echo ""
          echo "✓ Index created for ${vcf}"
          echo ""
          echo " === ORIGINAL STATISTICS === "
          echo "Variant Number: \$nvariant"
          echo "Sample Number: \$nsample"
          echo "Genotype Number: \$ngenotypes"
        } > "${log_file}"
        """
    }
}

// =====================
// Process: SPLIT_MULTIALLELIC
// Purpose: Split multiallelic variants into separate lines (one ALT per record), reindex, and log summary.
// Inputs:
//   - vcf, tbi  : the (already indexed) VCF and its index
//   - log_file  : log path (value) to append progress
// Outputs:
//   - <name>_split-multiallelic.vcf.gz      : split VCF (bgzipped)
//   - <name>_split-multiallelic.vcf.gz.tbi  : tabix index for the split VCF
//   - log_file  : same log value passed through

process SPLIT_MULTIALLELIC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi), val(log_file)

    output:
    tuple path("${vcf.simpleName}_split-multiallelic.vcf.gz"),
          path("${vcf.simpleName}_split-multiallelic.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail
    {
      echo ""
      echo "=== Splitting multiallelic variants for ${vcf} ==="
    } >> "${log_file}"

    # bcftools norm -m -any: split multi-allelic records into multiple lines (one ALT allele per line)
    bcftools norm -m -any "${vcf}" -Oz -o "${vcf.simpleName}_split-multiallelic.vcf.gz"
    tabix -p vcf "${vcf.simpleName}_split-multiallelic.vcf.gz"

    {
      echo "✓ Split + index done"
      echo "Output: ${vcf.simpleName}_split-multiallelic.vcf.gz"
      echo "Variants now: \$(bcftools index -n "${vcf.simpleName}_split-multiallelic.vcf.gz")"
    } >> "${log_file}"
    """
}

// =====================
// Process: GENOTYPE_QC
// Goal: Mask low‑quality per‑genotype calls to missing (./.) using FORMAT-based rules.
// Inputs:
//   - vcf, tbi : indexed input VCF
//   - log_file : path (value) to append progress & decisions
// Outputs:
//   - <name>-GTmasked.vcf.gz(.tbi) : same variants, but per‑genotype GT set to ./.
//     when any QC rule fails
// Notes:
//   - Uses bcftools +setGT plugin with a per‑genotype expression (-e).
//   - Rules are auto-enabled only if the corresponding FORMAT field exists in the VCF header.


process GENOTYPE_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi), val(log_file)

    output:
    tuple path("${vcf.simpleName}-GTmasked.vcf.gz"),
          path("${vcf.simpleName}-GTmasked.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail

    VCF_IN="${vcf}"
    VCF_OUT="${vcf.simpleName}-GTmasked.vcf.gz"
    OP="||"

    # Define candidate per‑genotype rules (only applied if the FORMAT tag exists)
    #  - GQ: genotype quality below threshold
    #  - DP: per‑genotype depth below threshold
    #  - AD: allele balance (ALT fraction) below threshold when total AD>0
    #        ALT fraction = AD[1] / (AD[0] + AD[1])
    
    declare -A gt_conditions=(
      [GQ]='FMT/GQ < ${params.qc.genotype.gq_threshold}'
      [DP]='FMT/DP < ${params.qc.genotype.dp_threshold}'
      [AD]='GT="het" && (FMT/AD[*:0]+FMT/AD[*:1])>0 && (FMT/AD[*:1]/(FMT/AD[*:0]+FMT/AD[*:1])) < ${params.qc.genotype.ab_ratio_threshold}'
    )

    {
      echo ""
      echo "=== GENOTYPE_QC on: \$VCF_IN ==="
    } >> "${log_file}"
    
    # Build the final mask expression only with tags present in the header
    gt_expr_parts=()
    for tag in "\${!gt_conditions[@]}"; do
      if bcftools view -h "\$VCF_IN" | grep -q "^##FORMAT=<ID=\${tag},"; then
        echo "✓ \${tag} (FORMAT) found — adding rule" >> "${log_file}"
        gt_expr_parts+=("\${gt_conditions[\$tag]}")
      else
        echo "x \${tag} (FORMAT) not found — skipping" >> "${log_file}"
      fi
    done

    # If no rules apply (none of the tags are present), pass-through the file unchanged

    if (( \${#gt_expr_parts[@]} == 0 )); then
      {
        echo "x No FORMAT-based rules available; no masking performed."
      } >> "${log_file}"
      cp -f "\$VCF_IN" "\$VCF_OUT"
      cp -f "\${VCF_IN}.tbi" "\${VCF_OUT}.tbi" 2>/dev/null || tabix -p vcf "\$VCF_OUT"
      exit 0
    fi

    # Combine rules with OR (mask if ANY condition is true)
    gt_expr="\${gt_expr_parts[0]}"
    for cond in "\${gt_expr_parts[@]:1}"; do gt_expr+=" \$OP \$cond"; done

    {
      echo "Final per-genotype mask expression:"
      echo "\$gt_expr"
      echo "Running bcftools +setGT (mask to ./.)"
    } >> "${log_file}"

    # Mask failing genotypes to ./.
    # -t q : expression applies per-sample/per-genotype (FORMAT context)
    # -n . : set GT to missing when -e expr evaluates to true
    # -e   : mask expression (constructed above)

    bcftools +setGT "\$VCF_IN" -Oz -o "\$VCF_OUT" -- -t q -n . -e "\$gt_expr"
    tabix -p vcf "\$VCF_OUT"

    {
      echo "✓ Genotype masking complete."
      echo "Output: \$VCF_OUT"
    } >> "${log_file}"
    """
}

// =====================
// Process: VARIANT_QC (hard filter)
// Goal: REMOVE non-passing variants (only PASSing sites remain in output).
// Strategy: Build a site-level boolean expression; `bcftools view -e <expr>`
//           EXCLUDES variants where the expression is TRUE.
// Inputs:
//   - vcf, tbi  : indexed input VCF
//   - log_file  : path (value) to append progress & decisions
// Outputs:
//   - <name>-variantQC.vcf.gz(.tbi) : VCF containing only variants that pass all active rules
// Notes:
//   - An INFO rule is added only if that tag exists in the header (robust to missing annotations).
//   - By default we OR the rules (fail if ANY rule is true). Switch OP to "&&" to require ALL.

process VARIANT_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi), val(log_file)

    output:
    tuple path("${vcf.simpleName}-variantQC.vcf.gz"),
          path("${vcf.simpleName}-variantQC.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail

    VCF_IN="${vcf}"
    VCF_OUT="${vcf.simpleName}-variantQC.vcf.gz"
    before_count=\$(bcftools index -n "\$VCF_IN") # original number of variants
    OP="||" # Combine conditions with OR: EXCLUDE variant if any condition is true

    # Candidate site-level rules (only enabled if the INFO tag is present)
    # QD  : Qual/Depth below threshold
    # DP  : Site depth below threshold (requires INFO/DP to be present)
    # MQ  : Mapping quality below threshold
    # FS  : Phred-scaled strand bias above threshold
    # ReadPosRankSum : Read position bias less than threshold
    declare -A site_conditions=(
      [QD]='INFO/QD < ${params.qc.variant.qd_threshold}'
      [DP]='INFO/DP < ${params.qc.variant.dp_threshold}'
      [MQ]='INFO/MQ < ${params.qc.variant.mq_threshold}'
      [FS]='INFO/FS > ${params.qc.variant.fs_threshold}'
      [ReadPosRankSum]='INFO/ReadPosRankSum < ${params.qc.variant.read_pos_rank_sum_threshold}'
    )

    has_info() { bcftools view -h "\$1" | grep -q "^##INFO=<ID=\$2,"; }

    {
      echo ""
      echo "=== VARIANT_QC on: \$VCF_IN ==="
    } >> "${log_file}"

    expr_parts=("QUAL < ${params.qc.variant.qual_threshold}") # Always include QUAL threshold (QUAL is a core VCF field, not in INFO)
    echo "✓ QUAL — adding: QUAL < ${params.qc.variant.qual_threshold}" >> "${log_file}"

    # Add INFO-based rules only if that tag exists
    for tag in "\${!site_conditions[@]}"; do
      if has_info "\$VCF_IN" "\$tag"; then
        echo "✓ \$tag (INFO) found — adding: \${site_conditions[\$tag]}" >> "${log_file}"
        expr_parts+=("\${site_conditions[\$tag]}")
      else
        echo "x \$tag (INFO) not found — skipping" >> "${log_file}"
      fi
    done

    expr_str="\${expr_parts[0]}"
    for cond in "\${expr_parts[@]:1}"; do expr_str+=" \$OP \$cond"; done

    {
      echo ""
      echo "Final filter expression:"
      echo "\$expr_str"
    } >> "${log_file}"

    # HARD FILTER: remove variants where expr is TRUE
    bcftools view -e "\$expr_str" "\$VCF_IN" -Oz -o "\$VCF_OUT"
    tabix -p vcf "\$VCF_OUT"
    after_count=\$(bcftools index -n "\$VCF_OUT")

    {
      echo ""
      echo "✓ Filtering complete. Output: \$VCF_OUT"
      echo "Removed: \$(( before_count - after_count ))"
    } >> "${log_file}"
    """
}

// =====================
// Process: SAMPLE_QC
// Goal: Identify and remove samples that fail per-sample QC metrics, then reindex the VCF.
// Metrics:
//   1) Mean coverage (FORMAT/DP)               -> below threshold
//   2) Call rate (FORMAT/GT present)           -> below threshold
//   3) Het/Hom ratio (from GT on SNVs)         -> above threshold (suggesting potential issues)
//   4) Singletons count                        -> above threshold
//   5) Contamination (sceVCF, requires FORMAT/AD and sceVCF available) -> above threshold
//
// Inputs:
//   - vcf, tbi         : indexed input VCF
//   - log_file         : path (value) to append progress & decisions
//   - sceVCF_path (val): "", or a directory containing `sceVCF`, or a full path to the `sceVCF` binary
//   - seq_type   (val) : "WES" or "WGS" (to pick thresholds)
// Outputs:
//   - <name>-sampleQC.vcf.gz(.tbi) with failing samples removed
//
// Notes:
//   - Uses only PASS SNVs for per-sample stats via SITE_SUBSET_CMD.
//   - If no samples are flagged, we copy the input to the output (safer than renaming).
//   - The contamination step is skipped if AD is missing or sceVCF is unavailable/unspecified.

process SAMPLE_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi), val(log_file)
    val sceVCF_path
    val seq_type

    output:
    tuple path("${vcf.simpleName}-sampleQC.vcf.gz"),
          path("${vcf.simpleName}-sampleQC.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    # Fail fast + propagate errors in pipelines
    set -euo pipefail

    INPUT_VCF="${vcf}"
    OUTPUT_VCF="${vcf.simpleName}-sampleQC.vcf.gz"

    # Workspace for intermediate per-sample metrics
    TMP="qc_tmp"; mkdir -p "\$TMP"

    {
      echo ""
      echo "=== SAMPLE_QC on: \$INPUT_VCF ==="
    } >> "${log_file}"

    # -----------------------------
    # Choose thresholds by seq type
    # -----------------------------
    if [[ "${seq_type}" == "WGS" ]]; then
      COV_THRESHOLD=${params.qc.sample.wgs.coverage_threshold}
      HET_HOM_THRESHOLD=${params.qc.sample.wgs.het_hom_threshold}
      CALL_RATE_THRESHOLD=${params.qc.sample.wgs.call_rate_threshold}
      SINGLETONS_THRESHOLD=${params.qc.sample.wgs.singletons_threshold}
      CONTAM_THRESHOLD=${params.qc.sample.wgs.contamination_threshold}
      echo "✓ Seq type: WGS thresholds applied" >> "${log_file}"
    else
      COV_THRESHOLD=${params.qc.sample.wes.coverage_threshold}
      HET_HOM_THRESHOLD=${params.qc.sample.wes.het_hom_threshold}
      CALL_RATE_THRESHOLD=${params.qc.sample.wes.call_rate_threshold}
      SINGLETONS_THRESHOLD=${params.qc.sample.wes.singletons_threshold}
      CONTAM_THRESHOLD=${params.qc.sample.wes.contamination_threshold}
      echo "✓ Seq type: WES thresholds applied" >> "${log_file}"
    fi

    # Generate bcftools stats for the input VCF
    echo "Generating bcftools stats..." >> "${log_file}"
    bcftools query -l "\$INPUT_VCF" > list-of-samples.txt
    bcftools stats -S list-of-samples.txt "\$INPUT_VCF" > bcftools-stats.txt

    # --- detect FORMAT presence in the VCF header ---
    DP_FOUND=\$(bcftools view -h "\$INPUT_VCF" | grep -c '^##FORMAT=<ID=DP,' || true)
    GT_FOUND=\$(bcftools view -h "\$INPUT_VCF" | grep -c '^##FORMAT=<ID=GT,' || true)


    if [[ "\$DP_FOUND" -gt 0 ]]; then
      echo "✓ DP format found - Performing coverage filtering" >> "${log_file}"
    else
      echo "x DP format NOT found - Coverage filtering NOT performed" >> "${log_file}"
    fi
    if [[ "\$GT_FOUND" -gt 0 ]]; then
      echo "✓ GT format found - Performing call rate, het/hom ratio and singleton filtering" >> "${log_file}"
    else
      echo "x GT format NOT found - Call rate, het/hom ratio and singleton filtering NOT performed" >> "${log_file}"
    fi

    # 1) Extract SN + PSC sections from bcftools stats
    awk '
      printing && /^#/ { printing=0 }
      /^# SN/ { printing=1 }
      /^# PSC[[:space:]]+\\[2\\]id[[:space:]]+\\[3\\]sample/ { printing=1 }
      printing
    ' bcftools-stats.txt > bcftools-stats_min.txt

    # 2) Read number of samples and number of records from SN section
    read nsamples_orig nvariants < <(
      awk '
        \$3=="number" && \$4=="of" && \$5=="samples:" { ns=\$6 }
        \$3=="number" && \$4=="of" && \$5=="records:" { nr=\$6 }
        END { print ns, nr }
      ' bcftools-stats_min.txt
    )

    # 3) Compute total genotypes
    ngenotypes=\$(( nsamples_orig * nvariants ))

    # 4) Add rHetHom [15] and CallRate [16] to PSC section
    awk -v OFS="\\t" -v ngenotypes="\$ngenotypes" '
      # Extend PSC header
      /^# PSC[[:space:]]+\\[2\\]id[[:space:]]+\\[3\\]sample[[:space:]]+\\[4\\]nRefHom/ {
        in_psc=1
        print \$0, "[15] rHetHom", "[16] CallRate"
        next
      }

      # Any other header stops PSC mode
      /^#/ { in_psc=0; print; next }

      # PSC rows: compute new metrics
      in_psc && \$1=="PSC" {
        r = (\$5+0)==0 ? "NA" : sprintf("%.6f", \$6/\$5)        # [15]
        cr = sprintf("%.6f", (1 - (\$14/ngenotypes)))         # [16], using [14]=nMissing
        print \$0, r, cr
        next
      }

      { print }
    ' bcftools-stats_min.txt > sample-qc-stats.txt

    # 5) Contamination check 
    CHARR_TSV="sceVCF-results.tsv"
    : > "\$CHARR_TSV"  # create empty file by default

    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=AD,"; then
      if [[ -n "${sceVCF_path}" ]]; then
        SCE_CMD=""
        if [[ -x "${sceVCF_path}" ]]; then
          SCE_CMD="${sceVCF_path}"
        elif [[ -x "${sceVCF_path}/sceVCF" ]]; then
          SCE_CMD="${sceVCF_path}/sceVCF"
        fi

        if [[ -n "\$SCE_CMD" ]]; then
          echo "✓ sceVCF found (\$SCE_CMD) — running contamination check" >> "${log_file}"
          "\$SCE_CMD" -o "\$CHARR_TSV" "\$INPUT_VCF"
          awk -F'\t' -v OFS='\t' '
          FNR==NR {
            # Read sceVC results into array (skip meta)
            if (\$1 ~ /^#/ || \$1 ~ /^##/) next
            charr[\$1] = \$10
            next
          }
          # When we hit the #PSC header, append [17] CHARR
          /^# PSC/ {
            print \$0, "[17] CHARR"
            next
          }
          # For PSC data lines, append CHARR value
          /^PSC/ {
            val = (\$3 in charr ? charr[\$3] : "NA")
            print \$0, val
            next
          }
          # Print all other lines unchanged
          { print }
        ' charr_full.tsv sample-qc-stats.txt > sample-qc-stats.txt.tmp && mv sample-qc-stats.txt.tmp sample-qc-stats.txt

        else
          echo "x sceVCF not found or not executable at: ${sceVCF_path} — skipping" >> "${log_file}"
        fi
      else
        echo "x Contamination check not running (sceVCF_path empty)" >> "${log_file}"
      fi
    else
      echo "x AD (FORMAT) not found — skipping contamination check (required for sceVCF)" >> "${log_file}"
    fi


    # --- now filter samples from the PSC section according to your rules ---
    # Rules:
    # - If DP present:        [10] average depth < COV_THRESHOLD           -> fail "low_cov"
    # - If GT present:        [15] rHetHom > HET_HOM_THRESHOLD             -> fail "high_rHetHom"
    #                         [16] CallRate < CALL_RATE_THRESHOLD          -> fail "low_callrate"
    # - Always check:         [11] nSingletons > SINGLETONS_THRESHOLD      -> fail "high_singletons"

    awk -v OFS="\\t" \\
        -v dp_found="\$DP_FOUND" \\
        -v gt_found="\$GT_FOUND" \\
        -v cov_thr="\$COV_THRESHOLD" \\
        -v het_hom_thr="\$HET_HOM_THRESHOLD" \\
        -v cr_thr="\$CALL_RATE_THRESHOLD" \\
        -v sing_thr="\$SINGLETONS_THRESHOLD" \\
        -v contam_thr="\$CONTAM_THRESHOLD" '
      BEGIN {
        print "sample","reason(s)","avg_depth","rHetHom","call_rate","nSingletons", "CHARR" > "sample-qc-fails.tsv"
      }

      # Detect PSC header and whether [17] CHARR exists
      /^# PSC[[:space:]]+\\[2\\]id[[:space:]]+\\[3\\]sample/ {
        in_psc = 1
        # Check if header line includes CHARR (robust to extra spacing)
        if (\$0 ~ /\\[17\\][[:space:]]+CHARR/) has_charr = 1
        next
      }

      # Any other header ends PSC block
      /^#/ { in_psc = 0; next }

    
      in_psc && \$1=="PSC" {
        sample = \$3
        avgd   = \$10+0
        nsing  = \$11+0
        rHH    = (\$15=="NA" ? "NA" : \$15+0)
        cr     = \$16+0
        chval = (\$17=="" ? "nan" : \$17+0)
        
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
        charr_raw = (has_charr && NF>=17 ? \$17 : "NA")
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

    echo "QC filtering done." >> "${log_file}"
    echo "- Failing samples: sample-qc-fails.tsv" >> "${log_file}"
    echo "- Passing samples: sample-keep.txt" >> "${log_file}"

    # Filter VCF to keep only passing samples
    if [[ -s sample-keep.txt ]]; then
        echo "Filtering VCF to keep passing samples..." >> "${log_file}"
        bcftools view -S sample-keep.txt -Oz -o "\$OUTPUT_VCF" "\$INPUT_VCF"
        bcftools index -t "\$OUTPUT_VCF"
        echo "✓ Filtered VCF created: \$OUTPUT_VCF" >> "${log_file}"
    else
        echo "WARNING: No samples passed QC filters!" >> "${log_file}"
        # Create empty VCF with just header
        bcftools view -h "\$INPUT_VCF" | bgzip > "\$OUTPUT_VCF"
        bcftools index -t "\$OUTPUT_VCF"
    fi

    # add final stats
    failed_samples=\$(wc -l < sample-qc-fails.tsv) # tsv contains header
    nsamples_qc=\$((nsamples_orig - failed_samples + 1)) # add one for the 1 substracted because of the header
    ngenotypes=\$(( nsamples_qc * nvariants ))
    {
    echo ""
    echo " === STATISTICS AFTER QC === "
    echo "Variant Number: \$nvariants"
    echo "Sample Number: \$nsamples_qc"
    echo "Genotype Number: \$ngenotypes" 
    } >> "${log_file}"
    # Clean up temporary files
    rm -rf "\$TMP"
    
    """
}

// =====================
// Process: FIX_PLOIDY
// Goal: Correct ploidy for haploid regions (X and Y chromosomes) based on sample sex
//       to ensure accurate hemizygous genotype calls and allele frequency calculations.
// Inputs:
//   - vcf, tbi      : indexed input VCF
//   - log_file (val): log path to append progress
//   - metadata_csv  : CSV with header; columns: SAMPLE, SEX (M/F or 1/2), ANCESTRY
// Outputs:
//   - <name>-ploidy_fixed.vcf.gz(.tbi): VCF with corrected ploidy for X/Y chromosomes
//
// Notes:
//   - Creates a gender.txt file from metadata (sample<TAB>sex)
//   - Uses bcftools +fixploidy with default ploidy rules for human X/Y chromosomes
//   - Depending on the reference_genome provided by the user, coordinates from GRCh37 or GRCh38 will be used. Check /assets to for further details

process FIX_PLOIDY {
    tag "$vcf"
    
    input:
    tuple path(vcf), path(tbi), val(log_file)
    path metadata_csv
    path ploidy_rules

    output:
    tuple path("${vcf.simpleName}-ploidy_fixed.vcf.gz"),
          path("${vcf.simpleName}-ploidy_fixed.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail

    INPUT_VCF="${vcf}"
    OUTPUT_VCF="${vcf.simpleName}-ploidy_fixed.vcf.gz"

    {
      echo ""
      echo "=== FIX_PLOIDY: Creating gender.txt from metadata ==="
    } >> "${log_file}"

    # Create gender.txt file: sample<TAB>sex
    # Convert M/F or 1/2 to M/F format for bcftools +fixploidy
    awk -F, 'NR>1 {
      sample = \$1
      sex_code = \$2
      
      # Remove whitespace
      gsub(/^[[:space:]]+|[[:space:]]+\$/, "", sex_code)
      gsub(/^[[:space:]]+|[[:space:]]+\$/, "", sample)
      
      # Convert to M/F format
      sex = ""
      if (sex_code == "M" || sex_code == "1") sex = "M"
      else if (sex_code == "F" || sex_code == "2") sex = "F"
      
      # Output: sample<TAB>sex
      if (sex != "") print sample "\\t" sex
    }' "${metadata_csv}" > gender.txt

    echo "✓ Gender file created" >> "${log_file}"
    
    # Show sample count
    n_samples=\$(wc -l < gender.txt)
    n_males=\$(awk '\$2=="M"' gender.txt | wc -l)
    n_females=\$(awk '\$2=="F"' gender.txt | wc -l)
    
    {
      echo "  Total samples in metadata file: \$n_samples"
      echo "  Males: \$n_males"
      echo "  Females: \$n_females"
      echo "  Only samples that passed the QC will be used for the groupings "
      echo ""
      echo "Example gender assignments:"
      head -n 5 gender.txt
      echo ""
      echo "=== Fixing ploidy for X/Y chromosomes ==="
    } >> "${log_file}"

    # Fix ploidy 
    # If GRCh37, ploidy_cmd points to the provided rules for GRCh37
    # If GRCh38, ploidy_cmd points to the provided rules for GRCh38 
    bcftools +fixploidy "\${INPUT_VCF}" -Oz -o "\$OUTPUT_VCF" -- -s gender.txt -p ${ploidy_rules}

    # Index the output
    tabix -p vcf "\$OUTPUT_VCF"

    {
      echo "✓ Ploidy correction complete"
      echo "  Output: \$OUTPUT_VCF"
      echo "  Males will have haploid genotypes (0 or 1) on X and Y"
      echo "  Females will have diploid genotypes (0/0, 0/1, 1/1) on X"
    } >> "${log_file}"

   
    """
}


// =====================
// Process: ADD_AF
// Goal: Recalculate allele-frequency tags (AF/AC/AN/etc.), stratified by:
//       1) Overall (all samples)
//       2) By sex (MALE, FEMALE)
//       3) By ancestry (EUR, AFR, etc.)
//       4) By ancestry+sex combinations (EUR_MALE, EUR_FEMALE, etc.)
//       Then drop genotypes/FORMAT to emit an INFO-only VCF.
// Inputs:
//   - vcf,tbi       : indexed input VCF
//   - log_file (val): log path to append progress
//   - metadata_csv  : CSV with header; columns used here: SAMPLE, SEX(1/2), ANCESTRY
// Outputs:
//   - <name>-AF_recalc.vcf.gz(.tbi): INFO-only VCF with stratified AF/AC/AN tags
//
// Notes:
//   - bcftools +fill-tags can compute tags per population when given a groups file
//     (sample<TAB>comma-separated-groups).
//   - We create groups for: SEX alone, ANCESTRY alone, and ANCESTRY_SEX combinations

process ADD_AF {
    tag "$vcf"
    publishDir "${params.output}", mode: 'move', pattern: "*-AF_recalc.vcf.gz*"
    
    input:
    tuple path(vcf), path(tbi), val(log_file)
    path metadata_csv

    output:
    tuple path("${vcf.simpleName}-AF_recalc.vcf.gz"),
          path("${vcf.simpleName}-AF_recalc.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail

    INPUT_VCF="${vcf}"
    OUTPUT_VCF="${vcf.simpleName}-AF_recalc.vcf.gz"

    {
      echo ""
      echo "=== AF recalculating: Creating groups.txt from metadata ==="
    } >> "${log_file}"

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
      sample = \$1
      sex_code = \$2
      ancestry = \$3
      
      # Remove any whitespace
      gsub(/^[[:space:]]+|[[:space:]]+\$/, "", sex_code)
      gsub(/^[[:space:]]+|[[:space:]]+\$/, "", ancestry)
      
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
      print sample "\\t" groups
    }' "${metadata_csv}" > groups.txt

    echo "✓ Groups file created with ancestry, sex, and ancestry+sex combinations" >> "${log_file}"
    
    # Show a few example lines for verification
    {
      echo "Example group assignments:"
      head -n 5 groups.txt
    } >> "${log_file}"

    {
      echo ""
      echo "=== Adding stratified allele frequencies to VCF (dropping all FORMAT/GT columns) ==="
    } >> "${log_file}"

    # Recalculate AF/AC/AN (and other tags) across:
    #   - All samples (default)
    #   - Per-group: MALE, FEMALE, EUR, AFR, EUR_MALE, EUR_FEMALE, etc.
    # Then drop genotypes (-G) and strip any FORMAT header remnants (-x FORMAT).
    
    bcftools +fill-tags "\${INPUT_VCF}" -Ou -- -S groups.txt \
      | bcftools view -G -Ou \
      | bcftools annotate -x FORMAT \
      -Oz -o "\$OUTPUT_VCF"

    # Index the INFO-only VCF
    tabix -p vcf "\$OUTPUT_VCF"

    {
      echo "✓ AF annotation complete with stratified groups"
      echo "  Output: \$OUTPUT_VCF"
      echo ""
      echo "AF tags created for:"
      echo "  - Overall: AF, AC, AN"
      echo "  - By sex: AF_MALE, AF_FEMALE, AC_MALE, AC_FEMALE, AN_MALE, AN_FEMALE"
      echo "  - By ancestry: AF_EUR, AF_AFR, etc."
      echo "  - By ancestry+sex: AF_EUR_MALE, AF_EUR_FEMALE, etc."
    } >> "${log_file}"
    """
}



workflow {
    // Channel: all .vcf.gz in params.input 
    vcf_ch = Channel.fromPath("${params.input}/*.vcf.gz", checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            def log_file = "${params.output}/${vcf.simpleName}.log"
            if (tbi.exists()) {
                tuple(vcf, true, tbi, log_file)
            } else {
                tuple(vcf, false, file("NO_FILE"), log_file)
            }
        }

    // Run INDEX_VCF to create indexed VCF files
    indexed = INDEX_VCF(vcf_ch)

    // Run SPLIT_MULTIALLELIC using the indexed VCF files
    split_multiallelic = SPLIT_MULTIALLELIC(indexed)

    genotype_qc = GENOTYPE_QC(split_multiallelic)

    // Run VARIANT_QC using the indexed VCF with the multiallelic variants splitted
    variant_qc = VARIANT_QC(genotype_qc)

    sample_qc = SAMPLE_QC(variant_qc, params.sceVCF_path, params.seq_type)

    // Create metadata channel if provided
    if (params.metadata_csv) {
        metadata_ch = Channel.fromPath(params.metadata_csv, checkIfExists: true)
    } else {
        error "ERROR: metadata_csv parameter is required for ADD_AF process"
    }
    
    // Select the correct ploidy file based on the config
    def ploidy_file = params.reference_genome == 'GRCh38' ? 
        file("${projectDir}/assets/ploidy_grch38.txt") : 
        file("${projectDir}/assets/ploidy_grch37.txt")

    // Fix ploidy for X/Y chromosomes before AF calculation
    ploidy_fixed = FIX_PLOIDY(sample_qc, metadata_ch, ploidy_file)
    
    // Add AF annotations
    af_annotated = ADD_AF(ploidy_fixed, metadata_ch)
}