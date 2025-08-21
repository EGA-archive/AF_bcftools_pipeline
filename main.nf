nextflow.enable.dsl=2

// Ensure the output dir exists (Groovy)
file(params.output_stats).mkdirs()

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
        {
          echo " === EGA BCFTools Pipeline ==="
          echo "Developed by: Mireia Marin Ginestar (mireia.marin@crg.eu)"
          echo "version 1.0.0"
          echo ""
          echo "✓ Index already exists for ${vcf}"
          echo ""
          echo " === ORIGINAL VARIANT AND SAMPLE NUMBER === "
          echo "Variant Number: \$(bcftools index -n ${vcf})"
          echo "Sample Number: \$(bcftools query -l ${vcf} | wc -l)"
        } > "${log_file}"
        """
    } else {
        """
        set -euo pipefail

        touch "${log_file}"
        echo "✗ No index found for ${vcf} — creating..." >> "${log_file}"
        tabix -p vcf "${vcf}"
        {
          echo " === EGA BCFTools Pipeline ==="
          echo "Developed by: Mireia Marin Ginestar (mireia.marin@crg.eu)"
          echo "version 1.0.0"
          echo ""
          echo "✓ Index created for ${vcf}"
          echo ""
          echo " === ORIGINAL VARIANT AND SAMPLE NUMBER === "
          echo "Original Variant Number: \$(bcftools index -n ${vcf})"
          echo "Original Sample Number: \$(bcftools query -l ${vcf} | wc -l)"
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
      [AD]='(FMT/AD[0]+FMT/AD[1])>0 && (FMT/AD[1]/(FMT/AD[0]+FMT/AD[1])) < ${params.qc.genotype.ad_ratio_threshold}'
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
      cp -a "\$VCF_IN" "\$VCF_OUT"
      cp -a "\${VCF_IN}.tbi" "\${VCF_OUT}.tbi" 2>/dev/null || tabix -p vcf "\$VCF_OUT"
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

    # Subset to PASS SNVs only for per-sample metrics (faster and cleaner)
    SITE_SUBSET_CMD=(bcftools view --threads ${task.cpus} -f PASS -v snps "\$INPUT_VCF")

    # -----------------------------
    # 1) Mean coverage (FORMAT/DP)
    # -----------------------------
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=DP,"; then
      echo "✓ DP (FORMAT) found — calculating mean coverage" >> "${log_file}"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%DP\\n]' \
      | awk '{sum[\$1]+=\$2; n[\$1]++} END{for(s in sum){if(n[s]>0) printf "%s\\t%.6f\\n",s,sum[s]/n[s]}}' \
      > "\$TMP/mean_dp.txt"
      awk -v thr="\$COV_THRESHOLD" '\$2 < thr {print \$1}' "\$TMP/mean_dp.txt" > "\$TMP/low_cov_samples.txt"
    else
      echo "x DP (FORMAT) not found — skipping mean coverage" >> "${log_file}"
      : > "\$TMP/mean_dp.txt"; : > "\$TMP/low_cov_samples.txt"
    fi

    # -----------------------------
    # 2) Call rate (FORMAT/GT)
    # -----------------------------
    GT_FOUND=\$(bcftools view -h "\$INPUT_VCF" | grep -c '^##FORMAT=<ID=GT,')
    if [[ "\$GT_FOUND" -gt 0 ]]; then
      echo "✓ GT (FORMAT) found — calculating call rate" >> "${log_file}"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%GT\\n]' \
      | awk '{tot[\$1]++; if(\$2=="./."||\$2==".|.") miss[\$1]++} END{for(s in tot){cr=1-((miss[s]+0)/tot[s]); printf "%s\\t%.6f\\n",s,cr}}' \
      > "\$TMP/call_rate.txt"
      awk -v thr="\$CALL_RATE_THRESHOLD" '\$2 < thr {print \$1}' "\$TMP/call_rate.txt" > "\$TMP/low_call_rate.txt"
    else
      echo "x GT (FORMAT) not found — skipping call rate" >> "${log_file}"
      : > "\$TMP/call_rate.txt"; : > "\$TMP/low_call_rate.txt"
    fi

    # -----------------------------
    # 3) Het/Hom ratio (from GT)
    #    het / hom_alt over PASS SNVs
    # -----------------------------
    if [[ "\$GT_FOUND" -gt 0 ]]; then
      echo "✓ GT (FORMAT) found — calculating Het/Hom ratio" >> "${log_file}"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%GT\\n]' \
      | awk '{
          g=\$2
          if(g ~ /[01][\\/|][01]/) {
            split(g,a,/[/|]/)
            if(a[1]!=a[2]) het[\$1]++
            else if(a[1]=="1") hom[\$1]++
          }
        }
        END{
          # het/hom_alt, handling zeros and pure-het edge cases
          PROCINFO["sorted_in"]="@ind_str_asc"
          for(s in het){
            if(hom[s]>0) printf "%s\\t%.6f\\n",s,het[s]/hom[s];
            else if(het[s]>0) printf "%s\\tinf\\n",s;
            else printf "%s\\t0\\n",s;
          }
          for(s in hom) if(!(s in het)) printf "%s\\t0\\n",s;
        }' \
      > "\$TMP/het_hom.txt"
      awk -v thr="\$HET_HOM_THRESHOLD" '\$2=="inf" || \$2+0 > thr {print \$1}' "\$TMP/het_hom.txt" > "\$TMP/bad_het_hom.txt"
    else
      echo "x GT (FORMAT) not found — skipping Het/Hom ratio" >> "${log_file}"
      : > "\$TMP/het_hom.txt"; : > "\$TMP/bad_het_hom.txt"
    fi

    # -----------------------------
    # 4) Singletons per sample
    #    Count heterozygous occurrences at AC==1 sites.
    # -----------------------------
    echo "✓ Counting singletons" >> "${log_file}"
    "\${SITE_SUBSET_CMD[@]}" \
    | bcftools view -i 'AC==1' -Ou \
    | bcftools query -f '[%SAMPLE\\t%GT\\n]' \
    | awk '{
        g=\$2
        if(g!="./." && g!=".|.") {
          if(g=="0/1" || g=="0|1" || g=="1/0" || g=="1|0") count[\$1]++
        }
      }
      END{for(s in count) printf "%s\\t%d\\n",s,count[s]}' \
    > "\$TMP/singletons.txt"
    awk -v thr="\$SINGLETONS_THRESHOLD" '\$2 > thr {print \$1}' "\$TMP/singletons.txt" > "\$TMP/high_singletons.txt"

    # -----------------------------
    # 5) Contamination via sceVCF (optional)
    #    Requirements: FORMAT/AD exists AND sceVCF is resolvable.
    #    We accept:
    #       - sceVCF_path == ""                     -> skip
    #       - sceVCF_path == "/path/to/sceVCF"      -> call that binary
    #       - sceVCF_path == "/path/to/dir"         -> use "/path/to/dir/sceVCF" if executable
    # -----------------------------
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=AD,"; then
      if [[ -n "${sceVCF_path}" ]]; then
        SCE_CMD=""
        if [[ -x "${sceVCF_path}" ]]; then
          # Full path to the binary
          SCE_CMD="${sceVCF_path}"
        elif [[ -x "${sceVCF_path}/sceVCF" ]]; then
          # Directory containing the binary
          SCE_CMD="${sceVCF_path}/sceVCF"
        fi

        if [[ -n "\$SCE_CMD" ]]; then
          echo "✓ sceVCF found (\$SCE_CMD) — running contamination check" >> "${log_file}"
          "\$SCE_CMD" -o "\$TMP/charr_full.tsv" "\$INPUT_VCF"
          # Expecting a TSV with at least: SAMPLE<TAB>CONTAM_VALUE on column 2
          awk -v thr="\$CONTAM_THRESHOLD" '\$2 > thr {print \$1}' "\$TMP/charr_full.tsv" > "\$TMP/high_contam.txt"
        else
          echo "x sceVCF not found or not executable at: ${sceVCF_path} — skipping" >> "${log_file}"
          : > "\$TMP/high_contam.txt"
        fi
      else
        echo "x Contamination check not running (sceVCF_path empty)" >> "${log_file}"
        : > "\$TMP/high_contam.txt"
      fi
    else
      echo "x AD (FORMAT) not found — skipping contamination check (required for sceVCF)" >> "${log_file}"
      : > "\$TMP/high_contam.txt"
    fi

    # -------------------------------------
    # Aggregate failing samples across tests
    # -------------------------------------
    cat "\$TMP"/low_cov_samples.txt "\$TMP"/low_call_rate.txt "\$TMP"/bad_het_hom.txt "\$TMP"/high_singletons.txt "\$TMP/high_contam.txt" \
    2>/dev/null | sort -u > "\$TMP/samples_to_remove.txt"

    # Log per-test counts (nice for quick diagnostics)
    for f in mean_dp.txt low_cov_samples.txt call_rate.txt low_call_rate.txt het_hom.txt bad_het_hom.txt singletons.txt high_singletons.txt high_contam.txt; do
      [[ -f "\$TMP/\$f" ]] || continue
      n=\$(wc -l < "\$TMP/\$f"); echo "  - \${f}: \${n} lines" >> "${log_file}"
    done

    # -------------------------------------
    # Remove failing samples (if any)
    # -------------------------------------
    if [[ -s "\$TMP/samples_to_remove.txt" ]]; then
      rm_count=\$(wc -l < "\$TMP/samples_to_remove.txt")
      echo "✗ Removing \${rm_count} samples failing QC" >> "${log_file}"
      bcftools view --threads ${task.cpus} -S ^"\$TMP/samples_to_remove.txt" "\$INPUT_VCF" -Oz -o "\$OUTPUT_VCF"
      tabix -p vcf "\$OUTPUT_VCF"
    else
      {
        echo "✓ No samples flagged — copying input to output"
      } >> "${log_file}"
      # Use copy (not move) to avoid altering the staged input
      cp -a "\$INPUT_VCF" "\$OUTPUT_VCF"
      if [[ -f "\${INPUT_VCF}.tbi" ]]; then
        cp -a "\${INPUT_VCF}.tbi" "\$OUTPUT_VCF.tbi"
      else
        tabix -p vcf "\$OUTPUT_VCF"
      fi
    fi

    {
      echo "✓ SAMPLE_QC complete. Output: \$OUTPUT_VCF"
    } >> "${log_file}"
    """
}

// =====================
// Process: ADD_AF
// Goal: Recalculate allele-frequency tags (AF/AC/AN/etc.), optionally per group,
//       from a CSV metadata file, then drop genotypes/FORMAT to emit an INFO-only VCF.
// Inputs:
//   - vcf,tbi       : indexed input VCF
//   - log_file (val): log path to append progress
//   - metadata_csv  : CSV with header; columns used here: SAMPLE, SEX(1/2), ANCESTRY
// Outputs:
//   - <name>-AF_recalc.vcf.gz(.tbi): INFO-only VCF with (per-group) AF/AC/AN, etc.
//
// Notes:
//   - bcftools +fill-tags can compute tags per population when given a groups file
//     (sample<TAB>comma-separated-groups).
//   - We stream through: +fill-tags → view -G (drop all genotypes) → annotate -x FORMAT (clean header).

process ADD_AF {
    tag "$vcf"

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
    # Here: groups = {SEX, ANCESTRY} if present. SEX in metadata is 1=MALE, 2=FEMALE.
    awk -F, 'NR>1{
      s=\$1
      sex=tolower(\$2)
      anc=\$3
      g=""
      if (sex=="2") g=g"FEMALE"
      else if (sex=="1") g=g"MALE"
      if (length(anc)) g=(g?g"," anc:anc)
      print s "\\t" g
    }' "${metadata_csv}" > groups.txt

    echo "✓ Groups file created" >> "${log_file}"

    {
      echo "=== Adding allele frequencies to VCF (dropping all FORMAT/GT columns) ==="
    } >> "${log_file}"

    # Recalculate AF/AC/AN (and other tags) across all samples AND per-group (via -S groups.txt),
    # then drop genotypes (-G) and strip any FORMAT header remnants (-x FORMAT).  :contentReference[oaicite:3]{index=3}
    bcftools +fill-tags "\${INPUT_VCF}" -Ou -- -S groups.txt \
      | bcftools view -G -Ou \
      | bcftools annotate -x FORMAT \
      -Oz -o "\$OUTPUT_VCF"

    # Index the INFO-only VCF
    tabix -p vcf "\$OUTPUT_VCF"

    {
      echo "✓ AF annotation complete. Output: \$OUTPUT_VCF"
    } >> "${log_file}"
    """
}


workflow {
    // Channel: all .vcf.gz in params.input 
    vcf_ch = Channel.fromPath("${params.input}/*.vcf.gz", checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            def log_file = "${params.output_stats}/${vcf.simpleName}.log"
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

    // Add AF annotations
    af_annotated = ADD_AF(sample_qc, metadata_ch)
}
