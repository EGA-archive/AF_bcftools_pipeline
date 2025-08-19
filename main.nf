nextflow.enable.dsl=2

process INDEX_VCF {
    tag "$vcf"

    input:
    tuple path(vcf), val(has_index), path(vcf_tbi, stageAs: "existing.tbi")
    
    output:
    tuple path(vcf), path("${vcf}.tbi")
    
    debug true
    script:

    if (has_index) {
        """
        echo ""
        echo "✓ Index already exists for ${vcf}"
        # Copy existing index to expected output name
        cp existing.tbi "${vcf}.tbi"
        """
    } else {
        """
        set -euo pipefail 
        echo "✗ No index found for ${vcf} — creating..."
        bcftools index --threads ${task.cpus} -t -f "${vcf}"
        echo ""
        echo "✓ Index created for ${vcf}" 
        """
    }
}

process SPLIT_MULTIALLELIC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${vcf.simpleName}_split-multiallelic.vcf.bgz"),
          path("${vcf.simpleName}_split-multiallelic.vcf.bgz.tbi")

    debug true
    script:
    """
    bcftools norm -m -any ${vcf} -Oz -o ${vcf.simpleName}_split-multiallelic.vcf.bgz
    tabix -p vcf ${vcf.simpleName}_split-multiallelic.vcf.bgz
    """
}

process GENOTYPE_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${vcf.simpleName}-masked.vcf.gz"), path("${vcf.simpleName}-masked.vcf.gz.tbi")

    debug true
    script:
    """
    set -euo pipefail

    VCF_IN="${vcf}"
    VCF_OUT="${vcf.simpleName}-masked.vcf.gz"

    # Join operator for per-genotype rules:
    # Use "||" to mask a genotype if it fails ANY rule (typical); use "&&" to mask only if ALL rules fail.
    OP="||"

    # --------- FORMAT-based masking rules only (per-genotype) ----------
    # NOTE: Only FORMAT/* expressions make sense for +setGT masking.
    # AD rule includes zero-division guard.
    declare -A gt_conditions=(
      [GQ]='FMT/GQ < ${params.qc.genotype.gq_threshold}'
      [DP]='FMT/DP < ${params.qc.genotype.dp_threshold}'
      [AD]='(FMT/AD[0]+FMT/AD[1])>0 && (FMT/AD[1]/(FMT/AD[0]+FMT/AD[1])) < ${params.qc.genotype.ad_ratio_threshold}'
    )

    echo "=== Checking FORMAT fields in \$VCF_IN for genotype masking ==="
    gt_expr_parts=()

    # Add a condition only if its FORMAT tag exists in the header
    for tag in "\${!gt_conditions[@]}"; do
      if bcftools view -h "\$VCF_IN" | grep -q "^##FORMAT=<ID=\${tag},"; then
        echo "✓ \${tag} (FORMAT) found — adding: \${gt_conditions[\$tag]}"
        gt_expr_parts+=("\${gt_conditions[\$tag]}")
      else
        echo "x \${tag} (FORMAT) not found — skipping"
      fi
    done

    # Build final per-genotype expression
    if (( \${#gt_expr_parts[@]} == 0 )); then
      echo "x No FORMAT-based rules available; no masking performed."
      # still emit an output identical to input, with index
      cp -a "\$VCF_IN" "\$VCF_OUT"
      cp -a "\${VCF_IN}.tbi" "\${VCF_OUT}.tbi" 2>/dev/null || tabix -p vcf "\$VCF_OUT"
      exit 0
    fi

    gt_expr="\${gt_expr_parts[0]}"
    for cond in "\${gt_expr_parts[@]:1}"; do
      gt_expr+=" \$OP \$cond"
    done

    echo "=== Final per-genotype mask expression (joined with '\$OP') ==="
    echo "\$gt_expr"
    echo "=== Running bcftools +setGT (mask to ./.) ==="

    # If you have a sex-aware ploidy file, add:  --ploidy-file ploidy.txt
    bcftools +setGT "\$VCF_IN" -Oz -o "\$VCF_OUT" -- \
      -t q -n . -e "\$gt_expr"

    tabix -p vcf "\$VCF_OUT"
    echo "✓ Genotype masking complete. Output: \$VCF_OUT"
    """
}

process VARIANT_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("${vcf.simpleName}-filtered.vcf.gz"), path("${vcf.simpleName}-filtered.vcf.gz.tbi")

    debug true
    script:
    """
    set -euo pipefail

    VCF_IN="${vcf}"
    VCF_OUT="${vcf.simpleName}-filtered.vcf.gz"
    OP="||"

    declare -A site_conditions=(
      [QD]='INFO/QD < ${params.qc.variant.qd_threshold}'
      [DP]='INFO/DP < ${params.qc.variant.dp_threshold}'
      [MQ]='INFO/MQ < ${params.qc.variant.mq_threshold}'
      [FS]='INFO/FS > ${params.qc.variant.fs_threshold}'
      [ReadPosRankSum]='INFO/ReadPosRankSum < ${params.qc.variant.read_pos_rank_sum_threshold}'
    )

    has_info() {
      bcftools view -h "\$1" | grep -q "^##INFO=<ID=\$2,"
    }

    echo "=== Checking site-level fields in \$VCF_IN ==="

    expr_parts=("QUAL < ${params.qc.variant.qual_threshold}")
    echo "✓ QUAL — adding: QUAL < ${params.qc.variant.qual_threshold}"

    for tag in "\${!site_conditions[@]}"; do
      if has_info "\$VCF_IN" "\$tag"; then
        echo "✓ \$tag (INFO) found — adding: \${site_conditions[\$tag]}"
        expr_parts+=("\${site_conditions[\$tag]}")
      else
        echo "x \$tag (INFO) not found — skipping"
      fi
    done

    expr_str="\${expr_parts[0]}"
    for cond in "\${expr_parts[@]:1}"; do
      expr_str+=" \$OP \$cond"
    done

    echo "=== Final filter expression ==="
    echo "\$expr_str"

    bcftools filter -e "\$expr_str" "\$VCF_IN" -Oz -o "\$VCF_OUT"
    tabix -p vcf "\$VCF_OUT"

    echo "✓ Filtering complete. Output: \$VCF_OUT"
    echo "Variants after QC: \$(bcftools index -n "\$VCF_OUT")"
    """
}

process SAMPLE_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi)
    path(sceVCF_path)
    val seq_type

    output:
    tuple path("${vcf.simpleName}.sample_qc.vcf.bgz"),
          path("${vcf.simpleName}.sample_qc.vcf.bgz.tbi")

    debug true
    script:
    """
    set -euo pipefail
    INPUT_VCF="${vcf}"
    OUTPUT_VCF="${vcf.simpleName}.sample_qc.vcf.bgz"
    TMP="qc_tmp"; mkdir -p "\$TMP"

    echo "=== Starting SAMPLE QC on: \$INPUT_VCF ==="

    if [[ "${seq_type}" == "WGS" ]]; then
      COV_THRESHOLD=${params.qc.sample.wgs.coverage_threshold}
      HET_HOM_THRESHOLD=${params.qc.sample.wgs.het_hom_threshold}
      CALL_RATE_THRESHOLD=${params.qc.sample.wgs.call_rate_threshold}
      SINGLETONS_THRESHOLD=${params.qc.sample.wgs.singletons_threshold}
      CONTAM_THRESHOLD=${params.qc.sample.wgs.contamination_threshold}
      echo "✓ Seq type: WGS thresholds applied"
    else
      COV_THRESHOLD=${params.qc.sample.wes.coverage_threshold}
      HET_HOM_THRESHOLD=${params.qc.sample.wes.het_hom_threshold}
      CALL_RATE_THRESHOLD=${params.qc.sample.wes.call_rate_threshold}
      SINGLETONS_THRESHOLD=${params.qc.sample.wes.singletons_threshold}
      CONTAM_THRESHOLD=${params.qc.sample.wes.contamination_threshold}
      echo "✓ Seq type: WES thresholds applied"
    fi

    SITE_SUBSET_CMD=(bcftools view --threads ${task.cpus} -f PASS -v snps "\$INPUT_VCF")

    # 1) Mean coverage
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=DP,"; then
      echo "✓ DP (FORMAT) found — calculating mean coverage"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%DP\\n]' \
      | awk '{sum[\$1]+=\$2; n[\$1]++} END{for(s in sum){if(n[s]>0) printf "%s\\t%.6f\\n",s,sum[s]/n[s]}}' \
      > "\$TMP/mean_dp.txt"
      awk -v thr="\$COV_THRESHOLD" '\$2 < thr {print \$1}' "\$TMP/mean_dp.txt" > "\$TMP/low_cov_samples.txt"
    else
      echo "x DP (FORMAT) not found — skipping"
      : > "\$TMP/mean_dp.txt"
      : > "\$TMP/low_cov_samples.txt"
    fi

    # 2) Call rate
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=GT,"; then
      echo "✓ GT (FORMAT) found — calculating call rate"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%GT\\n]' \
      | awk '{tot[\$1]++; if(\$2=="./."||\$2==".|.") miss[\$1]++} END{for(s in tot){cr=1-((miss[s]+0)/tot[s]); printf "%s\\t%.6f\\n",s,cr}}' \
      > "\$TMP/call_rate.txt"
      awk -v thr="\$CALL_RATE_THRESHOLD" '\$2 < thr {print \$1}' "\$TMP/call_rate.txt" > "\$TMP/low_call_rate.txt"
    else
      echo "x GT (FORMAT) not found — skipping call rate"
      : > "\$TMP/call_rate.txt"
      : > "\$TMP/low_call_rate.txt"
    fi

    # 3) Het/Hom ratio
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=GT,"; then
      echo "✓ GT (FORMAT) found — calculating Het/Hom ratio"
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
          for(s in het) {
            if(hom[s]>0) printf "%s\\t%.6f\\n",s,het[s]/hom[s];
            else if(het[s]>0) printf "%s\\tinf\\n",s;
            else printf "%s\\t0\\n",s;
          }
          for(s in hom) if(!(s in het)) printf "%s\\t0\\n",s;
        }' \
      > "\$TMP/het_hom.txt"
      awk -v thr="\$HET_HOM_THRESHOLD" '\$2=="inf" || \$2+0 > thr {print \$1}' "\$TMP/het_hom.txt" > "\$TMP/bad_het_hom.txt"
    else
      echo "x GT (FORMAT) not found — skipping Het/Hom ratio"
      : > "\$TMP/het_hom.txt"
      : > "\$TMP/bad_het_hom.txt"
    fi

    # 4) Singletons
    echo "✓ Counting singletons"
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

    # 5) Contamination
    if [[ -z "${sceVCF_path}" || "${sceVCF_path}" == "" ]]; then
        echo "x Contamination check not running (sceVCF path not provided)"
        : > "\$TMP/high_contam.txt"
    else
        # Check if AD format field is present
        if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=AD,"; then
            echo "✓ AD (FORMAT) found — checking sceVCF availability"
            if [[ -f "${sceVCF_path}" && -x "${sceVCF_path}" ]]; then
                echo "✓ sceVCF found at ${sceVCF_path} — running contamination check"
                export PATH="${sceVCF_path}:\$PATH" # add sceVCF executable to path
                ./sceVCF -o "\$TMP/charr_full.tsv" "\$INPUT_VCF"
                awk -v thr="\$CONTAM_THRESHOLD" '\$2 > thr {print \$1}' "\$TMP/charr_full.tsv" > "\$TMP/high_contam.txt"
            else
                echo "x sceVCF not found at ${sceVCF_path} or not executable — skipping contamination check"
                : > "\$TMP/high_contam.txt"
            fi
        else
            echo "x AD (FORMAT) not found — skipping contamination check (required for sceVCF)"
            : > "\$TMP/high_contam.txt"
        fi
    fi

    echo "=== Merging flagged samples ==="
    cat "\$TMP"/low_cov_samples.txt "\$TMP"/low_call_rate.txt "\$TMP"/bad_het_hom.txt "\$TMP"/high_singletons.txt "\$TMP/high_contam.txt" \
      2>/dev/null | sort -u > "\$TMP/samples_to_remove.txt"
    : > "\$TMP/samples_to_remove.txt" || true

    if [[ -s "\$TMP/samples_to_remove.txt" ]]; then
      echo "✗ Removing \$(wc -l < "\$TMP/samples_to_remove.txt") samples failing QC"
      bcftools view --threads ${task.cpus} -S ^"\$TMP/samples_to_remove.txt" "\$INPUT_VCF" -Oz -o "\$OUTPUT_VCF"
      bcftools index -ft "\$OUTPUT_VCF"
    else
    echo "✓ No samples flagged — renaming input VCF as output"
    mv "\$INPUT_VCF" "\$OUTPUT_VCF"
    if [[ -f "\${INPUT_VCF}.tbi" ]]; then
        mv "\${INPUT_VCF}.tbi" "\$OUTPUT_VCF.tbi"
    else
        tabix -p vcf "\$OUTPUT_VCF"
    fi
    fi

    echo "✓ SAMPLE QC complete. Output: \$OUTPUT_VCF"
    """
}

process ADD_AF {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi)
    path metadata_csv

    output:
    tuple path("with_AF.vcf.gz"), path("with_AF.vcf.gz.tbi")

    debug true 
    script:
    """
    set -euo pipefail

    echo "=== AF recalculating: Creating groups.txt from metadata ==="
    
    # Create groups.txt with 2 columns: SAMPLE \t comma-separated GROUPS
    awk -F, 'NR>1{
      s=\$1
      sex=tolower(\$2)
      anc=\$3
      g=""
      if (sex=="2") g=g"FEMALE"
      else if (sex=="1") g=g"MALE"
      if (length(anc)) g=(g?g"," anc:anc)
      print s "\\t" g
    }' ${metadata_csv} > groups.txt

    echo "✓ Groups file created"

    echo "=== Adding allele frequencies to VCF ==="
    
    # Write a new VCF with total AF + per-group AFs in INFO
    bcftools +fill-tags "${vcf}" -Oz -o with_AF.vcf.gz -- \\
      -S groups.txt

    # Index the output VCF
    tabix -p vcf with_AF.vcf.gz

    echo "✓ AF annotation complete. Output: with_AF.vcf.gz"
    """
}

workflow {
    // Channel: all .vcf.gz or .vcf.bgz in params.input 
    vcf_ch = Channel.fromPath("${params.input}/*.{vcf.gz,vcf.bgz}", checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            if (tbi.exists()) {
                // Both VCF and index exist
                tuple(vcf, true, tbi)
            } else {
                // Only VCF exists, use dummy file for index
                tuple(vcf, false, file("NO_FILE"))
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