nextflow.enable.dsl=2

// Ensure the output dir exists (Groovy)
file(params.output_stats).mkdirs()

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

    bcftools norm -m -any "${vcf}" -Oz -o "${vcf.simpleName}_split-multiallelic.vcf.gz"
    tabix -p vcf "${vcf.simpleName}_split-multiallelic.vcf.gz"

    {
      echo "✓ Split + index done"
      echo "Output: ${vcf.simpleName}_split-multiallelic.vcf.gz"
      echo "Variants now: \$(bcftools index -n "${vcf.simpleName}_split-multiallelic.vcf.gz")"
    } >> "${log_file}"
    """
}

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

    declare -A gt_conditions=(
      [GQ]='FMT/GQ < ${params.qc.genotype.gq_threshold}'
      [DP]='FMT/DP < ${params.qc.genotype.dp_threshold}'
      [AD]='(FMT/AD[0]+FMT/AD[1])>0 && (FMT/AD[1]/(FMT/AD[0]+FMT/AD[1])) < ${params.qc.genotype.ad_ratio_threshold}'
    )

    {
      echo ""
      echo "=== GENOTYPE_QC on: \$VCF_IN ==="
    } >> "${log_file}"

    gt_expr_parts=()
    for tag in "\${!gt_conditions[@]}"; do
      if bcftools view -h "\$VCF_IN" | grep -q "^##FORMAT=<ID=\${tag},"; then
        echo "✓ \${tag} (FORMAT) found — adding rule" >> "${log_file}"
        gt_expr_parts+=("\${gt_conditions[\$tag]}")
      else
        echo "x \${tag} (FORMAT) not found — skipping" >> "${log_file}"
      fi
    done

    if (( \${#gt_expr_parts[@]} == 0 )); then
      {
        echo "x No FORMAT-based rules available; no masking performed."
      } >> "${log_file}"
      cp -a "\$VCF_IN" "\$VCF_OUT"
      cp -a "\${VCF_IN}.tbi" "\${VCF_OUT}.tbi" 2>/dev/null || tabix -p vcf "\$VCF_OUT"
      exit 0
    fi

    gt_expr="\${gt_expr_parts[0]}"
    for cond in "\${gt_expr_parts[@]:1}"; do gt_expr+=" \$OP \$cond"; done

    {
      echo "Final per-genotype mask expression:"
      echo "\$gt_expr"
      echo "Running bcftools +setGT (mask to ./.)"
    } >> "${log_file}"

    bcftools +setGT "\$VCF_IN" -Oz -o "\$VCF_OUT" -- -t q -n . -e "\$gt_expr"
    tabix -p vcf "\$VCF_OUT"

    {
      echo "✓ Genotype masking complete."
      echo "Output: \$VCF_OUT"
    } >> "${log_file}"
    """
}

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
    OP="||"

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

    expr_parts=("QUAL < ${params.qc.variant.qual_threshold}")
    echo "✓ QUAL — adding: QUAL < ${params.qc.variant.qual_threshold}" >> "${log_file}"

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

    bcftools filter -e "\$expr_str" "\$VCF_IN" -Oz -o "\$VCF_OUT"
    tabix -p vcf "\$VCF_OUT"

    {
      echo ""
      echo "✓ Filtering complete. Output: \$VCF_OUT"
      echo "Variants after QC: \$(bcftools index -n "\$VCF_OUT")"
    } >> "${log_file}"
    """
}

process SAMPLE_QC {
    tag "$vcf"

    input:
    tuple path(vcf), path(tbi), val(log_file)
    path(sceVCF_path)
    val seq_type

    output:
    tuple path("${vcf.simpleName}-sampleQC.vcf.gz"),
          path("${vcf.simpleName}-sampleQC.vcf.gz.tbi"),
          val(log_file)

    debug true
    script:
    """
    set -euo pipefail
    INPUT_VCF="${vcf}"
    OUTPUT_VCF="${vcf.simpleName}-sampleQC.vcf.gz"
    TMP="qc_tmp"; mkdir -p "\$TMP"

    {
      echo ""
      echo "=== SAMPLE_QC on: \$INPUT_VCF ==="
    } >> "${log_file}"

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

    SITE_SUBSET_CMD=(bcftools view --threads ${task.cpus} -f PASS -v snps "\$INPUT_VCF")

    # 1) Mean coverage
    if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=DP,"; then
      echo "✓ DP (FORMAT) found — calculating mean coverage" >> "${log_file}"
      "\${SITE_SUBSET_CMD[@]}" \
      | bcftools query -f '[%SAMPLE\\t%DP\\n]' \
      | awk '{sum[\$1]+=\$2; n[\$1]++} END{for(s in sum){if(n[s]>0) printf "%s\\t%.6f\\n",s,sum[s]/n[s]}}' \
      > "\$TMP/mean_dp.txt"
      awk -v thr="\$COV_THRESHOLD" '\$2 < thr {print \$1}' "\$TMP/mean_dp.txt" > "\$TMP/low_cov_samples.txt"
    else
      echo "x DP (FORMAT) not found — skipping" >> "${log_file}"
      : > "\$TMP/mean_dp.txt"; : > "\$TMP/low_cov_samples.txt"
    fi

    # 2) Call rate
    GT_FOUND=\$(bcftools view -h "\$INPUT_VCF" | grep -c 'ID=GT')
    if [ "\$GT_FOUND" -gt 0 ]; then
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

    # 3) Het/Hom ratio
    if [ "\$GT_FOUND" -gt 0 ]; then
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
      echo "x GT (FORMAT) not found — skipping Het/Hom ratio" >> "${log_file}"
      : > "\$TMP/het_hom.txt"; : > "\$TMP/bad_het_hom.txt"
    fi

    # 4) Singletons
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

    # 5) Contamination
    if [[ -z "${sceVCF_path}" || "${sceVCF_path}" == "" ]]; then
        echo "x Contamination check not running (sceVCF path not provided)" >> "${log_file}"
        : > "\$TMP/high_contam.txt"
    else
        if bcftools view -h "\$INPUT_VCF" | grep -q "^##FORMAT=<ID=AD,"; then
            echo "✓ AD (FORMAT) found — checking sceVCF availability" >> "${log_file}"
            if [[ -f "${sceVCF_path}" && -x "${sceVCF_path}" ]]; then
                echo "✓ sceVCF found at ${sceVCF_path} — running contamination check" >> "${log_file}"
                export PATH="${sceVCF_path}:\$PATH"
                ./sceVCF -o "\$TMP/charr_full.tsv" "\$INPUT_VCF"
                awk -v thr="\$CONTAM_THRESHOLD" '\$2 > thr {print \$1}' "\$TMP/charr_full.tsv" > "\$TMP/high_contam.txt"
            else
                echo "x sceVCF not found at ${sceVCF_path} or not executable — skipping" >> "${log_file}"
                : > "\$TMP/high_contam.txt"
            fi
        else
            echo "x AD (FORMAT) not found — skipping contamination check (required for sceVCF)" >> "${log_file}"
            : > "\$TMP/high_contam.txt"
        fi
    fi

    cat "\$TMP"/low_cov_samples.txt "\$TMP"/low_call_rate.txt "\$TMP"/bad_het_hom.txt "\$TMP"/high_singletons.txt "\$TMP/high_contam.txt" \
    2>/dev/null | sort -u > "\$TMP/samples_to_remove.txt"

    if [[ -s "\$TMP/samples_to_remove.txt" ]]; then
      echo "✗ Removing \$(wc -l < "\$TMP/samples_to_remove.txt") samples failing QC" >> "${log_file}"
      bcftools view --threads ${task.cpus} -S ^"\$TMP/samples_to_remove.txt" "\$INPUT_VCF" -Oz -o "\$OUTPUT_VCF"
      bcftools index -ft "\$OUTPUT_VCF"
    else
      {
        echo "✓ No samples flagged — renaming input VCF as output"
      } >> "${log_file}"
      mv "\$INPUT_VCF" "\$OUTPUT_VCF"
      if [[ -f "\${INPUT_VCF}.tbi" ]]; then
          mv "\${INPUT_VCF}.tbi" "\$OUTPUT_VCF.tbi"
      else
          tabix -p vcf "\$OUTPUT_VCF"
      fi
    fi

    {
      echo "✓ SAMPLE_QC complete. Output: \$OUTPUT_VCF"
    } >> "${log_file}"
    """
}

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

    echo "✓ Groups file created" >> "${log_file}"

    {
      echo "=== Adding allele frequencies to VCF (dropping FORMAT cols) ==="
    } >> "${log_file}"

    bcftools +fill-tags "${vcf}" -- -S groups.txt \
      | bcftools view -G \
      | bcftools annotate --remove 'FORMAT' \
      -Oz -o "\$OUTPUT_VCF"

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
