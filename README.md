# BCFtools-Pipeline (Nextflow DSL2)

This repository contains a **Nextflow DSL2 pipeline** for variant calling quality control (QC) and downstream allele frequency annotation stratified by sex and ancestry.  

It integrates several QC steps (genotype-, variant-, and sample-level filtering) and recalculates allele frequencies (AFs) based on sample metadata using BCFTools. 



## Workflow Overview

The pipeline performs the following steps:

1. **INDEX_VCF**  
   - Ensures each input VCF has an index (`.tbi`).  
   - If missing, it creates one with `bcftools index`.

2. **SPLIT_MULTIALLELIC**  
   - Splits multiallelic variants into biallelic records using `bcftools norm`.

3. **GENOTYPE_QC**  
   - Masks low-quality genotypes (`./.`) according to thresholds:
     - Genotype Quality (GQ)
     - Depth of Coverage (DP)
     - Allelic Depth ratio (AD)  
   - Automatically detects if FORMAT tags exist before applying filters.

4. **VARIANT_QC**  
   - Filters out low-quality variants based on INFO and QUAL fields.  
   - Supported filters: QUAL, QD, DP, MQ, FS, ReadPosRankSum.

5. **SAMPLE_QC**  
   - Removes samples failing quality thresholds:
     - Mean coverage
     - Call rate
     - Het/Hom ratio
     - Number of singletons
     - Contamination (using [sceVCF](https://github.com/HTGenomeAnalysisUnit/SCE-VCF))  

   **Thresholds differ depending on `seq_type` (WES vs WGS).**

6. **ADD_AF**  
   - Uses a metadata CSV file to define groups (sex + ancestry).  
   - Annotates the final VCF with total and group-specific allele frequencies.

## Technical requirements 

- [Nextflow](https://nextflow.io/docs/latest/install.html) ≥ 24.04

- [Java](https://www.oracle.com/es/java/technologies/downloads/) ≥ 17

- [bcftools](https://samtools.github.io/bcftools/howtos/install.html) ≥ 1.10

- [tabix](https://www.htslib.org/doc/tabix.html) (htslib) ≥ 1.13

- [sceVCF binaries](https://github.com/HTGenomeAnalysisUnit/SCE-VCF/releases/tag/v0.1.3) for contamination checks

## Input Requirements

- **VCF files** (`.vcf.gz` or `.vcf.bgz`) in the input directory.  

    ⚠️ **Important: For allele frequency (AF) recalculations to be valid, the input VCFs should contain only unrelated individuals.**

- **Index files** (`.tbi`) [OPTIONAL], will be created if missing.  
- **Metadata CSV** with columns SAMPLE_ID, SEX and ANCESTRY.

    Example: 
    ```
    SAMPLE,SEX,ANCESTRY
    sample1,1,EUR
    sample2,2,AFR
    sample3,1,EAS
    ```
    **Some important considerations about the submitted metadata:** 
    
    - Ensure this is the order of your metadata columns.
    - Sex should be coded as `1` = MALE, `2` = FEMALE.  
    - Ancestry codes are free-text (e.g., EUR, AFR, EAS). But, take into consideration that the naming on the CSV will be used to annotate the ancestry AF fields in the VCF. 


- **sceVCF binary path** 

## Configuration 

Edit *nextflow.config* to adjust parameters. 

```groovy
params {
    // Paths
    input        = "/path/to/folder/with/vcfs/"
    sceVCF_path  = "/path/to/sceVCF_binaries"
    metadata_csv = "/path/to/metadata.csv"
    seq_type     = "WES" // or "WGS"
    threads      = 4

     // Pipeline behavior
    seq_type     = 'WGS' // or 'WES'
    threads      = 4

     // -----------------------
    // QC Thresholds
    // -----------------------
    qc {

        // Genotype QC defaults 
        genotype {
        gq_threshold         = 20
        dp_threshold         = 10
        ad_ratio_threshold   = 0.2
        }

        // Variant QC defaults
        variant {
        qual_threshold                 = 30
        qd_threshold                   = 2.0
        dp_threshold                   = 10
        mq_threshold                   = 40
        fs_threshold                   = 60
        read_pos_rank_sum_threshold    = -8.0
        }

        // Sample QC (choose set via `seq_type`)
        sample {
        // Whole-Exome defaults
        wes {
            coverage_threshold       = 10
            het_hom_threshold        = 10
            call_rate_threshold      = 0.95
            singletons_threshold     = 5000
            contamination_threshold  = 0.00015
        }
        // Whole-Genome defaults
        wgs {
            coverage_threshold       = 15
            het_hom_threshold        = 3.3
            call_rate_threshold      = 0.95
            singletons_threshold     = 100000
            contamination_threshold  = 0.05
        }
        }
    }
}

```
**Notes**

- `seq_type` selects which Sample QC thresholds are applied (wes vs wgs).

- metadata_csv is required for the ADD_AF step.

## Running the Pipeline

### Basic run: 

```bash
nextflow run main.nf
```

### With custom parameters

```bash
nextflow run main.nf \
  --input "/data/vcfs/" \
  --metadata_csv "/data/metadata.csv" \
  --sceVCF_path "/tools/sceVCF" \
  --seq_type "WGS" \
  --threads 8
```

**NOTE:** This will overwrite the configuration parameters set in nextflow.config

## Outputs 

After successful execution, you'll find inside the */work* folder:

**Intermediate outputs**

- Indexed VCFs (.tbi)

- Split multiallelic VCFs

- Masked genotypes VCFs

- Variant-filtered and masked genotypes VCFs

- Variant and Sample-filtered with masked genotypes VCFs

**Final output**

* with_AF.vcf.gz (TODO: make the name the same as the original VCF with -AF_recalc) → fully filtered VCF with allele frequency annotations and witout sample level information. 

- with_AF.vcf.gz.tbi → index for the final VCF.


---
💡 How to find the files ? 

Examples of nextflow output: 
