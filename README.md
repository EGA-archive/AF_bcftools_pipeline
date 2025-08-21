# BCFtools-Pipeline (Nextflow DSL2)

This repository contains a **Nextflow DSL2 pipeline** for **SNP** variant calling quality control (QC) and downstream allele frequency annotation stratified by sex and ancestry.  

It integrates several QC steps (genotype-, variant-, and sample-level filtering) and recalculates allele frequencies (AFs) based on sample metadata using BCFTools. 



## Workflow Overview

The pipeline performs the following steps:

1. **INDEX_VCF**  
   - Ensures each input VCF has an index (`.tbi`).  
   - If missing, it creates one with `tabix`.

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

- **VCF files** (`.vcf.gz`) in the input directory.  

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
    
    - Ensure the order of your metadata columns is the one shown above.
    - Sex should be coded as `1` = MALE, `2` = FEMALE.  
    - Ancestry codes are free-text (e.g., EUR, AFR, EAS). But, take into consideration that the naming on the CSV will be used to annotate the ancestry AF fields in the VCF. 


- **sceVCF binary path** 

## Configuration 

Edit *nextflow.config* to adjust parameters. 

```groovy
params {
    // Paths
    input        = "/path/to/folder/with/vcfs/"
    output_stats = "/home/mireia/GitHub/bcftools-pipeline/bcftools-pipeline"
    sceVCF_path  = "/path/to/sceVCF_binaries"
    metadata_csv = "/path/to/metadata.csv"
    seq_type     = "WGS" // or "WES"
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

- metadata_csv is required.

- The fields INFO/QD, INFO/DP, INFO/MQ, INFO/FS, INFO/ReadPosRankSum,  FORMAT/GQ, FORMAT/GT and FORMAT/AD must be correctly described in the header for the filtering to work correctly.

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

### HTML report and timeline

```
nextflow run main.nf -with-report report.html -with-timeline timeline.html
```

This command will automatically generate two reports, report.html and timeline.html. 

**report.html** will contain: 

- Process run counts and statuses (cached/failed/succeeded)

- CPU, memory, time usage per process

- I/O stats and container info (if used)

- Cache hit ratio and overall runtime

**timeline.html** will contain:

- Each task’s start/finish times

- Parallelism and queueing

- Retries/restarts and durations


## Outputs 

After a successful execution, you'll find inside the */work* folder:

**Intermediate outputs**

- Indexed VCFs (.tbi)

- Split multiallelic VCFs

- Masked genotypes VCFs

- Variant-filtered and masked genotypes VCFs

- Variant and Sample-filtered with masked genotypes VCFs

**Final output**

* *input_vcf_baseName*.vcf.gz (TODO: make the name the same as the original VCF with -AF_recalc) → fully filtered VCF with allele frequency annotations and witout sample level information. 

- *input_vcf_baseName*.vcf.gz.tbi → index for the final VCF.

- In `$output_stats` there will be a log file with information about each step of the pipeline run per file. 


---
💡 **How to find the files inside the /work folder?**

Example a of nextflow output: 

```
[34/34b365] process > INDEX_VCF (test.vcf.bgz)                                           [100%] 1 of 1 ✔
[bc/56f27a] process > SPLIT_MULTIALLELIC (test.vcf.bgz)                                  [100%] 1 of 1 ✔
[46/bf10e7] process > GENOTYPE_QC (test_split-multiallelic.vcf.bgz)                      [100%] 1 of 1 ✔
[88/8b101c] process > VARIANT_QC (test_split-multiallelic-masked.vcf.gz)                 [100%] 1 of 1 ✔
[33/a6c930] process > SAMPLE_QC (test_split-multiallelic-masked-filtered.vcf.gz)         [100%] 1 of 1 ✔
[c7/065171] process > ADD_AF (test_split-multiallelic-masked-filtered.sample_qc.vcf.bgz) [100%] 1 of 1 ✔
```
Go to /work/bc/56f27a[-->] and you'll find the VCF with the multiallelic variants splitted.

Go to /work/c7/065171[-->] for the final files of the pipeline. 

## How to generate metadata_csv

If you don't have the necessary metadata to create `metadata.csv` you can infer it using tools such as [GRAF tools](https://github.com/ncbi/graf) and [Hail](https://hail.is/docs/0.2/index.html). 

### Graf tools
GRAF provides functions to infer sex (graf sex, PLINK input), detect related samples (graf rel, PLINK input), and assign ancestry ([graf anc](https://github.com/jimmy-penn/grafanc/tree/master), VCF input). Note that, for ancestry inference, only super-population calls (e.g. European, Asian, African...) are recommended, since finer-level predictions are not sufficiently accurate. 

### Hail
Hail offers similar capabilities with `impute_sex` for sex inference and `king` and `pc_relate` for relatedness filtering.

## References

Lu W, Gauthier LD, Poterba T, Giacopuzzi E, Goodrich JK, Stevens CR, King D, Daly MJ, Neale BM, Karczewski KJ. CHARR efficiently estimates contamination from DNA sequencing data. bioRxiv [Preprint]. 2023 Jun 28:2023.06.28.545801. doi: 10.1101/2023.06.28.545801. Update in: Am J Hum Genet. 2023 Dec 7;110(12):2068-2076. doi: 10.1016/j.ajhg.2023.10.011. PMID: 37425834; PMCID: PMC10327099.