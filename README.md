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

6. **FIX_PLOIDY**

    - Metadata Conversion: Generates a gender.txt mapping (Sample ID → Sex) from the input metadata.

    - Ploidy Correction: Uses bcftools +fixploidy to ensure sex chromosomes are handled correctly. This converts diploid "heterozygous" calls in male non-PAR regions into proper haploid calls, which is essential for accurate Allele Frequency (AF) calculation.

    - Genome-Aware: Automatically selects PAR (Pseudoautosomal Region) coordinates for GRCh37 or GRCh38 based on the reference_genome parameter.

    💡 Customization: If your data uses non-standard chromosome names or specific PAR boundaries, you can modify the coordinate files in the /assets directory.



7. **ADD_AF**  
   - Uses a metadata CSV file to define groups (sex + ancestry).  
   - Annotates the final VCF with total and group-specific allele frequencies.

## Installation

### OPTION 1: Clone the repository

Technical requirements:

- [Nextflow](https://nextflow.io/docs/latest/install.html) ≥ 26.10

- [Java](https://www.oracle.com/es/java/technologies/downloads/) ≥ 17

- [bcftools](https://samtools.github.io/bcftools/howtos/install.html) ≥ 1.24

- [tabix](https://www.htslib.org/doc/tabix.html) (htslib) ≥ 1.24

- [sceVCF binaries](https://github.com/HTGenomeAnalysisUnit/SCE-VCF/releases/tag/v0.1.3) for contamination checks

### OPTION 2: Docker

Clone all the contents of the folder and build the container: 

```
docker build --tag af_bcftools_pipeline:latest .
```

To check the installation was correct: 

```
docker image ls
```

In the image list you should see: `af_bcftools_pipeline:latest`

## Input Requirements

- **VCF files** (`.vcf.gz`) in the input directory.  

    ⚠️ **Important: For allele frequency (AF) recalculations to be valid, the input VCFs should contain only unrelated individuals.**

- **Index files** (`.tbi`) [OPTIONAL], will be created if missing.  
- **Metadata CSV** with columns SAMPLE_ID, SEX and ANCESTRY.

    Example: 
    ```
    SAMPLE,SEX,ANCESTRY
    sample1,M,EUR
    sample2,F,AFR
    sample3,M,EAS
    ```
    **Some important considerations about the submitted metadata:** 
    
    - Ensure the order of your metadata columns is the one shown above.
    - Sex should be coded as `M` = MALE, `F` = FEMALE.  
    - Ancestry codes are free-text (e.g., EUR, AFR, EAS). But, take into consideration that the naming on the CSV will be used to annotate the ancestry AF fields in the VCF. 


- **sceVCF binary path** 

## Running the pipeline

### Configure the pipeline

Edit `nextflow.config` to adjust parameters. 

```groovy
// Nextflow configuration file
nextflow.enable.dsl = 2

// Default parameters
params {
    // Input/Output paths
    input = "/data/input" // move input VCFs to this folder.  
    output = "/data/output" // output directory
    sceVCF_path = "/usr/local/bin/sceVCF" // File downloaded by clonning. 
    metadata_csv =  "/data/input/mock_metadata_pop13_1.csv" // header: sample/sex/ancestry. 
    
    // Pipeline parameters
    seq_type = 'WGS'  // or 'WES'
    threads = 4
    reference_genome = 'GRCh37' // or 'GRCh38'
    
    // Quality control thresholds
    qc {
        // Genotype QC thresholds
        genotype {
            gq_threshold = 20
            dp_threshold = 10
            ab_ratio_min_threshold = 0.2 // if no lower theshold is needed, set ab_ratio_max_threshold = 1
            ab_ratio_max_threshold = 0.8 // if no upper theshold is needed, set ab_ratio_max_threshold = 1
        }
        // Variant QC thresholds  
        variant {
            qual_threshold = 30
            qd_threshold = 2.0
            dp_threshold = 10
            mq_threshold = 40
            fs_threshold = 60
            read_pos_rank_sum_threshold = -8.0
        }
        
        // Sample QC thresholds (WES)
        sample {
            wes {
                coverage_threshold = 10
                het_hom_threshold = 10
                call_rate_threshold = 0.95
                singletons_threshold = 5000
                contamination_threshold = 0.00015
            }
            // Sample QC thresholds (WGS)
            wgs {
                coverage_threshold = 15
                het_hom_threshold = 3.3
                call_rate_threshold = 0.95
                singletons_threshold = 100000
                contamination_threshold = 0.05
            }
        }
    }
}

// Process-specific configurations
process {
    // Default process settings
    cpus = params.threads
    memory = '4.GB'
    
    withName: INDEX_VCF {
        cpus = params.threads
        memory = '2.GB'
    }
    
    withName: SPLIT_MULTIALLELIC {
        cpus = params.threads
        memory = '4.GB'
    }
    
    withName: GENOTYPE_QC {
        cpus = params.threads
        memory = '4.GB'
    }
    
    withName: VARIANT_QC {
        cpus = params.threads
        memory = '4.GB'
    }
    
    withName: SAMPLE_QC {
        cpus = params.threads
        memory = '8.GB'
    }
    
    withName: ADD_AF {
        cpus = params.threads
        memory = '6.GB'
    }
}

// Executor configuration
executor {
    name = 'local'
    cpus = params.threads * 2
    memory = '16.GB'
}

// Prohibit docker-in-docker
docker.enabled = false
singularity.enabled = false

process.stageInMode  = 'symlink'
process.stageOutMode = 'symlink'

```
**Notes**

- `seq_type` selects which Sample QC thresholds are applied (wes vs wgs).

- metadata_csv is required.

- The fields INFO/QD, INFO/DP, INFO/MQ, INFO/FS, INFO/ReadPosRankSum,  FORMAT/GQ, FORMAT/GT and FORMAT/AD must be correctly described in the header for the filtering to work correctly.

### OPTION 1: After clonning the repository 

#### Basic run: 

Set up the correct input/output paths in `nextflow.config` and then:

```bash
nextflow run main.nf
```

#### With custom parameters

Provide the input/output path as custom parameters: 

```bash
nextflow run main.nf \
  --input "/data/vcfs/" \
  --metadata_csv "/data/metadata.csv" \
  --sceVCF_path "/tools/sceVCF" \
  --seq_type "WGS" \
  --threads 8
```

**NOTE:** This will overwrite the configuration parameters set in nextflow.config


### OPTION 2: From a Docker Container

1. Move the VCF and the CSV to the input folder. Do not modify the paths in the nextflow.config

2. Run:
```
docker compose up
```

When using the Docker set up, please ensure the paths in `nextflow.config` are: 

``` 
params {
    input = "/data/input" // move input VCFs to this folder.  
    output = "/data/output" // output directory
    sceVCF_path = "/usr/local/bin/sceVCF" // File downloaded by clonning. 
    metadata_csv =  "/data/input/metadata_file.csv" // header: sample/sex/ancestry. 
```
## Outputs 

After a successful execution, you'll find inside the */work* folder:

**Intermediate outputs**

- Indexed VCFs (.tbi)

- Split multiallelic VCFs

- Masked genotypes VCFs

- Variant-filtered and masked genotypes VCFs

- Variant and Sample-filtered with masked genotypes VCFs

- sample-qc-fails.tsv file with the samples deleted and which filter step failed

**Final output**

In the folder /output you'll find: 

* *input_vcf_baseName*.vcf.gz  → fully filtered VCF with allele frequency annotations and witout sample level information. 

- *input_vcf_baseName*.vcf.gz.tbi → index for the final VCF.

## Bonus Tracks


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

The intermediate files are saved inside the nextflow folder for each step. 

Go to /work/bc/56f27a[-->] and you'll find the VCF with the multiallelic variants splitted.

Go to /work/33/a6c930[-->] and you'll find all the intermediate files created during the sample QC.

Go to /work/c7/065171[-->] for the final files of the pipeline. 

**Recomendation:** Use at least the argument `-with-report report.html`, this will always show you the folders where the files have been saved.



 
---
💡 **HTML report and timeline**

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


---
💡 **How to generate metadata for the samples**

Three categories of metadata are required per sample to accurately calculate allele frequencies: sex, ancestry, and relatedness. If you don't have this metadata available there are several tools that can infer it using the genomic data. 

Two examples that have been tested and have good accuracy are [GRAF tools](https://github.com/ncbi/graf) and [Hail](https://hail.is/docs/0.2/index.html). 

### Graf tools
GRAF provides functions to infer sex (`graf sex`, PLINK input), detect related samples (`graf rel`, PLINK input), and assign ancestry ([graf anc](https://github.com/jimmy-penn/grafanc/tree/master), VCF input). Note that, for ancestry inference, only super-population calls (e.g. European, Asian, African...) are recommended, since finer-level predictions are not sufficiently accurate. 

⚠️ If you choose to use Graf Tools, please use  [SetRsIdsInBimFile_EGAversion.pl](https://github.com/EGA-archive/AF_bcftools_pipeline/blob/main/SetRsIdsInBimFile_EGAversion.pl) script to perform the annotation step, as we have fixed several bugs present in the original version.

### Hail
Hail offers similar capabilities with `impute_sex` for sex inference and `king` and `pc_relate` for [relatedness filtering](https://hail.is/docs/0.2/guides/genetics.html#remove-related-individuals-from-a-dataset). Another interesting feature is `maximal_independent_set` which outputs the largest subset of unrelated individuals from a dataset.

## References

Lu W, Gauthier LD, Poterba T, Giacopuzzi E, Goodrich JK, Stevens CR, King D, Daly MJ, Neale BM, Karczewski KJ. CHARR efficiently estimates contamination from DNA sequencing data. bioRxiv [Preprint]. 2023 Jun 28:2023.06.28.545801. doi: 10.1101/2023.06.28.545801. Update in: Am J Hum Genet. 2023 Dec 7;110(12):2068-2076. doi: 10.1016/j.ajhg.2023.10.011. PMID: 37425834; PMCID: PMC10327099.

NCBI. (2024). GRAF: Genetic relationship and fingerprinting [Software]. GitHub. https://github.com/ncbi/graf
