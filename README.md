# Ultra-rare transmitted variant pipeline
Pipeline to generate ultra-rare transmitted callsets for genomes and exomes

## Overview: 
This pipeline starts with left-normalized and reheadered VCFs from FreeBayes and GATK (see genome workflow for details on generating these VCFs). It will:
1. Merge these two variant callsets together (per family), 
2. Identify the parent of origin and recipient children (if any) (per family), 
3. Count alleles across parents by genotype (e.g. HOM REF, HET, HOM ALT) and output into a VCF (criteria detailed below) (across families/per cohort),
4. Extract ultra-rare variants (criteria detailed below) (per family),
5. Annotate variants with SNPEFF (Gene name, transcript IDs, variant effect, etc.) (per family),
6. Concatenate annotated ultra-rare variants (across families/per cohort),
7. "Filter" sites against recent repeats, gaps, centromeres, PARs, and LCRs
8. Extract protein-coding sites

## Included scripts:
- `inheritance_exome.snake`
- `inheritance_genome.snake`
- `phase.json`
- `scripts/inheritance.py`
- `scripts/exomes/batch_private_list.R`
- `scripts/genomes/batch_private_list.R`
- `scripts/gather_private_list_batches.R`
- `scripts/genomes/extract_rare.R`
- `scripts/exomes/extract_rare.R`
- `scripts/extract_coding_rare.R`
- `scripts/combine_sets_new.R`

This pipeline assumes a file called `family_ids.txt` is present in the directory the snakefile is called. This file should contain the list of unique fmailies to be included in the analysis.

**NOTES:**
This pipeline has been fully tested on centOS6. To run on centOS7 the versions of vcflib, bedtools, and htslib specified in the shell prefix will need to be updated to the following: `vcflib/202002`, `bedtools/2.29.0`, and `htslib/1.9-20`

## Example
Run this command within the example directory on a high memory server (e.g. lynx, ocelot, or a qlogin with at least 300GB of memory). Do ensure you have the `DRMAA_LIBRARY_PATH` specified either in your `.bash_profile` or prior to executing the command.
`snakemake -j 100 --jobname "{rulename}.{jobid}" --drmaa "-V -cwd -e ./log -o ./log {params.sge_opts} -w n -S /bin/bash" -k -w 120 --rerun-incomplete -s Snakefile -p`

## Details:
### **Step 1: Merge variant callsets together**
This step merges the left-normalized and reheadered per family VCFs generated by FreeBayes and GATK using GATK-3.8 CombineVariants. To generate these files refer to the genome workflow github. No additional filters are applied at this step, however multiallelic sites are split and an additional info flag called "set" is be added to the VCF during the merge. This "set" flag will be used to select for high quality sites in downstream steps.

Input is two VCFs. Output is one VCF.

Dependencies: GATK3.8, java8, genome reference (reference genomes for gh19 and hg38 can be found here: `/net/eichler/vol27/projects/autism_inheritance/nobackups/reference/`)

**NOTES:** The CombineVariants utility has been depreciated/removed in newer versions of GATK, so this step is version sensitive.

### **Step 2: Identify parent-of-origin and recipient children**
**THIS SCRIPT IS WRITTEN IN PYTHON2**
This step will annotate the per family merged FreeBayes-GATK VCFs with the following information: 
1. sites are checked for Mendelian inheritance
2. parent-of-origin for informative sites
3. all children are checked for transmission - children that recieve the informative event from parents are noted.
Two INFO flags are added to the VCF from this step "MENDEL" and "INHER". The "MENDEL" flag has a binary True/False value indicating whether the site follows stric rules of Mendelian inheritance. The "INHER" flag contains information regarding transmission of the event to children and is formatted as "parent-of-origin"-"recipient children" (e.g. fa-p1 or mo-p1,s1). 

Input and output are one VCF

Dependencies: python2, inheritance.py

**NOTES:** 
Sex chromosome handling has been built-in as follows: any site with a variant call (e.g. 1/., 0/1, or 1/1) in males on the X chromosome is considered variant, all Y chromosome sites on the female are ignored, X chromosome sites in females are treated like an autosome

This script has not been tested on extended families, so it is unclear how it will perform. Adjustments may need to be made to accomodate family structures beyond mother, father, children.

This script will produce an empty file for families missing a parent or if it runs incorrectly. This will cause a problem with down stream rules

### **Step 3: Count alleles across parents**
**THIS IS SPECIFIED AS A LOCAL RULE**
This step is precariously stable (meaning it requires ~350-400GB of RAM to run for a cohort of 2300 families, larger cohorts have not been tested). Due to the nature of this step a large amount of memory is required and should be run on a high memory node (i.e. lynx, ocelot, or qlogin with > 350GB of memory). A batch-and-gather approach is used to count variants across parents within a cohort. Additional functionality to combine variant counts across cohorts is provided in the "Accessory functions" section below. The cohort is processed in batches of 150 families (this number is speicied in `n_batches` parameter in the `phase.json` and can be changed as needed). The batch processing iteratively reads in the family level data from the `combined.nomulti.vcf.gz` files in the `snvs/` directory, filters loci for variant quality (QUAL > 50) and read depth (DP > 10 for genomes and DP > 20 for exomes), determines the parental genotypes (e.g. het, hom, ref) and then combines the genotype counts with the previously processed families in the batch (i.e. genotype counts are added to existing counts for previously observed sites and new sites are added to the growing variant list). Each batch outputs a file with the genotype counts for that batch. The batch genotype count files are then gathered/aggregated together by variant, the HOM ALT genotype counts are multipled by 2 to get the allele counts (2 alleles per HOM genotype), and the the HOM REF allele counts are determined by subtracting the HET and HOM ALT allele counts from the 4n available parental alleles. All observed variants in the cohort will be output into a VCF file which will contain two flags in the INFO field "POP_GT_COUNT" and "FAMID". The "POP_GT_COUNT" contains the HOM REF, HET, and HOM ALT allele counts, respecitvely, for the variant. The "FAMID" field will contain the family ID for the variant only if it is unique to that family (i.e. HET=1). Genotypes have been removed, and QUAL and FILTER information have been set to ".".

Input is all family VCFs from the merge rule (step 1). Output is one VCF called `variant_counts.vcf.gz` as well as intermediate files called `batch{n}_counts.txt`. These intermediate files can be removed when this step is complete.

Dependencies: R/3.5.1 or later, tidyverse, parallel, exomes: exomes/batch_private_list.R, genomes: genomes/batch_private_list.R, gather_private_list_batches.R

**NOTES:** 
Sex chromosome handling has been built-in as follows: any site with a variant call (e.g. 1/., 0/1, or 1/1) in males on the X chromosome is considered variant and counted only once, all Y chromosome sites on the female are ignored, X chromosome sites in females are treated like an autosome.

### **Step 4: Extract ultra-rare variants**
This step is run per family and utilizes the VCF files generated in steps 2 and 3 to extract the ultra-rare variants observed in each family. Ultra-rare variants must meet the following criteria: 
	1) Observed in both the FreeBayes and GATK callset (FreeBayes genotypes are prioritized when genotypes between the callers are discordant) OR is an MNV
	2) QUAL > 50
	3) Variant follows rules of strict Mendelian inheritance
	4) Variant has a read depth > 10 for genome and > 20 for exomes
This step will add the "CARRIER" INFO flag to a family VCF which contains only ultra-rare variants. Ultra-rare is defined as parent HET < 10 and parent HOM ALT == 0.

This is the step where the allele count cutoff is used to extracted only alleles with parent HET and HOM ALT allele counts less than or equal to the user specified threshold. This can be changed in the `phase.json` for the `max_het_ac` and `max_hom_ac` parameters. For the anlayses in Wilfert et al. the values used were `"max_het_ac": 9` and `"max_hom_ac": 0`.

Dependencies: R, tidyverse

### **Step 5: Annotate**
This step will annotate the ultra-rare varint family VCF with gene names, variant effect, etc. using SNPEFF. It will then annotate with CADD, ExAC allele frequencies and dbSNP - these VCFs need to be provided inthe config file. The path to SNPEFF and the genome build also need to be included in the config file. This step utilizes two SNPEFF accessory functions. The first splits SNPEFF annotations into separate lines per annotation if multiple are present. The second converts the SNPEFF annotations from a VCF format to a tab-delimited text file. **You will need to specify you genome build and version here. The json is currently set to hg38.**

Dependencies: SNPEFF, java8, CADDv1.4, ExAC, dbSNP

**NOTES:** CADDv1.4, ExAC, dbSNP VCFs will need to be downloaded/present (these are currently located here for both hg19 and hg38: `/net/eichler/vol27/projects/autism_inheritance/nobackups/db/`). Make sure you have alread cached the genome data necessary to run SNPEFF (refer to SNPEFF documentation if using for the first time).

Step 6: Concatenate per family ultra-rare variants
Per family SNPEFF text files are concatenated and a header is added.

### **Step 7: "Filter"**
BEDTools is used to create a file of "clean" positions. This file is used to select protein-coding sites which are not located within a recent repeat, centromere, gap, PAR, or LCR. Paths to these files are currently hardcoded and should be added to the config if they ever move or change. BEDTools 2.24.0 was used here but newer versions should also be fine.

Dependencies: BEDTools/2.24.0 or later, centromere and gaps, LCRs, recent repeats, PARs (these data are located here for hg38: `/net/eichler/vol27/projects/autism_inheritance/nobackups/reference/`)

**NOTES:** 
This step does not explicity remove any sites from the input file, but rather creates a BED file that contains "clean" variant calls. This BED file is used in the next step when extracting protein-coding varaints. 

The centromere, gap, LCR, PAR, and recent repeat files are currently for hg38 and are hardcoded into the pipeline.

### **Step 8: Extract protein-coding variants**
This step loads the ultra-rare callset into memory (run this on lynx/ocelot/high memory node), removes low quality/problem regions, extracts protein-coding sites according to standard Eichler lab criteria (see below), and parses transmitted sites by person (e.g. a site that is transmitted will be observed three times - once for each child). This will create four additional columns: var_type, transmitted, person, carrier. The var_type column indicates whether the site is lgd, mis or syn according to the Eichler lab rule. The transmitted column indicates to which child the site was transmitted (e.g. p1, s1, p2, etc.). If the site was not transmitted you will observe a "." or NA. the carrier column will indicate which parent is the carrier - this column should match with the CARRIER column. The person column will contain the sample ID for the child carrying the allele. Sites with no carrier will have a person column value of `{FAMID}.NA`.

Dependencies: R/3.5.1 or later, tidyverse

Additional utilities:
Combining ultra-rare callsets across cohorts
