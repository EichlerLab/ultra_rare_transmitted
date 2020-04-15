args = commandArgs(trailingOnly = T)
file = args[1]
cohort = args[2]
counts = args[3]
het_ac = as.numeric(args[4])
hom_ac = as.numeric(args[5])

# Usage statement
if(length(args) != 5){
    cat("Usage: Rscript combine_sets_new.R <manifest_file> <out_file> <counts_file> <het_ac> <hom_ac> \nUse -h, --help for additional help.")
    stop()
}else{
    if(args[1] == "-h" | args[1] == "--help"){
        cat("Usage: Rscript combined_sets_new.R <manifest_file> <out_file> <counts_file> <het_ac> <hom_ac> \nManifest file must contain three tab-delimited columns and no header. The first column must contain the cohort name, second contains the full path to the variant counts VCF and results directory, and the thrid contains the name (with full path) to the pedigree file). The counts file is expected to be in VCF format and contains variant counts from across all cohorts.")
        stop()
    }
}


# LOAD PAKAGES & MANIFEST
require(tidyverse)

manifest = read.table(file, sep = "\t", stringsAsFactors = F)
variant_counts = read_delim(counts, delim = "\t", comment = "#", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))

variant_list = variant_counts %>% 
    separate(INFO, c("POP_GT_COUNT", NA), sep = ";") %>%
    separate(POP_GT_COUNT, c(NA, "HET", "HOM_ALT"), sep = ",", convert = T) %>% 
    filter(HET <= het_ac & HOM_ALT <= hom_ac) %>%
    mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
    select(HET, HOM_ALT, SNP_ID)

rm(variant_counts)
gc()

# CREATE COMBINED PRIVATE DATASET
dat_file = paste(manifest$V2[which(manifest$V1 == cohort)], "results/rare/all_families_rare.snpeff.txt.gz", sep = "")
dat = read_delim(dat_file, delim = "\t", col_names = T, col_type = "cicccnccnccccciiicccccc")

# EXTRACT ULTRA-RARE VARIANT
out = dat %>% 
  mutate(
    CHROM = ifelse(grepl("chr", CHROM), CHROM, paste("chr", CHROM, sep = "")),
    SNP_ID = paste(CHROM, POS, REF, ALT, sep = ":"),
    LOC_ID = paste(CHROM, POS, sep = "_")
  ) %>%
  inner_join(., variant_list, by = c("SNP_ID"))

# output
write.table(out, paste("results/rare/", cohort, "_rare.snpeff.txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
