require(tidyverse)

args = commandArgs(trailingOnly = T)
manifest = args[1]

# Usage statement
if(length(args) != 1){
    cat("Usage: Rscript combine_sets_new.R <manifest_file> <out_file> \nUse -h, --help for additional help.")
    stop()
}
if(args[1] == "-h" | args[1] == "--help"){
    cat("Usage: Rscript combined_sets_new.R <manifest_file> <out_file> \nManifest file must contain two tab-delimited columns with no header. The first column must contain the cohort name and the second contains the full path to the variant counts VCF and results directory. The counts file is expected to be in VCF format and contains variant counts from across all cohorts.")
    stop()
}

manifest_dat = read_delim(manifest, delim = "\t", col_names = F)
cohorts = manifest_dat$X1

family_ids = c()
counts = data.frame()
for (cohort in cohorts) {
	cohort_families = scan(paste(manifest_dat$X2[which(manifest_dat$X1 == cohort)], "family_ids.txt", sep = ""), "c")
	family_ids = c(family_ids, cohort_families)
	
	count_file = paste(manifest_dat$X2[which(manifest_dat$X1 == cohort)], "variant_counts.vcf.gz", sep = "")
	cohort_counts = read_delim(count_file, delim = "\t", comment = "#", col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")) 
	
	tmp = cohort_counts %>%
    separate(INFO, c("POP_GT_COUNT", NA), sep = ";") %>%
	separate(POP_GT_COUNT, c(NA, "het_count", "hom_count"), sep = ",", convert = T) %>%
	bind_rows(counts, .) %>%
	group_by(CHROM, POS, ID, REF, ALT, QUAL, FILTER) %>%
	summarize(het_count = sum(het_count), hom_count = sum(hom_count))
	
	counts = tmp
}

n_alleles = length(family_ids) * 4

out = counts %>% 
	mutate(
		ref_count = n_alleles - het_count - hom_count,
		INFO = paste("POP_GT_COUNT=", ref_count, ",", het_count, ",", hom_count, sep = "")
	) %>%
	select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

cat(c("##fileFormat=VCFv4.1", "##fileDate=", "##source=generate_private_list.R", "##INFO=<ID=POP_GT_COUNT,Number=3,Type=Integer,Description=Parental population genotype counts ordered by hom REF, het ALT, hom ALT>", "##INFO=<ID=FAMID,Number=1,Type=String,Description=Family ID for private variants>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), sep = "\n", file = "variant_counts.vcf")
write_delim(out[c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")], path = "variant_counts.vcf", delim = "\t", col_names = F, append = T)

write.table(family_ids, "family_ids.txt", sep = "\n", col.names = F, row.names = F, quote = F)
