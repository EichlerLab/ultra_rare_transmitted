#!/usr/env/bin R
## Extract private variants ##
## load packages ##
require(readr)
require(tidyr)
require(dplyr)
require(parallel)

# command line arguments
args = commandArgs(trailingOnly = T)
ped_file = args[1]
counts_file = args[2]
file = args[3]
het_ac = args[4]
hom_ac = args[5]

## define functions ##
# extract and format rare variants from a family
extract_rare = function(file){
	famid = gsub("inheritance/", "", gsub("[a-z.]*.vcf.gz", "", file))
	raw = suppressWarnings(read_delim(file, delim = "\t", comment = "#", col_types = "cicccnccccc--", col_names = c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "FA", "MO"), progress = F))
	info = c(famid, nrow(raw))
	info_cols = names(raw)[1:8]
	
	pro_sex = paste(unique(ped$V6[which(grepl("p", ped$V2) & ped$V1 == famid)]), collapse = ",")
	if(any(grepl("s", ped$V2) & ped$V1 == famid)){
		sib_sex = paste(unique(ped$V6[which(grepl("s", ped$V2) & ped$V1 == famid)]), collapse =",")
	}else{
		sib_sex = "NA"
	}

	# extract SNVs called by GATK and FreeBayes
	inherited = raw %>% filter(grepl(";DP=[2-9][0-9];|;DP=[1-9][0-9][0-9][0-9]*;", INFO) & QUAL > 50 & grepl("MENDEL=True", INFO) & (grepl("set=Intersect", INFO) | nchar(REF) > 1 & nchar(REF) == nchar(ALT)))
	info = c(info, nrow(inherited))
	rm(raw)
	# extract rare variants
	snv_id = inherited %>% unite(SNVID, c(X.CHROM, POS, REF, ALT), sep = ":")
	rm(inherited)
	family = semi_join(snv_id, rare, by = "SNVID") %>% separate(SNVID, c("X.CHROM", "POS", "REF", "ALT"), sep = ":")
	info = c(info, nrow(family))
	rm(snv_id)
	print(info)
	if(nrow(family) > 0){
		family$CARRIER = ifelse(grepl("^1/0|^0/1", family$FA) & !family$X.CHROM %in% c("chrX", "chrY"), "fa", ifelse(family$X.CHROM %in% c("chrX", "chrY") & grepl("^1/1", family$FA), "fa", ifelse(family$X.CHROM != "chrY" & grepl("^1/0|^0/1", family$MO), "mo", NA)))
		family$FAMILY = famid
		family$PROBAND_SEX = pro_sex
		family$SIBLING_SEX = sib_sex
		family$INFO = paste(family$INFO, ";CARRIER=", family$CARRIER, ";PROBAND_SEX=", family$PROBAND_SEX, ";SIBLING_SEX=", family$SIBLING_SEX, ";FAMID=", family$FAMILY, sep = "")
		write.table(family[which(is.na(family$CARRIER) == F), info_cols], file = paste("results/rare/families/", famid, ".txt", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F)
		rm(family)
		gc()
	}
	return("done")
}

## initalize data ##
# read pedigree file
ped = read.table(ped_file, sep =" ", stringsAsFactors = F)

# calculate minor allele count
families = scan("family_ids.txt", "c")
n_parents = length(families) * 2

print(c(het_ac, hom_ac, n_parents))
# read variant counts VCF
print("reading variant counts")
counts = suppressWarnings(read_delim(counts_file, delim = "\t", comment = "#", progress = F, col_names = F))
rare = counts %>% 
  separate(X8, c("POP_GT_COUNT", NA), sep = ";") %>% 
  mutate(POP_GT_COUNT = gsub("POP_GT_COUNT=", "", POP_GT_COUNT)) %>% 
  separate(POP_GT_COUNT, c("hom_ref", "het", "hom_alt"), sep = ",", convert = T) %>% 
  filter(het <= het_ac & hom_alt <= hom_ac) %>%
  transmute(SNVID = paste(X1, X2, X4, X5, sep = ":"))

rm(counts)
## load file names ##
files = dir("inheritance", ".inheritance.vcf.gz$", full.names = T)

# load header
n = as.numeric(system(paste("zcat ", files[1], " | head -n 7000 | grep -c '#'", sep = ""), intern = T)) - 1
header = read.table(files[1], sep = "\n", comment.char = "@", nrow = n)

info = as.character(header[grep("INFO=|FILTER=", header$V1),])
new_info = gsub(">", '">', gsub("Description=", 'Description="', info))
new_info = c(new_info, '##INFO=<ID=FAMID,Number=1,Type=String,Description="Family ID">', '##INFO=<ID=CARRIER,Number=1,Type=String,Description="Carrier parent">', '##INFO=<ID=PROBAND_SEX,Number=1,Type=String,Description="Proband sex">', '##INFO=<ID=SIBLING_SEX,Number=1,Type=String,Description="Sibling sex">')

tmp = as.character(header[grep("FORMAT|INFO", header$V1, invert = T),])
new_header = c(tmp[1:2], sort(new_info), tmp[3:length(tmp)])

# output new header
print("outputting VCF header")
famid = gsub("inheritance/", "", gsub("[a-z.]*.vcf.gz", "", file))
write.table(new_header, file = paste("results/rare/", famid, ".header", sep = ""), sep = "\n", col.names = F, quote = F, row.names = F)
cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", sep = "\n", file = paste("results/rare/", famid, ".header", sep = ""), append = T)

# extract rare variants from each family (runs 10 families in parallel)
print("extracting rare variants")
extract_rare(file)
