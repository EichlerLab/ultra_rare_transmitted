require(readr)
require(dplyr)

files = dir(".", "batch[0-9]*_variant_counts.txt")

snvs.tab = tibble()
for(file in files){
    cat(file, sep = "\n", file = stderr())
    dat = suppressWarnings(read_delim(file, delim = "\t", comment = "#", col_types = "cicciic", col_names = c("X.CHROM", "POS", "REF", "ALT", "n_het", "n_hom", "famid"), progress = T))
    tmp = bind_rows(snvs.tab, dat)
    rm(dat)
    
    cat("counting genotypes by variant", sep = "\n", file = stderr())
    tmp.tab = tmp %>% group_by_at(vars(X.CHROM, POS, REF, ALT)) %>% summarize_at(vars(n_het, n_hom), sum)
    fam.tab = tmp %>% group_by_at(vars(X.CHROM, POS, REF, ALT)) %>% summarize_at(vars(famid), paste, collapse = ",")
    rm(tmp)
    snvs.tab = ungroup(inner_join(tmp.tab, fam.tab, by = c("X.CHROM", "POS", "REF", "ALT")))
    rm(list = c("tmp.tab", "fam.tab"))
    snvs.tab$famid[which(snvs.tab$n_het > 1 | snvs.tab$n_hom > 0)] = "."
    
    cat(length(unique(snvs.tab$famid[which(snvs.tab$n_het == 1 & snvs.tab$n_hom == 0)])), sep = "\n", file = stderr())
}

n_genomes = as.numeric(system("ls snvs/*.vcf.gz | wc -l", intern = T))

pop_gt_info = paste("POP_GT_COUNT=", paste((n_genomes * 4) - snvs.tab$n_het - (snvs.tab$n_hom*2), snvs.tab$n_het, snvs.tab$n_hom*2, sep = ","), sep = "")
famid_info = paste("FAMID=", snvs.tab$famid, sep = "")

cat("creating VCF", sep = "\n", file = stderr())
snvs.tab$QUAL = "."
snvs.tab$ID = "."
snvs.tab$FILTER = "."
snvs.tab$INFO = paste(pop_gt_info, famid_info, sep = ";")

cat("outputting to file", sep = "\n", file = stderr())
cat(c("##fileFormat=VCFv4.1", "##fileDate=", "##source=generate_private_list.R", "##INFO=<ID=POP_GT_COUNT,Number=3,Type=Integer,Description=Parental population genotype counts ordered by hom REF, het ALT, hom ALT>", "##INFO=<ID=FAMID,Number=1,Type=String,Description=Family ID for private variants>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), sep = "\n", file = "variant_counts.vcf")
write_delim(snvs.tab[c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")], path = "variant_counts.vcf", delim = "\t", col_names = F, append = T)
