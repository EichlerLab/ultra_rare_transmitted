args = commandArgs(trailingOnly = T)
batch = as.numeric(args[1])

require(dplyr)
require(readr)
require(parallel)

genomes = dir("snvs", ".combined.nomulti.vcf.gz$", full.names = T)
genomes = genomes[grep("12051", genomes, invert = T)]

extract_snvs = function(file){
    famid = gsub("snvs/", "", gsub("[a-z.]*.vcf.gz", "", file))
    raw_exome = suppressWarnings(read_delim(file, delim = "\t", comment = "#", col_types = "ci_ccn_c_cc__", col_names = c("X.CHROM", "POS", "REF", "ALT", "QUAL", "INFO", "FA", "MO"), progress = F))

    filtered = raw_exome[which(grepl(";DP=[1-9][0-9];|;DP=[1-9][0-9][0-9]", raw_exome$INFO) & raw_exome$QUAL > 50 & !grepl("^AC=0;|;AC=0;", raw_exome$INFO) & raw_exome$X.CHROM %in% c(1:22, "X", "Y", paste("chr", c(1:22, "X", "Y", "M"), sep = ""), "MT")),]
    rm(raw_exome)
    filtered$n_het = grepl("^1/0|^0/1", filtered$FA) + grepl("^1/0|^0/1", filtered$MO)
    fa_x = filtered$FA[which(filtered$X.CHROM %in% c("X", "chrX"))]
    mo_x = filtered$MO[which(filtered$X.CHROM %in% c("X", "chrX"))]
    fa_y = filtered$FA[which(filtered$X.CHROM %in% c("Y", "chrY"))]
    filtered$n_het[which(filtered$X.CHROM %in% c("X", "chrX"))] = grepl("^1/0|^0/1|^1/1", fa_x) + grepl("^1/0|^0/1", mo_x)
    filtered$n_het[which(filtered$X.CHROM %in% c("Y", "chrY"))] = 0 + grepl("^1/0|^0/1|^1/1", fa_y)
    filtered$n_hom = grepl("^1/1", filtered$FA) + grepl("^1/1", filtered$MO)
    filtered$n_hom[which(filtered$X.CHROM %in% c("X", "chrX"))] = 0 + grepl("^1/1", mo_x)
    filtered$n_hom[which(filtered$X.CHROM %in% c("Y", "chrY"))] = 0
    filtered$famid = famid

    out = filtered[which(filtered$n_het > 0 | filtered$n_hom > 0),]
    rm(filtered)
    info = paste(famid, nrow(out), sep = "\t")
    cat(info, sep = "\n")
    return(out[, c("X.CHROM", "POS", "REF", "ALT", "n_het", "n_hom", "famid")])
}

batch_size = 150
start = seq(1, length(genomes), by = batch_size)
end = seq(batch_size, length(genomes), by = batch_size)
if(length(end) < length(start)){
    end = c(end, length(genomes))
}

cat(paste("looping over genomes in batch ", batch, sep = ""), sep ="\n")
tmp = mclapply(genomes[start[batch]:end[batch]], extract_snvs, mc.cores = 20)
sum(unlist(lapply(tmp, is.null)))
snvs = bind_rows(tmp)
rm(tmp)

cat("counting genotypes by variant", sep ="\n")
tmp.tab = snvs %>% group_by_at(vars(X.CHROM, POS, REF, ALT)) %>% summarize_at(vars(n_het, n_hom), sum)
fam.tab = snvs %>% group_by_at(vars(X.CHROM, POS, REF, ALT)) %>% summarize_at(vars(famid), paste, collapse = ",")
rm(snvs)
snvs.tab = inner_join(tmp.tab, fam.tab, by = c("X.CHROM", "POS", "REF", "ALT"))
rm(list = c("tmp.tab", "fam.tab"))
snvs.tab$famid[which(snvs.tab$n_het > 1 | snvs.tab$n_hom > 0)] = "."
print(length(unique(snvs.tab$famid[which(snvs.tab$n_het == 1 & snvs.tab$n_hom == 0)])))
        
cat("writing batch to file", sep = "\n")
write_delim(snvs.tab, path = paste("batch", batch, "_variant_counts.txt", sep = ""), delim = "\t", col_names = T)
