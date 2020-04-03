require(tidyverse)
args = commandArgs(trailingOnly = T)
variants = args[1]
counts_vcf = args[2]
filter = args[3]
outfile = args[4]

vcf = read_delim(counts_vcf, comment = "#", delim = "\t", col_names = c("CHROM", "POS", "REF", "ALT", "INFO"), col_types = "ci-cc--c")
allele_counts = vcf %>% 
   mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = ":")) %>% 
   separate(INFO, c("POP_GT_COUNT", NA), sep = ";") %>% 
   separate(POP_GT_COUNT, c(NA, "HET", "HOM_ALT"), sep = ",", convert = T)

filtered = read_delim(filter, delim = "\t", col_names = c("CHROM", "START", "END"))
dat = read_delim(variants, delim = "\t", col_names = T)
rare_var = dat %>% 
   separate(INHER, c("carrier", "transmitted"), sep = "-", drop = F) %>% 
   separate_rows(transmitted, sep = ",") %>% 
   mutate(SNP_ID = paste(CHROM, POS, REF, ALT, sep = ":"), id = paste(CHROM, POS, REF, ALT, `EFF[*].GENE`, sep = "_"), person = paste(FAMID, transmitted, sep = ".")) %>% 
   left_join(., allele_counts[c("SNP_ID", "HET", "HOM_ALT")], by = "SNP_ID") %>%
   inner_join(., filtered, by = c("CHROM", "POS" = "END"))

# Extract coding region
MIS_EFF = c("NON_SYNONYMOUS_CODING", "missense_variant")
LGD_EFF = c("SPLICE_SITE_DONOR", "splice_donor_variant", "SPLICE_SITE_ACCEPTOR", "splice_acceptor_variant", "FRAME_SHIFT", "frameshift_variant", "STOP_LOST", "stop_lost", "STOP_GAINED", "stop_gained")
SYN_EFF = c("SYNONYMOUS_CODING", "synonymous_variant")

lgd = rare_var %>% 
   filter(grepl(paste(LGD_EFF, collapse = "|"), `EFF[*].EFFECT`)) %>% 
   mutate(var_type = "lgd")
mis = rare_var %>% 
   filter(grepl(paste(MIS_EFF, collapse = "|"), `EFF[*].EFFECT`) & !id %in% lgd$id) %>% 
   mutate(var_type = "mis")
syn = rare_var %>% 
   filter(grepl(paste(SYN_EFF, collapse = "|"), `EFF[*].EFFECT`) & !id %in% c(lgd$id, mis$id)) %>% 
   mutate(var_type = "syn")

info = bind_rows(list(lgd, mis, syn))

write.table(info, outfile, sep = "\t", col.names = T, row.names = F, quote = F)
