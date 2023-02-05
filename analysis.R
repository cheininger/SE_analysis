# Libraries ####

library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(WebGestaltR)

# Define functions ####

subset_quants <- function(df) {
    df <- df[,c(1:3,7:9)]
    colnames(df) <- c("chr", "start", "end", "strand", "ID", "quant")
    return(df)
}

remove_chr_double <- function(df) {
    names(df)[names(df) == "CHROM"] <- "chr"
    df$chr <- gsub("(chrchr)", "chr", df$chr)
    return(df)
}

filter_chr_col <- function(df) {
    df <- df %>% filter(grepl("^(chr){1}[[:digit:]{1,2}XY]", df$chr))
    return(df)
}

rose_enhancer_id <- function(df) {
    df <- df %>% mutate('REGION_ID' = paste0("PROK_LIV_ROSE_", chr, "_", START, "_", STOP))
    return(df)
}

enhancer_quant_id <- function(df) {
    df <- df %>% mutate('quant_ID' = paste0("PROK_LIV_", 
                                            select(df, matches("^enhancer.seqnames$")),
                                            "_",
                                            select(df, matches("^enhancer.start$")),
                                            "_",
                                            select(df, matches("^enhancer.end$"))))
    return(df)
}

sort_df <- function(df) {
    df[order(df$chr),]
}


# First, read all ROSE output files into a list and make GRanges ####

rose_out_list <- list.files(path = "./data", pattern = "*_AllStitched.table.txt", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE, skip = 5, col_types = c("cciiiidii")) 


rose_out_list <- lapply(rose_out_list, remove_chr_double)           # some entries in the chromosome name contain "chrchr"
                                                                    # this function deals with that

rose_out_list <- lapply(rose_out_list, filter_chr_col)              # filter so that only chromosomes 1-19 and X and Y are included 

rose_out_list <- lapply(rose_out_list, rose_enhancer_id)            # generate ID for each stitched region from ROSE

rose_out_list <- lapply(rose_out_list, sort_df)                     # sort by chromosome name

rose_out_gr <- lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)


rose_se_extended <- lapply(rose_out_list, function(df) { df$START <- as.integer(df$START - 50000)             
return(df)}) %>% lapply(., function(df) { df$STOP <- as.integer(df$STOP + 50000)                                # extend start and stop positions by 50k
return(df)}) %>% lapply(., function(df) { df$START[df$START < 0] <- 0                                           # set negative values for start to 0
df$START <- as.integer(df$START)
return(df)}) %>% lapply(., function(df) { df <- df[c("chr", "START", "STOP", "REGION_ID", "stitchedPeakRank", "isSuper")]
return(df)}) %>% lapply(., function(df) { df <- filter(df, isSuper == 1)                                        # continue only with super-enhancers
return(df)})


rose_se_extended_gr <- lapply(rose_se_extended, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)


# Read ranges for fasting-activated enhancers ####

fasting_activated_enhancer_list <- list.files(path = "./data/", pattern = "fasting_activated_*", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE)

fasting_activated_enhancer_gr_list <- lapply(fasting_activated_enhancer_list, makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = TRUE)


# Subset ROSE ranges with fasting-activated enhancers ####


se_w_fasting_activated_list <- vector(mode = "list", length = 5)

se_w_fasting_activated_list[[1]] <- as.data.frame(subsetByOverlaps(rose_se_extended_gr$`1h_FAST`, fasting_activated_enhancer_gr_list$`1h`))
se_w_fasting_activated_list[[2]] <- as.data.frame(subsetByOverlaps(rose_se_extended_gr$`3h_FAST`, fasting_activated_enhancer_gr_list$`3h`))
se_w_fasting_activated_list[[3]] <- as.data.frame(subsetByOverlaps(rose_se_extended_gr$`6h_FAST`, fasting_activated_enhancer_gr_list$`6h`))
se_w_fasting_activated_list[[4]] <- as.data.frame(subsetByOverlaps(rose_se_extended_gr$`12h_FAST`, fasting_activated_enhancer_gr_list$`12h`))
se_w_fasting_activated_list[[5]] <- as.data.frame(subsetByOverlaps(rose_se_extended_gr$`24h_FAST`, fasting_activated_enhancer_gr_list$`24h`))

names(se_w_fasting_activated_list) <- c("1h", "3h", "6h", "12h", "24h")


# Generate homerpeak files for motif analysis with HOMER ####

homerpeak_1h <- se_w_fasting_activated_list$`1h`[c("REGION_ID", "seqnames", "start", "end")]
homerpeak_1h <- homerpeak_1h[rep(seq_len(nrow(homerpeak_1h)), each = 2), ]
homerpeak_1h <- homerpeak_1h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_1h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_1h, "./output/se_activated_1h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_3h <- se_w_fasting_activated_list$`3h`[c("REGION_ID", "seqnames", "start", "end")]
homerpeak_3h <- homerpeak_3h[rep(seq_len(nrow(homerpeak_3h)), each = 2), ]
homerpeak_3h <- homerpeak_3h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_3h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_3h, "./output/se_activated_3h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_6h <- se_w_fasting_activated_list$`6h`[c("REGION_ID", "seqnames", "start", "end")]
homerpeak_6h <- homerpeak_6h[rep(seq_len(nrow(homerpeak_6h)), each = 2), ]
homerpeak_6h <- homerpeak_6h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_6h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_6h, "./output/se_activated_6h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_12h <- se_w_fasting_activated_list$`12h`[c("REGION_ID", "seqnames", "start", "end")]
homerpeak_12h <- homerpeak_12h[rep(seq_len(nrow(homerpeak_12h)), each = 2), ]
homerpeak_12h <- homerpeak_12h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_12h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_12h, "./output/se_activated_12h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_24h <- se_w_fasting_activated_list$`24h`[c("REGION_ID", "seqnames", "start", "end")]
homerpeak_24h <- homerpeak_24h[rep(seq_len(nrow(homerpeak_24h)), each = 2), ]
homerpeak_24h <- homerpeak_24h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_24h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_24h, "./output/se_activated_24h.homerpeak", delim = "\t", col_names = FALSE)


# Make GRanges objects from filtered super-enhancers, resetting positions in the process

se_w_activated_extended <- lapply(se_w_fasting_activated_list, function(df) { df$start <- as.integer(df$start + 50000)
return(df)}) %>% lapply(., function(df) { df$end <- as.integer(df$end - 50000)                                  # reset positions for super-enhancers
return(df)}) %>% lapply(., function(df) { df <- df[c("seqnames", "start", "end", "REGION_ID", "stitchedPeakRank", "isSuper")]
return(df)}) %>% lapply(., function(df) { names(df)[names(df) == "seqnames"] <- "chr"
return(df)}) %>% lapply(., makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = TRUE)


# Read gene expression TPM data ####

tpm_matrix <- read_delim("./data/tpm_matrix_gene.txt", delim = "\t", col_names = TRUE)

tpm_matrix <- column_to_rownames(tpm_matrix, var = "gene")

# Get TSS per gene from Biomart ####

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "mmusculus_gene_ensembl",
                      mirror = "www")

gene_tss <- getBM(attributes = c("external_gene_name", 
                                 "ensembl_gene_id",
                                 "chromosome_name",
                                 "start_position",
                                 "end_position",
                                 "transcription_start_site",
                                 "strand",
                                 "gene_biotype"),
                  filters = "external_gene_name",
                  values = rownames(tpm_matrix),
                  mart = ensembl)


is_unique <- function(x) length(unique(x)) == 1

# All transcripts of a gene should have the same start_position, i.e. the 5´
# most recorded position for any of the transcripts
stopifnot(tapply(gene_tss$start_position, gene_tss$ensembl_gene_id, is_unique))

# All transcripts of a gene should have the same end_position, i.e. the 3´
# most recorded position for any of the transcripts
stopifnot(tapply(gene_tss$end_position, gene_tss$ensembl_gene_id, is_unique))

# Therefore start_position should always be less than end_position, no matter
# the strand of the gene
stopifnot(mapply(function(x, y) x < y,
                 gene_tss$start_position,
                 gene_tss$end_position))

# Reduce to most upstream TSS per gene -----------------------------------------

# Detach biomaRt before loading dplyr because they both define `select`
# detach("package:biomaRt")

# Added dplyr:: before select to avoid conflict and detaching biomaRt
suppressPackageStartupMessages(library("dplyr"))

tss_per_gene <- gene_tss %>%
    group_by(ensembl_gene_id, external_gene_name, chromosome_name,
             strand, gene_biotype) %>%
    summarize(tss = if (all(strand == 1)) min(start_position) else max(end_position)) %>%
    ungroup() %>%
    arrange(ensembl_gene_id)

# Format -----------------------------------------------------------------------

tss_out <- tss_per_gene %>%
    mutate(strand = ifelse(strand == 1, "+", "-")) %>%
    dplyr::select(gene_id = ensembl_gene_id,
           gene_name = external_gene_name,
           chr = chromosome_name,
           tss,
           strand,
           biotype = gene_biotype)

tss_out$tss_1 <- tss_out$tss + 1                        # add "end" position for TSS to be able to make GRanges object 

tss_out$chr <- paste0("chr", tss_out$chr)

tss_out <- tss_out[c("chr", "tss", "tss_1", "gene_name", "biotype")]

tss_gr <- makeGRangesFromDataFrame(tss_out, seqnames.field = "chr", start.field = "tss", end.field = "tss_1",
                                   strand.field = "strand", keep.extra.columns = TRUE, ignore.strand = FALSE)


# Find genes within 50 kbp of SE by overlapping tss_gr with extended rose output ####


tss_se_anno_list <- list(
    '1h' = as.data.frame(mergeByOverlaps(tss_gr, se_w_activated_extended$`1h`)),
    '3h' = as.data.frame(mergeByOverlaps(tss_gr, se_w_activated_extended$`3h`)),
    '6h' = as.data.frame(mergeByOverlaps(tss_gr, se_w_activated_extended$`6h`)),
    '12h' = as.data.frame(mergeByOverlaps(tss_gr, se_w_activated_extended$`12h`)),
    '24h' = as.data.frame(mergeByOverlaps(tss_gr, se_w_activated_extended$`24h`))
    )


tss_se_anno_list <- lapply(tss_se_anno_list, function(df) {
    names(df) <- c("tss_chr", "tss_start", "tss_end", "tss_width", "tss_strand",
                   "tss_gene_name", "tss_biotype", "gene_name", "biotype",
                   "se_chr", "se_start", "se_end", "se_width", "se_strand", 
                   "se_ID", "se_stitchedPeakRank", "se_isSuper",
                   "ID", "stitchedPeakRank", "isSuper")
    return(df)}) %>% lapply(., function(df) {
        df <- df[c("gene_name", "se_chr", "se_start", "se_end",
                   "se_ID", "se_stitchedPeakRank", "se_isSuper")]
        return(df)}) %>% lapply(., function(df) { df$se_start <- as.integer(df$se_start + 50000)
        return(df)}) %>% lapply(., function(df) { df$se_end <- as.integer(df$se_end - 50000)
        return(df)})



# Add TPM information to TSS annotated SEs ####

samples <- names(tss_se_anno_list)


for (sample in samples) {
    
    for (i in 1:nrow(tss_se_anno_list[[sample]])) {
        
        gene <- tss_se_anno_list[[sample]]$gene_name
        
        tss_se_anno_list[[sample]]$TPM <- tpm_matrix[gene, paste0(sample, "_FAST")]
        
    }
    
}


tss_se_anno_list <- lapply(tss_se_anno_list, function(df) {
    names(df) <- c("gene_name", "chr", "start", "end", "ID", "stitchedPeakRank",
                   "isSuper", "TPM")
    return(df)
}) %>% lapply(., function(df) {
    df <- filter(df, TPM >= 5)
    return(df)
}) %>% lapply(., function(df) {
    df <- df %>% distinct(gene_name, .keep_all = TRUE)
    return(df)
})


write_delim(as.data.frame(tss_se_anno_list$`1h`$gene_name), "./output/annotated_genes_se_1h.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(tss_se_anno_list$`3h`$gene_name), "./output/annotated_genes_se_3h.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(tss_se_anno_list$`6h`$gene_name), "./output/annotated_genes_se_6h.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(tss_se_anno_list$`12h`$gene_name), "./output/annotated_genes_se_12h.txt", delim = "\t", col_names = FALSE)
write_delim(as.data.frame(tss_se_anno_list$`24h`$gene_name), "./output/annotated_genes_se_24h.txt", delim = "\t", col_names = FALSE)


# Run WebGestalt for each set of genes ####

WebGestaltR(enrichMethod = "ORA", organism = "mmusculus", interestGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process_noRedundant",
            interestGeneFile = "./output/annotated_genes_se_1h.txt", referenceSet = "genome",
            sigMethod = "fdr", outputDirectory = "./output/", isOutput = TRUE, projectName = "genes_1_up")

WebGestaltR(enrichMethod = "ORA", organism = "mmusculus", interestGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process_noRedundant",
            interestGeneFile = "./output/annotated_genes_se_3h.txt", referenceSet = "genome",
            sigMethod = "fdr", outputDirectory = "./output/", isOutput = TRUE, projectName = "genes_3_up")

WebGestaltR(enrichMethod = "ORA", organism = "mmusculus", interestGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process_noRedundant",
            interestGeneFile = "./output/annotated_genes_se_6h.txt", referenceSet = "genome",
            sigMethod = "fdr", outputDirectory = "./output/", isOutput = TRUE, projectName = "genes_6_up")

WebGestaltR(enrichMethod = "ORA", organism = "mmusculus", interestGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process_noRedundant",
            interestGeneFile = "./output/annotated_genes_se_12h.txt", referenceSet = "genome",
            sigMethod = "fdr", outputDirectory = "./output/", isOutput = TRUE, projectName = "genes_12_up")

WebGestaltR(enrichMethod = "ORA", organism = "mmusculus", interestGeneType = "genesymbol",
            enrichDatabase = "geneontology_Biological_Process_noRedundant",
            interestGeneFile = "./output/annotated_genes_se_24h.txt", referenceSet = "genome",
            sigMethod = "fdr", outputDirectory = "./output/", isOutput = TRUE, projectName = "genes_24_up")




# tss_se_anno_gr <- lapply(tss_se_anno_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)
















