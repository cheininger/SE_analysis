# Libraries ####

library(tidyverse)
library(GenomicRanges)
library(biomaRt)

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

extend_range <- function(mode, position, extend_by = 50000) {
    
    if (mode == "start") {
        
        if((position - extend_by) < 0) {
            start_extended <- 0
        } else {
            start_extended <- position - extend_by
        }
        
        return(start_extended)
    } else if (mode == "end") {
        
        end_extended <- position + extend_by
        
        return(end_extended)
    } else {
        
        stop("mode must be one of either 'start' or 'end'")
        
    }
    
}

# first read all ROSE output files into a list and make GRanges ####

rose_out_list <- list.files(path = "./data", pattern = "*_AllStitched.table.txt", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE, skip = 5, col_types = c("cciiiidii")) 


rose_out_list <- lapply(rose_out_list, remove_chr_double)

rose_out_list <- lapply(rose_out_list, filter_chr_col)

rose_out_list <- lapply(rose_out_list, rose_enhancer_id)

rose_out_list <- lapply(rose_out_list, sort_df)


rose_extended_list <- lapply(rose_out_list, function(df) { df$START <- as.integer(df$START - 50000) 
return(df)}) %>% lapply(., function(df) { df$STOP <- as.integer(df$STOP + 50000)
return(df)}) %>% lapply(., function(df) { df$START[df$START < 0] <- 0
df$START <- as.integer(df$START)
return(df)}) %>% lapply(., function(df) { df <- df[c("chr", "START", "STOP", "REGION_ID", "stitchedPeakRank", "isSuper")]})

rose_se_extended <- lapply(rose_extended_list, function(df) { df_sub <- filter(df, isSuper == 1)
return(df_sub)})


rose_extended_gr <- lapply(rose_extended_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)

rose_se_extended_gr <- lapply(rose_se_extended, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)



# read enhancer quantification data into a list and make GRanges ####

enhancer_quant_list <- list.files(path = "./data", pattern = "enhancer_quant_*", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = FALSE)

enhancer_quant_list <- lapply(enhancer_quant_list, subset_quants)

enhancer_quant_list <- lapply(enhancer_quant_list, filter_chr_col)


enhancer_quant_gr_list <- lapply(enhancer_quant_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)



# Merge all ranges per timepoint by overlap ####

merged_0h <- as.data.frame(mergeByOverlaps(enhancer_quant_gr_list[[1]], rose_gr_list[[1]]))


merged_list <- list(
    '0h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`0h_AL`, rose_gr_list$`0h_AL`)),
    '1h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`1h_AL`, rose_gr_list$`1h_AL`)),
    '1h_FAST' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`1h_FAST`, rose_gr_list$`1h_FAST`)),
    '3h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`3h_AL`, rose_gr_list$`3h_AL`)),
    '3h_FAST' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`3h_FAST`, rose_gr_list$`3h_FAST`)),
    '6h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`6h_AL`, rose_gr_list$`6h_AL`)),
    '6h_FAST' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`6h_FAST`, rose_gr_list$`6h_FAST`)),
    '12h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`12h_AL`, rose_gr_list$`12h_AL`)),
    '12h_FAST' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`12h_FAST`, rose_gr_list$`12h_FAST`)),
    '24h_AL' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`24h_AL`, rose_gr_list$`24h_AL`)),
    '24h_FAST' = as.data.frame(mergeByOverlaps(enhancer_quant_gr_list$`24h_FAST`, rose_gr_list$`24h_FAST`)))


merged_list <- lapply(merged_list, enhancer_quant_id)


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
                  values = tpm_matrix$gene,
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
detach("package:biomaRt")
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
    select(gene_id = ensembl_gene_id,
           gene_name = external_gene_name,
           chr = chromosome_name,
           tss,
           strand,
           biotype = gene_biotype)

tss_out$tss_1 <- tss_out$tss + 1 

tss_out$chr <- paste0("chr", tss_out$chr)

tss_out <- tss_out[c("chr", "tss", "tss_1", "gene_name", "biotype")]

tss_gr <- makeGRangesFromDataFrame(tss_out, seqnames.field = "chr", start.field = "tss", end.field = "tss_1",
                                   strand.field = "strand", keep.extra.columns = TRUE, ignore.strand = FALSE)


# Find genes within 50 kbp of SE by overlapping tss_gr with extended rose output ####

tss_anno_list <- list(
    '0h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`0h_AL`)),
    '1h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`1h_AL`)),
    '1h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`1h_FAST`)),
    '3h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`3h_AL`)),
    '3h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`3h_FAST`)),
    '6h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`6h_AL`)),
    '6h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`6h_FAST`)),
    '12h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`12h_AL`)),
    '12h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`12h_FAST`)),
    '24h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`24h_AL`)),
    '24h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_extended_gr$`24h_FAST`)))


tss_se_anno_list <- list(
    '0h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`0h_AL`)),
    '1h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`1h_AL`)),
    '1h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`1h_FAST`)),
    '3h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`3h_AL`)),
    '3h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`3h_FAST`)),
    '6h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`6h_AL`)),
    '6h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`6h_FAST`)),
    '12h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`12h_AL`)),
    '12h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`12h_FAST`)),
    '24h_AL' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`24h_AL`)),
    '24h_FAST' = as.data.frame(mergeByOverlaps(tss_gr, rose_se_extended_gr$`24h_FAST`)))


tss_se_anno_list <- lapply(tss_se_anno_list, function(df) {
    names(df) <- c("tss_chr", "tss_start", "tss_end", "tss_width", "tss_strand",
                   "tss_gene_name", "tss_biotype", "gene_name", "biotype",
                   "se_extended_chr", "se_extended_start", "se_extended_end",
                   "se_extended_width", "se_extended_strand", "se_extended_ID",
                   "se_extended_stitchedPeakRank", "se_extended_isSuper",
                   "ID", "stitchedPeakRank", "isSuper")
    return(df)}) %>% lapply(., function(df) {
        df <- df[c("gene_name", "se_extended_chr", "se_extended_start", "se_extended_end",
                   "se_extended_ID", "se_extended_stitchedPeakRank", "se_extended_isSuper")]
        return(df)})


# Add TPM information to TSS annotated SEs ####

samples <- names(tss_se_anno_list)


for (sample in samples) {
    
    for (i in 1:nrow(tss_se_anno_list[[sample]])) {
        
        gene <- tss_se_anno_list[[sample]]$gene_name
        
        tss_se_anno_list[[sample]]$TPM <- tpm_matrix[gene, sample]
        
    }
    
}


















