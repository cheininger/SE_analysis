# Libraries

library(tidyverse)
library(GenomicRanges)

# Define functions

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


# first read all ROSE output files into a list

rose_out_list <- list.files(path = "./data", pattern = "*_AllStitched.table.txt", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE, skip = 5, col_types = c("cciiiidii")) 


rose_out_list <- lapply(rose_out_list, remove_chr_double)

rose_out_list <- lapply(rose_out_list, filter_chr_col)

rose_out_list <- lapply(rose_out_list, rose_enhancer_id)


# rose_gr_list <- GRangesList(lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE))

rose_gr_list <- lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)

# read enhancer quantification data into a list

enhancer_quant_list <- list.files(path = "./data", pattern = "enhancer_quant_*", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = FALSE)

enhancer_quant_list <- lapply(enhancer_quant_list, subset_quants)

enhancer_quant_list <- lapply(enhancer_quant_list, filter_chr_col)

enhancer_quant_gr_list <- lapply(enhancer_quant_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)



# pseudo-code

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
















