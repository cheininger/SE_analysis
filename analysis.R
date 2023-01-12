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


# first read all ROSE output files into a list

rose_out_list <- list.files(path = "./data", pattern = "*_AllStitched.table.txt", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE, skip = 5) 


rose_out_list <- lapply(rose_out_list, remove_chr_double)

rose_out_list <- lapply(rose_out_list, filter_chr_col)


# rose_gr_list <- GRangesList(lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE))

rose_gr_list_df <- lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE, keep.extra.columns = TRUE)

# read enhancer quantification data into a list

enhancer_quant_list <- list.files(path = "./data", pattern = "enhancer_quant_*", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = FALSE)

enhancer_quant_list <- lapply(enhancer_quant_list, subset_quants)

enhancer_quant_list <- lapply(enhancer_quant_list, filter_chr_col)

enhancer_quant_gr_list <- lapply(enhancer_quant_list, makeGRangesFromDataFrame, ignore.strand = TRUE, keep.extra.columns = FALSE)



# pseudo-code

mergeByOverlaps(enhancer_quant_list_gr[[1]], rose_gr_list_df[[1]])

















