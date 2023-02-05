# Libraries ####

library(tidyverse)
library(GenomicRanges)


# Define functions ####

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


rose_out_list <- lapply(rose_out_list, function(df) {
    df <- df %>% mutate("width" = STOP - START)
    
    return(df)
})


rose_super_list <- lapply(rose_out_list, function(df) {
    df <- filter(df, isSuper == 1) %>% 
        mutate("width" = STOP - START)
    return(df)
})


# Make density plot for width of ROSE super-enhancer regions ####


list_names <- names(rose_super_list)

for (i in 1:length(rose_super_list)) {
    
    name <- list_names[[i]]
    
    df <- rose_super_list[[i]]
    
    plot <- ggplot() + 
        geom_density(data = df, 
                     mapping = aes(x = df$width)) +
        theme_classic() 
    
    ggsave(filename = paste0("./plots/", name, "_width_dens.png"), plot = plot, device = "png", dpi = 300)
    
    rm(plot)
    rm(i)
    rm(name)
    
}


# Make data frame with ROSE region IDs as rownames and samples as columns ####

issuper_df <- make_issuper_matrix(rose_out_list)




