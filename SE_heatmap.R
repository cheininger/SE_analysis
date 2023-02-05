# Libraries ####

library(tidyverse)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)

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

make_issuper_matrix <- function(list) {
    
    # Get names of data frames in list and add "_TPM" for later use of column names
    
    list_names <- names(list)
    list_names <- paste0(list_names, "_isSuper")
    
    # Take "REGION_ID" column as base for isSuper matrix
    
    new_df <- as.data.frame(unique(list[[1]]$REGION_ID))
    colnames(new_df) <- "REGION_ID"
    
    # Append isSuper values from every data frame to matrix
    
    for(df in list) {
        new_df <- new_df %>% inner_join(df[c("REGION_ID", "isSuper")], by = "REGION_ID")
    }
    
    # Rename isSuper columns according to data frame names
    
    colnames(new_df)[2:ncol(new_df)] <- list_names
    
    # Set "REGION_ID" column as rownames
    
    new_df <- column_to_rownames(new_df, var = "REGION_ID")
    
    # Output resulting data frame
    
    return(new_df)
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


# Make matrix of super-enhancers for heatmap ####

# Filter every data frame to not contain duplicate IDs

rose_unique_list <- lapply(rose_out_list, function(df) {
    
    df <- df[!base::duplicated(df$REGION_ID),]
    
    return(df)
    
})

# Make data frame of only isSuper values

issuper_df <- make_issuper_matrix(rose_unique_list)

issuper_matrix <- data.matrix(issuper_df)

issuper_df <- as.data.frame(issuper_matrix[rowSums(issuper_matrix) >= 1, ])

issuper_df <- issuper_df[, c("0h_AL_isSuper", "1h_AL_isSuper", "1h_FAST_isSuper", "3h_AL_isSuper", "3h_FAST_isSuper",
                             "6h_AL_isSuper", "6h_FAST_isSuper", "12h_AL_isSuper", "12h_FAST_isSuper", "24h_AL_isSuper", "24h_FAST_isSuper")]

# Make heatmap

issuper_heatmap <- Heatmap(as.matrix(issuper_df),
                           col = colorRamp2(c(0, 0.5, 1), c("gray", "black", "red")),
                           show_row_names = FALSE,
                           cluster_columns = FALSE)

draw(issuper_heatmap)






