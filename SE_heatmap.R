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

# Make data frame of only isSuper values ####

issuper_df <- make_issuper_matrix(rose_unique_list)

issuper_matrix <- data.matrix(issuper_df)

issuper_df <- as.data.frame(issuper_matrix[rowSums(issuper_matrix) >= 1, ])

issuper_df <- issuper_df[, c("0h_AL_isSuper", "1h_AL_isSuper", "1h_FAST_isSuper", "6h_AL_isSuper", 
                             "6h_FAST_isSuper", "12h_AL_isSuper", "12h_FAST_isSuper", "24h_AL_isSuper", "24h_FAST_isSuper")]

names(issuper_df) <- c("0h_AL", "1h_AL", "1h_FAST", "6h_AL", "6h_FAST", "12h_AL", "12h_FAST", "24h_AL", "24h_FAST")

# Make heatmap ####

png("./plots/issuper_heatmap.png", width = 30, height = 40, units = "cm", res = 300)

issuper_heatmap <- Heatmap(as.matrix(issuper_df),
                           col = colorRamp2(c(0, 0.5, 1), c("gray", "black", "red")),
                           show_row_names = FALSE,
                           cluster_columns = FALSE,
                           column_names_gp = gpar(fontsize = 36))

draw(issuper_heatmap)

dev.off()


# Filter for SEs that appear during fasting per time point ####

coming_se_1h <- issuper_df[, c("1h_AL_isSuper", "1h_FAST_isSuper")] %>% 
    filter(., `1h_FAST_isSuper` == 1 & `1h_AL_isSuper` == 0)

coming_se_6h <- issuper_df[, c("6h_AL_isSuper", "6h_FAST_isSuper")] %>% 
    filter(., `6h_FAST_isSuper` == 1 & `6h_AL_isSuper` == 0)

coming_se_12h <- issuper_df[, c("12h_AL_isSuper", "12h_FAST_isSuper")] %>% 
    filter(., `12h_FAST_isSuper` == 1 & `12h_AL_isSuper` == 0)

coming_se_24h <- issuper_df[, c("24h_AL_isSuper", "24h_FAST_isSuper")] %>% 
    filter(., `24h_FAST_isSuper` == 1 & `24h_AL_isSuper` == 0)


# Make homerpeak file for appearing SEs per time point ####

homerpeak_1h <- rownames_to_column(coming_se_1h, var = "REGION_ID")
homerpeak_1h <- left_join(homerpeak_1h, rose_out_list$`1h_FAST`, by = "REGION_ID")
homerpeak_1h <- homerpeak_1h[c("REGION_ID", "chr", "START", "STOP")]
homerpeak_1h <- homerpeak_1h[rep(seq_len(nrow(homerpeak_1h)), each = 2), ]
homerpeak_1h <- homerpeak_1h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_1h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_1h, "./output/coming_se_1h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_6h <- rownames_to_column(coming_se_6h, var = "REGION_ID")
homerpeak_6h <- left_join(homerpeak_6h, rose_out_list$`6h_FAST`, by = "REGION_ID")
homerpeak_6h <- homerpeak_6h[c("REGION_ID", "chr", "START", "STOP")]
homerpeak_6h <- homerpeak_6h[rep(seq_len(nrow(homerpeak_6h)), each = 2), ]
homerpeak_6h <- homerpeak_6h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_6h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_6h, "./output/coming_se_6h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_12h <- rownames_to_column(coming_se_12h, var = "REGION_ID")
homerpeak_12h <- left_join(homerpeak_12h, rose_out_list$`12h_FAST`, by = "REGION_ID")
homerpeak_12h <- homerpeak_12h[c("REGION_ID", "chr", "START", "STOP")]
homerpeak_12h <- homerpeak_12h[rep(seq_len(nrow(homerpeak_12h)), each = 2), ]
homerpeak_12h <- homerpeak_12h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_12h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_12h, "./output/coming_se_12h.homerpeak", delim = "\t", col_names = FALSE)

homerpeak_24h <- rownames_to_column(coming_se_24h, var = "REGION_ID")
homerpeak_24h <- left_join(homerpeak_24h, rose_out_list$`24h_FAST`, by = "REGION_ID")
homerpeak_24h <- homerpeak_24h[c("REGION_ID", "chr", "START", "STOP")]
homerpeak_24h <- homerpeak_24h[rep(seq_len(nrow(homerpeak_24h)), each = 2), ]
homerpeak_24h <- homerpeak_24h %>% 
    mutate(strand = rep(c("+", "-"), (nrow(homerpeak_24h) / 2))) %>% 
    mutate(REGION_ID = paste0(REGION_ID, "_", strand))
write_delim(homerpeak_24h, "./output/coming_se_24h.homerpeak", delim = "\t", col_names = FALSE)







