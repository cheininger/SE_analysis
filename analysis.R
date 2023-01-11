# Libraries

library(tidyverse)
library(GenomicRanges)

# first read all ROSE output files into a lits

rose_out_list <- list.files(path = "./data", pattern = "*_AllStitched.table.txt", full.names = TRUE) %>% 
    setNames(., str_extract(., "([:digit:]+h\\_[:upper:]+)")) %>% 
    lapply(., read_delim, delim = "\t", col_names = TRUE, skip = 5) 



rose_gr_list <- GRangesList(lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE))

rose_gr_list_df <- lapply(rose_out_list, makeGRangesFromDataFrame, ignore.strand = FALSE)

# 
