###setwd('/Users/gaga/Desktop/coverage_cutoff/cov_aoto_test') # will be removed

# get sample name from command line
Args <- commandArgs()

# read output from extractID.py into R
library(dplyr)
###input <- paste0('extracted_', 'ERR011103', '.txt') # will be removed
input <- paste0('extracted_', Args[6], '.txt')
ori <- read.table(input, header = FALSE)
colnames(ori) <- c('taxid', 'm_length', 'r_length')

# read genome summary into R
##idxstat <- paste0(Args[6], '.all.txt') # bacteria, virus and eukaryota
idxstat <- paste0('/public/users/liangminling/liangqx/genome_summary/All/', Args[6], '.all.txt') # bacteria only, for testing cutoff
###idxstat <- paste0('ERR011103', '.all.txt') # will be remove
gs <- read.table(idxstat, header = FALSE, fill = FALSE ,quote = '', skip = 1, sep = '\t')
colnames(gs) <- c('taxid', 'tax_name', 'avg_genome_length', 'reads_mapped', 'length_normalized_reads')

# column 4 - column3, add up the new column of the same taxid
then <- ori %>%
        mutate(real_map = r_length - m_length) %>%
        group_by(taxid) %>%
        summarize(map_length = sum(m_length), ref_length = sum(r_length), real_map_length = sum(real_map)) 

# join the table generated above with genome summary
big <- left_join(gs, then, by = 'taxid') %>%
        mutate(true_genome_length = rep(NA, n())) %>%
        filter(is.na(ref_length) == FALSE)

# determine genome length by the greater one between ref_length and avg_genome_length
for(i in 1:nrow(big)){
        if(big$ref_length[i] > big$avg_genome_length[i]){
                big$true_genome_length[i] <- big$ref_length[i]
        }
        if(big$ref_length[i] <= big$avg_genome_length[i]){
                big$true_genome_length[i] <- big$avg_genome_length[i]
        }
}

# calculate coverage by true_genome_length
cover <- big %>%
        mutate(real_cov = real_map_length / true_genome_length)
       

# output 7 table based on 7 kinds of cutoff
cutting <- function(data, cutoff){ # cutoff is a fraction：1X should be written as 0.01
        data %>%
                filter(real_cov > cutoff) %>%
                select(taxid, tax_name, avg_genome_length, reads_mapped, length_normalized_reads) %>%
                ###write.table(file = paste0('cut', as.character(cutoff), '_', 'ERR011103', '.txt')) # will be removed
                write.table(file = paste0('cut', as.character(cutoff * 100), '_', Args[6], '.txt'))
}

##source('/public/users/liangminling/liangqx/R_script/replace_genome_length_cov.R')

cutting(cover, 0) # no cutoff
cutting(cover, 0.0001) # 0.01% 
cutting(cover, 0.001) # 0.1%
cutting(cover, 0.01) # 1%
cutting(cover, 0.05) # 5%
cutting(cover, 0.1) # 10%
cutting(cover, 0.2) # 20%