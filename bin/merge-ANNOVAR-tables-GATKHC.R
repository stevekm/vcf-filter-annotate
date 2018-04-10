#!/usr/bin/env Rscript

# ~~~~~~ CUSTOM FUNCTIONS ~~~~~~~ #
read.ANNOVAR.vcf_txt <- function(file, na_string = '.'){
    # read the '_multianno.txt' file produced by ANNOVAR table_annovar.pl using the '--vcfinput' arg
    # tab-delimited text file has more columns than column names
    # output dataframe with colnames present
    
    # read in the file
    tmp_annovar_df <- read.delim(file = file, 
                                 header = FALSE, 
                                 sep = '\t', 
                                 check.names = FALSE, 
                                 fill = TRUE, stringsAsFactors = FALSE, 
                                 na.strings = na_string)
    # get the colnames that are present
    av_colnames <- tmp_annovar_df[1,][which(! is.na(tmp_annovar_df[1,]) & ! tmp_annovar_df[1,] == "" )]
    # create new df without colnames
    tmp_annovar_df2 <- tmp_annovar_df[2:nrow(tmp_annovar_df),]
    # add the colnames
    names(tmp_annovar_df2)[1:length(av_colnames)] <- av_colnames
    return(tmp_annovar_df2)
}

all_equal <- function(...){
    # returns TRUE or FALSE if all the elements passed are equal
    arguments <- list(...)
    TF <- vapply(1:(length(arguments)-1),
                 function(n) identical(arguments[[n]], arguments[[n+1]]),
                 logical(1))
    return(all(TF))
}


# ~~~~~~~ RUN ~~~~~~~ #
# merge ANNOVAR output with the .vcf recalculated .tsv file
args <- commandArgs(TRUE)
# sampleID <- "NC-HAPMAP"
# recalc_tsv_file <- "output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.recalc.tsv"
# annovar_file <- "output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.hg19_multianno.txt"
# avinput_file <- "output/NC-HAPMAP/HaplotypeCaller/NC-HAPMAP.avinput"
# output_file <- "NC-HAPMAP.annotation.tsv"
sampleID <- args[1]
recalc_tsv_file <- args[2]
annovar_file <- args[3]
avinput_file <- args[4]
output_file <- args[5]

recalc_df <- read.delim(file = recalc_tsv_file, header = TRUE, sep = '\t', check.names = FALSE)
avinput_df <- read.delim(file = avinput_file, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
annovar_df <- read.ANNOVAR.vcf_txt(file = annovar_file)


# check that all df's have the same number of rows
if(! all_equal(nrow(annovar_df), nrow(avinput_df), nrow(recalc_df)) ) {
    print("ERROR: dataframes have unequal numbers of rows")
    quit(status = 1)
}


avinput_colnames <-  c("Chr", "Start", "End", "Ref", "Alt", "AF", "QUAL", 
                       "AD", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                       "FILTER", "INFO", "FORMAT", sampleID)
names(avinput_df)[1:length(avinput_colnames)] <- avinput_colnames


# merge the recalc table against avinput to get the FREQ col
avinput_recalc <- merge(x = avinput_df, y = recalc_df[c("CHROM", "POS", "REF", "ALT", "FREQ")], by = c("CHROM", "POS", "REF", "ALT"))

# merge the FREQ col from avinput to annovar output
annovar_merged <- merge(x = avinput_recalc[c("Chr", "Start", "End", "Ref", "Alt", "FREQ")], y = annovar_df, by = c("Chr", "Start", "Ref", "Alt", "End"))

# remove extraneous columns
keep_cols <- names(annovar_merged)[which(! grepl(pattern = "^V[[:digit:]]*$", x = names(annovar_merged)) & ! grepl(pattern = "^Otherinfo$", x = names(annovar_merged)) )]
annovar_merged <- annovar_merged[keep_cols]

# check that all df's have the same number of rows
if(! all_equal(nrow(annovar_df), 
               nrow(avinput_df), 
               nrow(recalc_df), 
               nrow(avinput_recalc),
               nrow(annovar_merged)) ) {
    print("ERROR: dataframes have unequal numbers of rows")
    quit(status = 1)
}

write.table(x = annovar_merged, file = output_file, sep = '\t', quote = FALSE, na = '.', row.names = FALSE, col.names = TRUE)
