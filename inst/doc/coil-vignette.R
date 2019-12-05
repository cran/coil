## ----loadlib, echo=TRUE, results='hide', message=FALSE, warning=FALSE----
library(coil)

## ------------------------------------------------------------------------
output = coi5p_pipe(example_nt_string)
output

## ------------------------------------------------------------------------
#see the available components
names(output)
#retrieve only the amino acid sequence from the object
output$aaScore

## ------------------------------------------------------------------------
ex_table_to_use = which_trans_table("Scyliorhinidae")
ex_table_to_use

## ------------------------------------------------------------------------
output = coi5p_pipe(example_nt_string, trans_table = ex_table_to_use)
output

## ------------------------------------------------------------------------
  #build the coi5p object
  dat = coi5p(example_nt_string, name="example_sequence_1")
  #frame the sequence
  dat = frame(dat)
  #since we determined the genetic code above, we can use
  #the proper transaltion table as opposed to conducting 
  #the default censored translation
  dat = translate(dat, trans_table = 2)
  #check to see if an insertion or deletion is likely
  dat = indel_check(dat)
  dat

## ------------------------------------------------------------------------
#this is the example data 
dim(example_barcode_data)
names(example_barcode_data)
# to look at the full dataframe:
# example_barcode_data

## ------------------------------------------------------------------------

example_barcode_data$coi_output = lapply(1:length(example_barcode_data$id), function(i){
  coi5p_pipe(example_barcode_data$sequence[i], 
             name = example_barcode_data$id[i], 
             trans_table = example_barcode_data$genetic_code[i])
})

example_barcode_data$coi_output[[1]] #example of the first output

## ------------------------------------------------------------------------
example_barcode_data$framed_seq = unlist(lapply(example_barcode_data$coi_output, function(x){
  x$framed
}))

#has coi5p trimmed characters?
nchar(example_barcode_data$framed_seq[[5]]) < nchar(example_barcode_data$sequence[[5]])

## ------------------------------------------------------------------------
#extract only a single column
col_df = flatten_coi5p(example_barcode_data$coi_output, keep_cols = 'aaSeq')
#extract multiple columns
multi_df = flatten_coi5p(example_barcode_data$coi_output, keep_cols = c('framed','aaSeq'))
#extract all columns
full_coi5p_df = flatten_coi5p(example_barcode_data$coi_output)
#full_coi5p_df

## ------------------------------------------------------------------------
full_coi5p_df = data.frame(matrix(ncol = 9, nrow = 0),stringsAsFactors = FALSE )
colnames(full_coi5p_df) = c("name", "raw", "framed", "was_trimmed", "align_report",
                            "aaSeq", "aaScore", "indel_likely", "stop_codons")

for(i in 1:length(example_barcode_data$id)){
	out_data = coi5p_pipe(example_barcode_data$sequence[i], 
							name = example_barcode_data$id[i], 
							trans_table = example_barcode_data$genetic_code[i])
  #for extreme memory conservation - could write each line of output to a .csv
	#instead of binding it to an output dataframe.
	full_coi5p_df = rbind(full_coi5p_df, flatten_coi5p(list(out_data)))
}

## ------------------------------------------------------------------------
library(seqinr)
# load the example fasta file included with coil
# included in the file's header line:
# the name of the sample, its genetic code, taxonomic designation and some notes
ex_fasta_file = system.file("extdata/example_barcode_data.fasta", package = "coil")

#read in the example fasta file using seqinr
ex_data = seqinr::read.fasta(ex_fasta_file, as.string = TRUE)

#here is what the output from read.fasta looks like
#head(ex_data)

#parse the data in the header line by splitting the name on the | character
parsed_names_data = lapply(1:length(ex_data), function(i){
  unlist(strsplit(names(ex_data)[[i]],"\\|"))
})

# subset the components of the header line and build these and the sequence 
# into a dataframe matching the style used in the coi5p batch example
example_barcode_data_from_scratch = data.frame(
  id = sapply(parsed_names_data, function(x) x[[1]]),
  genetic_code = sapply(parsed_names_data, function(x) x[[2]]),
  taxa = sapply(parsed_names_data, function(x) x[[3]]),
  sequence = unname(unlist(ex_data)),
  notes = sapply(parsed_names_data, function(x) x[[4]])
)

#head(example_barcode_data_from_scratch)

