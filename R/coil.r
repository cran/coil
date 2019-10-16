

#' \pkg{coil}: evaluation of COI-5P barcode data
#'
#' \pkg{coil} is an R package designed for the cleaning, contextualization and assessment of cytochrome c oxidase I DNA
#'  barcode data (\href{https://en.wikipedia.org/wiki/Cytochrome_c_oxidase_subunit_I}{COI-5P}, or the five prime portion of COI).
#'  It contains functions for placing COI-5P barcode sequences into a common reading frame, translating DNA sequences
#'  to amino acids and for assessing the likelihood that a given barcode sequence includes an insertion or deletion error.
#'  These functions are provided as a single function analysis pipeline and are also available individually for efficient
#'  and targeted analysis of barcode data.
#'
#'@details
#'  \pkg{coil} is built around the custom `coi5p` object, which takes a COI-5P DNA barcode sequence as input.
#'  The package contains functions for: setting a sequence in reading frame, translating the sequence to amino acids
#'  and checking the sequence for evidence of insertion or deletion errors
#'
#'
#'@section Functions:
#'\itemize{
#'\item \code{\link{coi5p_pipe}} Run the entire coi5p analysis pipeline for an input sequence.
#'\item \code{\link{coi5p}} Builds a coi5p class object.
#'\item \code{\link{frame}} Sets the sequence into a common reading frame
#'\item \code{\link{which_trans_table}} Suggests a translation table for a taxonomic designation.
#'\item \code{\link{censored_translation}} Conducts translation, but ambigious codons are translated to placeholders.
#'\item \code{\link{translate}} Translate a DNA sequence to amino acids.
#'\item \code{\link{indel_check}} Check to see if an insertion or deletion error is likely.
#'}
#'
#'@section Data
#'\itemize{
#'\item \code{\link{example_nt_string}} String of DNA barcode data used in the package documentation's examples.
#'\item \code{\link{example_barcode_data}} A dataframe of coi5p barcode data, demonstrating different example cases.
#'}
#'
#'@author Cameron M. Nugent
#'
#'@docType package
#'@name coil
##################
NULL
