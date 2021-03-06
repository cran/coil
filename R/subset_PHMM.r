
#' Subset an existing PHMM.
#'
#' The \code{subsetPHMM} function allows an existing PHMM to be subset
#' by profile position. This eliminates the need for the training of an additional, smaller
#' model if query sequences should be compared to only a subsection of an existing PHMM.
#' The nt_coi_PHMM and aa_coi_PHMM can therefore be subset using this function to constrain
#' coil's framing and error evaluation to a subset of the COI-5P region
#'
#' @param x an object of class \code{"PHMM"} to be subset.
#' @param start The first PHMM position to be inclused in the output PHMM.
#' @param end The last PHMM position to be included in the output PHMM.
#' @return an object of class \code{"PHMM"}
#' @seealso \code{\link{derivePHMM}}
#'
#' @examples
#' ## subset positions 2-100 of the COI-5P PHMM
#' short_nt_PHMM <- subsetPHMM(nt_coi_PHMM, 2, 100)
#' @name subsetPHMM
subsetPHMM = function(x, start, end){

  if(end < start){
    stop("Index error, end position must match or exceed the start position.")
  }
  if(end > x$size){
    stop(paste0("Index error, end position out of bounds. Input PHMM only has a length of: ", x$size))
  }

  #get the map positions for subsetting the alignment-based fields
  map_start = x$map[[start]]
  map_end = x$map[[end]]
  #subset to get the new inserts
  new_inserts = x$inserts[map_start:map_end]

  #subset to get the new insert lengths
  new_insert_lengths = x$insertlengths[map_start:map_end]
  #reindex the insert length labels
  if(!is.null(new_insert_lengths)){
    names(new_insert_lengths) = 0:(length(new_insert_lengths)-1)
  }
  #subset the mask
  new_mask = x$mask[map_start:map_end]
  #subset the map
  new_map = x$map[start:end]
  #set the first remaining map entry equal to one,
  #make the necessary relative adjustment to all other positions
  new_map = new_map-(new_map[1]-1)
  #reindex the labels
  if(!is.null(new_map)){
    names(new_map) = 1:(length(new_map))
  }

  #Subset the A and E matrixies, reindex the positions
  #A is 0 indexed, E is 1 indexed
  new_A = x$A[,start:end]
  colnames(new_A) = 0:(length(colnames(new_A))-1)
  new_E = x$E[,start:end]
  colnames(new_E) = 1:length(colnames(new_A))

  #build the new PHMM object
  newPHMM = structure(list(name = x$name,
                           description = x$description,
                           size = (end - start)+1,
                           alphabet = x$alphabet,
                           A = new_A,
                           E = new_E,
                           qa = x$qa,
                           qe = x$qe,
                           inserts = new_inserts,
                           insertlengths = new_insert_lengths,
                           map = new_map,
                           date = x$date,
                           nseq = x$nseq,
                           weights = x$weights,
                           reference = x$reference,
                           mask = new_mask
  ), class =  "PHMM")

  return(newPHMM)
}
