% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deploy_PHMMs.r
\name{set_frame}
\alias{set_frame}
\title{Take an input sequence and get it into the reading frame.}
\usage{
set_frame(org_seq, path_out)
}
\description{
Uses the path of the ntPHMM to locate the first contiguous
set of 5 matching base pairs (5 sequential 1s) for both the
front and the back of the sequence. Sequence information outside of this
first set of matches is trimmed (low probability of being true barcode sequence).
}
\keyword{internal}
