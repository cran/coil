% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translation.r
\name{censored_translation}
\alias{censored_translation}
\title{Censored Translation of a DNA string.}
\usage{
censored_translation(dna_str, reading_frame = 1)
}
\arguments{
\item{dna_str}{The DNA string to be translated.}

\item{reading_frame}{reading frame = 1 means the first bp in the string is the start of the
first codon, can pass 1, 2 or 3. For 2 and 3 the first 1 and 2 bp will be
dropped from translation respectively.}
}
\description{
Translate a DNA sequence using the censored translation table,
this translates codons for which the amino acids is unambiguous across
mitochondrial genetic codes across the animal kingdom and does not
translate those for which the amino acid varies,
but rather outputs a ? in the string.
}
\details{
Censored translation table:
     AA  = FFLLSSSSYY?*CCWWLLLLPPPPHHQQRRRRII?MTTTTNN?KSS??VVVVAAAADDEEGGGG
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
}