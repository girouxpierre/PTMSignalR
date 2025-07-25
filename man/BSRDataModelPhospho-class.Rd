% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSRDataModel-class.R
\docType{class}
\name{BSRDataModelPhospho-class}
\alias{BSRDataModelPhospho-class}
\title{BulkSignalR Data Model Object}
\description{
  An S4 class to represent the expression data used for inferring
  ligand-receptor interactions.
}
\section{Slots}{
  
  \describe{
    \item{\code{ncounts}}{Normalized read count matrix. Row names must be set to HUGO
      official gene symbols.}
    
    \item{\code{log.transformed}}{Logical indicating whether values in
      \code{ncounts} were log2-transformed.}
    
    \item{\code{normalization}}{Name of the normalization method.}
    
    \item{\code{param}}{List containing the statistical model parameters.}
    
    \item{\code{initial.organism}}{Organism for which the data were obtained.}
    
    \item{\code{initial.orthologs}}{List of genes for which human
      orthologs exist.}
  }}

\examples{
  new("BSRDataModelPhospho", ncounts=matrix(1.5, nrow=2, ncol=2,
                                     dimnames=list(c("A","B"), c("C","D"))),
      log.transformed=TRUE,
      normalization="TC")
  
}
