library(methods)

#' BulkSignalR Data Model Compare PTM Object
#'
#' An S4 class to represent the expression data used for inferring
#' ligand-receptor interactions based on sample cluster comparisons.
#'
#' @slot comp   A named list of BSRClusterComp objects, one per
#' comparison.
#' @slot compPTM   A named list of BSRClusterComp objects, one per
#' comparison.
#' @export
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp, max.pval=1,"random.example")
#' 
setClass("BSRDataModelCompPTM",
         contains=c("BSRDataModelPTM"),
         slots=c(comp="list",
                 compPTM="list"),
         prototype=list(
           initial.organism="hsapiens",
           initial.orthologs=list("A","B","C"),
           ncounts=matrix(1.0,nrow=2,ncol=1,dimnames=list(c("A","B"),"C")),
           PTM=matrix(1.0,nrow=2,ncol=1,dimnames=list(c("A","B"),"C")),
           symPos=matrix(1.0,nrow=2,ncol=1,dimnames=list(c("A","B"),"C")),
           log.transformed=FALSE,
           single=FALSE,
           normalization="UQ",
           param=list(spatial.smooth=FALSE),
           comp=list(),
           compPTM=list()
         ))

setValidity("BSRDataModelCompPTM",
            function(object) {
              if (length(object@comp) > 0){
                if (is.null(names(object@comp)))
                  return("comp names must be set")
                for (c in names(object@comp))
                  if (!is(object@comp[[c]], "BSRClusterComp"))
                    return("comp must contain objects of class BSRClusterComp")
              }
              
              TRUE
            }
)

setMethod("show", "BSRDataModelCompPTM",
          function(object) {
            callNextMethod()
            cat("Defined comparisons:\n")
            utils::str(object@compPTM)
          }
)


# Accessors & setters ========================================================

if (!isGeneric("comp")) {
  if (is.function("comp"))
    fun <- comp
  else
    fun <- function(x) standardGeneric("comp")
  setGeneric("comp", fun)
}
#' Comparisons list accessor
#'
#' @name comp
#' @aliases comp,BSRDataModelCompPTM-method
#' @param x object BSRDataModelCompPTM 
#' @export
setMethod("comp", "BSRDataModelCompPTM", function(x) x@comp)

if (!isGeneric("comp<-")) {
  if (is.function("comp<-"))
    fun <- `comp<-`
  else
    fun <- function(x, value) standardGeneric("comp<-")
  setGeneric("comp<-", fun)
}
#' Comparisons list setter (internal use only, use addComparison() otherwise)
#'
#' @param x object BSRDataModelCompPTM 
#' @param value value to be set for BSRDataModelCompPTM
#' @keywords internal 
setMethod("comp<-", "BSRDataModelCompPTM", function(x,value){
  x@comp <- value
  methods::validObject(x)
  x
})


if (!isGeneric("compPTM")) {
  if (is.function("compPTM"))
    fun <- compPTM
  else
    fun <- function(x) standardGeneric("compPTM")
  setGeneric("compPTM", fun)
}
#' Comparisons list accessor
#'
#' @name compPTM
#' @aliases compPTM,BSRDataModelCompPTM-method
#' @param x object BSRDataModelCompPTM 
#' @export
setMethod("compPTM", "BSRDataModelCompPTM", function(x) x@compPTM)

if (!isGeneric("compPTM<-")) {
  if (is.function("compPTM<-"))
    fun <- `compPTM<-`
  else
    fun <- function(x, value) standardGeneric("compPTM<-")
  setGeneric("compPTM<-", fun)
}
#' Comparisons list setter (internal use only, use addComparison() otherwise)
#'
#' @param x object BSRDataModelCompPTM 
#' @param value value to be set for BSRDataModelCompPTM
#' @keywords internal 
setMethod("compPTM<-", "BSRDataModelCompPTM", function(x,value){
  x@compPTM <- value
  methods::validObject(x)
  x
})


# generation of BSRDataModelCompPTM from BSRDataModelPTM =====================

#' @title conversion of BSRDataModelPTM into BSRDataModelCompPTM 
#'
#' @description In case ligand-receptor inferences should be obtained
#' based on gene/protein regulation P-values comparing two clusters of
#' samples, it is necessary to first promote the BSRDataModelPTM object that
#' contains the count matrix into a BSRDataModelCompPTM object able to contain
#' a list of such cluster pairs comparisons. This function performs this
#' promotion adding an empty list of comparisons.
#' @param bsrdm    A BSRDataModelPTM object.
#'
#' @export
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' 
as.BSRDataModelCompPTM <- function(bsrdm){
  
  if (!is(bsrdm, "BSRDataModelPTM"))
    stop("bsrdm must be of class BSRDataModelPTM")
  new("BSRDataModelCompPTM", bsrdm, comp=list())
  
} # as.BSRDataModelCompPTM


# defining/adding a cluster comparison ========================================

if (!isGeneric("defineClusterComp")) {
  if (is.function("defineClusterComp"))
    fun <- defineClusterComp
  else
    fun <- function(obj, ...) standardGeneric("defineClusterComp")
  setGeneric("defineClusterComp", fun)
}
#' Definition of the comparison between two clusters of samples
#'
#' Define the columns of the expression matrix that belong to each cluster,
#' and store the result of the cluster differences statistical analysis
#' obtained by an external tool such as edgeR, DESeq2, etc.
#'
#' @name defineClusterComp
#' @aliases defineClusterComp,BSRDataModelCompPTM-method
#'
#' @param obj    A BSRDataModelCompPTM object output by
#'   \code{\link{as.BSRDataModelCompPTM}}.
#' @param colA   Cluster A column indices.
#' @param colB   Cluster B column indices.
#' @param stats  A data.frame containing statistics about the differential
#' analysis cluster A versus B. \code{stats} must contain at least the
#' columns 'pval' (for P-values) and 'logFC' for log-fold-changes A/B.
#' @param statsPTM  A data.frame containing statistics about the differential
#' analysis cluster A versus B. \code{stats} must contain at least the
#' columns 'pval' (for P-values) and 'logFC' for log-fold-changes A/B.
#'
#' @details Create a BSRClusterComp object describing a comparison
#' of two clusters of columns taken from the expression matrix
#' in the BSRDataModelCompPTM object \code{obj}. Such a cluster comparison
#' description is the basis for inferring LRIs from differential
#' expression P-values instead of correlation analysis.
#' 
#' The rows of \code{stats} must be in the same order as those of the count
#' matrix in \code{obj}. Alternatively, \code{stats}
#' rows can be named and a 1-1 correspondence must exist between these names
#' and those of the count matrix.
#'
#' @return A BSRClusterCompPTM object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp, max.pval=1,"random.example")
#' 
#' @importFrom methods new
setMethod("defineClusterComp", "BSRDataModelCompPTM", function(obj, colA,
                                                                   colB, stats, statsPTM){
  
  if (!is.integer(colA))
    stop("colA must contain integer indices")
  if (!is.integer(colB))
    stop("colB must contain integer indices")
  if (length(intersect(colA, colB)) > 0)
    stop("colA and colB must be disjoint")
  if (any(colA < 1 | colA > ncol(ncounts(obj))))
    stop("colA indices must fall in [1; ncol(ncounts)]")
  if (any(colB < 1 | colB > ncol(ncounts(obj))))
    stop("colB indices must fall in [1; ncol(ncounts)]")
  if (!is.data.frame(stats))
    stop("stats must be a data.frame")
  if (!all(c("pval","logFC") %in% names(stats)))
    stop("stats data.frame must contain columns named 'pval' and 'logFC'")
  if (nrow(stats) != nrow(ncounts(obj)))
    stop("stats and ncounts(obj) number of rows differ")
  if (!is.null(rownames(stats)) &&
      (sum(rownames(stats) %in% rownames(ncounts(obj))) != nrow(stats)))
    stop("stats rownames defined but do not all match ncounts(obs)")
  if (!is.null(rownames(statsPTM)) &&
      (sum(rownames(statsPTM) %in% rownames(PTM(obj))) != nrow(statsPTM)))
    stop("statsPTM rownames defined but do not all match ncounts(obs)")
  if (is.null(rownames(stats)))
    rownames(stats) <- rownames(ncounts(obj))
  
  new("BSRClusterCompPTM", colA=colA, colB=colB, stats=stats, statsPTM=statsPTM)
  
}) # defineClusterCompPTM


if (!isGeneric("addClusterComp")) {
  if (is.function("addClusterComp"))
    fun <- addClusterComp
  else
    fun <- function(obj, ...) standardGeneric("addClusterComp")
  setGeneric("addClusterComp", fun)
}
#' Add a comparison between two clusters of samples
#'
#' Add a comparison to a BSRDataModelCompPTM object.
#'
#' @name addClusterComp
#' @aliases addClusterComp,BSRDataModelCompPTM-method
#'
#' @param obj    A BSRDataModelCompPTM object output by
#'   \code{\link{as.BSRDataModelCompPTM}}.
#' @param cmp   A BSRClusterComp object to add.
#' @param cmp.name  The name of the comparison to add.
#'
#' @details Add \code{cmp} to the list of comparisons contained in
#' \code{obj}.
#' 
#' @return A BSRDataModelCompPTM object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' @importFrom methods new
setMethod("addClusterComp", "BSRDataModelCompPTM", function(obj, cmp,
                                                                cmp.name){
  
  if (!is(cmp, "BSRClusterCompPTM"))
    stop("cmp must be of class BSRClusterCompPTM")
  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (cmp.name %in% names(compPTM(obj)))
    stop("cmp.name is already in the list of comparisons")
  
  tmp <- c(compPTM(obj), list(cmp))
  names(tmp)[length(tmp)] <- cmp.name
  compPTM(obj) <- tmp
  obj
  
}) # addClusterComp


if (!isGeneric("addClusterCompPTM")) {
  if (is.function("addClusterCompPTM"))
    fun <- addClusterCompPTM
  else
    fun <- function(obj, ...) standardGeneric("addClusterCompPTM")
  setGeneric("addClusterCompPTM", fun)
}
#' Add a comparison between two clusters of samples
#'
#' Add a comparison to a BSRDataModelCompPTM object.
#'
#' @name addClusterCompPTM
#' @aliases addClusterCompPTM,BSRDataModelCompPTM-method
#'
#' @param obj    A BSRDataModelCompPTM object output by
#'   \code{\link{as.BSRDataModelCompPTM}}.
#' @param cmp   A BSRClusterComp object to add.
#' @param cmp.name  The name of the comparison to add.
#'
#' @details Add \code{cmp} to the list of comparisons contained in
#' \code{obj}.
#' 
#' @return A BSRDataModelCompPTM object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' @importFrom methods new
setMethod("addClusterCompPTM", "BSRDataModelCompPTM", function(obj, cmp,
                                                                    cmp.name){
  
  if (!is(cmp, "BSRClusterCompPTM"))
    stop("cmp must be of class BSRClusterCompPTM")
  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (cmp.name %in% names(compPTM(obj)))
    stop("cmp.name is already in the list of comparisons")
  
  tmp <- c(compPTM(obj), list(cmp))
  names(tmp)[length(tmp)] <- cmp.name
  compPTM(obj) <- tmp
  obj
  
}) # addClusterCompPTM


if (!isGeneric("removeClusterComp")) {
  if (is.function("removeClusterComp"))
    fun <- removeClusterComp
  else
    fun <- function(obj, ...) standardGeneric("removeClusterComp")
  setGeneric("removeClusterComp", fun)
}
#' Remove a comparison from a BSRDataModelCompPTM object.
#'
#' @name removeClusterComp
#' @aliases removeClusterComp,BSRDataModelCompPTM-method
#'
#' @param obj    A BSRDataModelCompPTM object output by
#'   \code{\link{as.BSRDataModelCompPTM}}.
#' @param cmp.name  The name of the comparison to remove.
#'
#' @details Remove the comparison with \code{cmp.name} from the list of
#' comparisons contained in \code{obj}.
#' 
#' @return A BSRDataModelCompPTM object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' 
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' bsrdm.comp
#' bsrdm.comp <- removeClustercompPTM(bsrdm.comp, "random.example")
#' bsrdm.comp
#' 
setMethod("removeClusterComp", "BSRDataModelCompPTM", function(obj, cmp.name){
  
  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (!(cmp.name %in% names(compPTM(obj))))
    stop("cmp.name must be in the list of comparisons")
  
  if (length(compPTM(obj)) == 1)
    compPTM(obj) <- list()
  else
    compPTM(obj) <- compPTM(obj)[names(compPTM(obj)) != cmp.name]
  
  obj
  
}) # removeClusterComp


if (!isGeneric("removeClusterCompPTM")) {
  if (is.function("removeClusterCompPTM"))
    fun <- removeClusterCompPTM
  else
    fun <- function(obj, ...) standardGeneric("removeClusterCompPTM")
  setGeneric("removeClusterCompPTM", fun)
}
#' Remove a comparison from a BSRDataModelCompPTM object.
#'
#' @name removeClusterCompPTM
#' @aliases removeClusterCompPTM,BSRDataModelCompPTM-method
#'
#' @param obj    A BSRDataModelCompPTM object output by
#'   \code{\link{as.BSRDataModelCompPTM}}.
#' @param cmp.name  The name of the comparison to remove.
#'
#' @details Remove the comparison with \code{cmp.name} from the list of
#' comparisons contained in \code{obj}.
#' 
#' @return A BSRDataModelCompPTM object.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' 
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' bsrdm.comp
#' bsrdm.comp <- removeClustercompPTM(bsrdm.comp, "random.example")
#' bsrdm.comp
#' 
setMethod("removeClusterCompPTM", "BSRDataModelCompPTM", function(obj, cmp.name){
  
  if (!is.character(cmp.name))
    stop("cmp.name must be of type character")
  if (length(cmp.name) == 0)
    stop("cmp.name must have length > 0")
  if (!(cmp.name %in% names(compPTM(obj))))
    stop("cmp.name must be in the list of comparisons")
  
  if (length(compPTM(obj)) == 1)
    compPTM(obj) <- list()
  else
    compPTM(obj) <- compPTM(obj)[names(compPTM(obj)) != cmp.name]
  
  obj
  
}) # removeClusterCompPTM


# performing initial inference ===========================================


if (!isGeneric("initialInference")) {
  if (is.function("initialInference"))
    fun <- initialInference
  else
    fun <- function(obj, ...) standardGeneric("initialInference")
  setGeneric("initialInference", fun)
}
#' Inference of ligand-receptor interactions based on regulation
#'
#' Computes putative LR interactions along with their statistical confidence.
#' The analysis is based on P-values associated to the regulation of
#' genes/proteins in a comparison of two clusters of samples.
#' In this initial inference, all the relevant pathways are reported,
#' see reduction functions to reduce this list.
#'
#' @name initialInference
#' @aliases initialInference,BSRDataModelCompPTM-method
#'
#' @param obj       A BSRDataModelCompPTM object.
#' @param cmp.name        The name of the cluster comparison that should be used for
#'   the inference.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param max.pval        The maximum P-value imposed to both the ligand
#'   and the receptor.
#' @param min.logFC       The minimum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
#' @param fdr.proc      The procedure for adjusting P-values according to
#' \code{\link[multtest]{mt.rawp2adjp}}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#'
#' @details Perform the initial ligand-receptor inference. Initial means that
#' no reduction is applied. All the (ligand, receptor, downstream pathway)
#' triples are reported, i.e., a given LR pair may appear multiple times
#' with different pathways downstream the receptor. Specific reduction
#' functions are available from the package to operate subsequent
#' simplifications based on the BSRInferenceComp object created by this method.
#' 
#' Here, ligand-receptor interactions are inferred based on gene or protein
#' regulation-associated P-values when comparing two clusters of samples. Since
#' a BSRDataModelCompPTM object can contain several such comparisons, the name
#' of the comparison to use must be specified (parameter \code{cmp.name}).
#'
#' @return A BSRInferenceComp object with initial inferences set.
#'
#' @export
#'
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' @importFrom methods new
setMethod("initialInference", "BSRDataModelCompPTM", function(obj, cmp.name, rank.p=0.55,
                                                                  max.pval=0.01, min.logFC=1, neg.receptors=FALSE,
                                                                  restrict.genes=NULL, reference=c("REACTOME-GOBP","REACTOME","GOBP"),
                                                                  max.pw.size=200, min.pw.size=5, min.positive=4, restrict.pw=NULL,
                                                                  with.complex=TRUE, fdr.proc=c("BH","Bonferroni","Holm","Hochberg",
                                                                                                "SidakSS","SidakSD","BY","ABH","TSBH")){
  
  if (!(cmp.name %in% names(compPTM(obj))))
    stop("cmp.name must exist in the names of comparisons contained in obj")
  reference <- match.arg(reference)
  fdr.proc <- match.arg(fdr.proc)
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  
  # retrieve the BSRClusterComp object
  cc <- compPTM(obj)[[cmp.name]]
  print(cc)
  
  inf.param <- list()
  inf.param$colA <- colA(cc)
  inf.param$colB <- colB(cc)
  inf.param$max.pval <- max.pval
  inf.param$min.logFC <- min.logFC
  inf.param$neg.receptors <- neg.receptors
  inf.param$restrict.genes <- restrict.genes
  lr <- .getRegulatedLRPTM(obj, cc, max.pval=max.pval, min.logFC=min.logFC,
                        neg.receptors=neg.receptors, restrict.genes=restrict.genes)
  
  inf.param$reference <- reference
  inf.param$min.pw.size <- min.pw.size
  inf.param$max.pw.size <- max.pw.size
  inf.param$with.complex <- with.complex
  inf.param$min.positive <- min.positive
  inf.param$restrict.pw <- restrict.pw
  pairs <- .checkRegulatedReceptorSignalingPTM(obj, cc, lr, reference=reference,
                                            min.pw.size=min.pw.size, max.pw.size=max.pw.size,
                                            min.positive=min.positive, with.complex=with.complex,
                                            restrict.pw=restrict.pw)
  
  
  inf.param$fdr.proc <- fdr.proc
  inf.param$rank.p <- rank.p
  #inter <- .pValuesRegulatedLRPTM(pairs, param(obj), rank.p=rank.p, fdr.proc=fdr.proc)
  inter <- .pValuesRegulatedLRPTM(pairs, param(obj), rank.p=rank.p, fdr.proc=fdr.proc)
  cat("\n")
  cat(unlist(colnames(inter)))
  
  ligands <- strsplit(inter$L, ";")
  receptors <- strsplit(inter$R, ";")
  tg <- strsplit(inter$target.genes, ";")
  PTMg <- strsplit(inter$PTM.genes, ";")
  pg <- strsplit(inter$PTM.genes, ";")
  dpg <- strsplit(inter$dePTM.genes, ";")
  tgpval <- lapply(strsplit(inter$target.pval, ";"),
                   function(x) as.numeric(x))
  
  # cat("\n protLFC \n")
  # cat(unlist(inter$target.logFC))
  tglogFC <- lapply(strsplit(inter$target.logFC, ";"),
                    function(x) as.numeric(x))
  tgcorr <- lapply(strsplit(inter$target.corr, ";"),
                   function(x) as.numeric(x))
  
  PTMglogFC <- NA
  # cat("\n PTMLFC \n")
  # cat(unlist(inter$PTM.logFC))
  PTMglogFC <- lapply(strsplit(inter$PTM.logFC[!is.na(inter$PTM.logFC)], ";"),
                    function(x) as.numeric(x))
  PTMgcorr <- NA
  PTMgcorr <- lapply(strsplit(inter$PTM.corr[!is.na(inter$PTM.corr)], ";"),
                   function(x) as.numeric(x))
  
  pglogFC <- NA
  # cat("\n PTMLFC \n")
  # cat(unlist(inter$PTM.logFC))
  # pglogFC <- lapply(strsplit(inter$PTM.logFC[!is.na(inter$PTM.logFC)], ";"),
  #                   function(x) as.numeric(x))
  pglogFC <- lapply(strsplit(as.character(inter$PTM.logFC), ";"),
                    function(x) as.numeric(x))
  pgcorr <- NA
  # pgcorr <- lapply(strsplit(inter$PTM.corr[!is.na(inter$PTM.corr)], ";"),
  #                  function(x) as.numeric(x))
  pgcorr <- lapply(strsplit(as.character(inter$PTM.corr), ";"),
                   function(x) as.numeric(x))
  
  dpglogFC <- NA
  # cat("\n PTMLFC \n")
  # cat(unlist(inter$PTM.logFC))
  # dpglogFC <- lapply(strsplit(inter$dePTM.logFC[!is.na(inter$dePTM.logFC)], ";"),
  #                   function(x) as.numeric(x))
  dpglogFC <- lapply(strsplit(as.character(inter$dePTM.logFC), ";"),
                     function(x) as.numeric(x))
  dpgcorr <- NA
  # dpgcorr <- lapply(strsplit(inter$dePTM.corr[!is.na(inter$dePTM.corr)], ";"),
  #                  function(x) as.numeric(x))
  dpgcorr <- lapply(strsplit(as.character(inter$dePTM.corr), ";"),
                    function(x) as.numeric(x))
  inf.param$ligand.reduced <- FALSE
  inf.param$receptor.reduced <- FALSE
  inf.param$pathway.reduced <- FALSE
  
  new("BSRInferenceCompPTM", LRinter=inter[,c("L","R","pw.id","pw.name","pval","qval","L.logFC","R.logFC","LR.corr","LR.pval","pvalLRT","pvalRPTM",
                                                  "rank","rank.corr","len","lenPTM", "lenp", "lendp", "PTM.genes","PTM.logFC","PTM.genes","PTM.logFC","dePTM.genes","dePTM.logFC","target.genes")], ligands=ligands,
      receptors=receptors, t.genes=tg, tg.corr=tgcorr, PTM.genes=PTMg, PTMg.corr=PTMgcorr, p.genes=pg, pg.corr=pgcorr, dp.genes=dpg, dpg.corr=dpgcorr,
      tg.logFC=tglogFC, pg.logFC=pglogFC, inf.param=inf.param, cmp.name=cmp.name)
  
  # new("BSRInferenceCompPTM", LRinter=inter[,c("L","R","pw.id","pw.name","pval","qval","L.logFC","R.logFC","LR.pval","LR.corr",
  #                                        "rank","len","rank.pval","rank.corr")],
  #     ligands=ligands, receptors=receptors, t.genes=tg, tg.corr=tgcorr,
  #     tg.pval=tgpval, tg.logFC=tglogFC, inf.param=inf.param, cmp.name=cmp.name)
  
}) # initialInference


# Scoring of gene signatures in a BSRSignature object ==========================

if (!isGeneric("scoreLRGeneSignatures")) {
  if (is.function("scoreLRGeneSignatures"))
    fun <- scoreLRGeneSignatures
  else
    fun <- function(obj, ...) standardGeneric("scoreLRGeneSignatures")
  setGeneric("scoreLRGeneSignatures", fun)
}
#' Score ligand-receptor gene signatures
#'
#' Compute ligand-receptor gene signature scores over a BSRDataModelCompPTM
#' specific comparison.
#'
#' @name scoreLRGeneSignatures
#' @aliases scoreLRGeneSignatures,BSRDataModelCompPTM-method
#'
#' @param obj           A BSRDataModelCompPTM object.
#' @param sig           A BSRSignatureComp object.
#' @param LR.weight    A number between 0 and 1 defining the relative weight
#' of the ligand and the receptor in the signature.
#' @param robust       A logical indicating that z-scores should be computed
#' with median and MAD instead of mean and standard deviation.
#' @param name.by.pathway     A logical indicating whether row names of the
#' resulting score matrix should be pathway names.
#' @param rownames.LRP A logical indicating, in case \code{name.by.pathway}
#' was set to TRUE, whether ligand and receptor names should be added on top.
#' No role if \code{name.by.pathway} was set to FALSE.
#' @param abs.z.score  A logical to use absolute z-scores (useful if the
#' activity of a paythway is reported by a mixture of up- and down-genes
#' whose z-score averages might hide actual activity).
#' @return A matrix containing the scores of each ligand-receptor gene
#' signature in each sample.
#'
#' @export
#' @examples
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelCompPTM(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClustercompPTM(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClustercompPTM(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1,  "random.example")
#' 
#' # reduction
#' bsrinf.red <- reduceToBestPathway(bsrinf)
#' 
#' # signature extraction and scoring
#' bsrsig.red <- getLRGeneSignatures(bsrinf.red, qval.thres=1e-6)
#' scores.red <- scoreLRGeneSignatures(bsrdm.comp, bsrsig.red,
#'                          name.by.pathway=TRUE, rownames.LRP=TRUE)
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom methods is
setMethod("scoreLRGeneSignatures", "BSRDataModelCompPTM", function(obj,
                                                                       sig, LR.weight=0.5, robust=FALSE,
                                                                       name.by.pathway=FALSE, abs.z.score=FALSE,rownames.LRP=FALSE){
  
  if (!is(sig, "BSRSignatureComp"))
    stop("sig must be a BSRSignature object")
  if (LR.weight<=0 || LR.weight>=1)
    stop("LRweight must reside in (0;1)")
  
  # retrieve the BSRClusterComp object
  cmp.name <- cmpName(sig)
  if (!(cmp.name %in% names(compPTM(obj))))
    stop("The comparison name in sig is not present in obj")
  cc <- compPTM(obj)[[cmp.name]]
  
  # species management  
  if( initialOrganism(obj)!="hsapiens" )
    all.genes <- unlist(initialOrthologs(obj))
  else
    all.genes <- rownames(ncounts(obj))
  
  # get the ncount matrix with proper columns
  ncounts <- ncounts(obj)[, c(colA(cc), colB(cc))]
  
  # intersect signature gene names with RNA-seq data
  ligands <- sapply(ligands(sig), function(x) intersect(x, all.genes))
  receptors <- sapply(receptors(sig), function(x) intersect(x, all.genes))
  t.genes <- sapply(tGenes(sig), function(x) intersect(x, all.genes))
  
  good <- sapply(ligands, length) > 0 & sapply(receptors, length) > 0 &
    sapply(t.genes, length) > 0
  ligands   <- ligands[good]
  receptors <- receptors[good]
  t.genes   <- t.genes[good]
  pathways  <- pathways(sig)[good]
  
  # scale ncounts
  if (logTransformed(obj))
    ncounts <- 2**ncounts
  if (robust)
    z <- (ncounts-apply(ncounts,1,stats::median))/apply(ncounts,1,stats::mad)
  else
    z <- (ncounts-rowMeans(ncounts))/apply(ncounts,1,stats::sd)
  if (abs.z.score)
    z <- abs(z)
  
  if( initialOrganism(obj) != "hsapiens" )
    rownames(z) <- all.genes
  
  # compute the LR gene signatures
  i <- NULL
  pwn <- foreach::foreach(i=seq_len(length(pathways)), .combine=c) %do% {
    
    if (name.by.pathway){
      if(rownames.LRP){
        paste0("{",paste(ligands[[i]], collapse=";") ,"} / {",
               paste(receptors[[i]], collapse=";"),"} | ",pathways[[i]])
      }
      else
        pathways[[i]]
    }
    else if (!name.by.pathway){
      paste0("{",paste(ligands[[i]], collapse=";") ,"} / {",
             paste(receptors[[i]], collapse=";"),"}")
    }
    
  }
  
  res <- matrix(0,nrow=length(pathways),ncol=ncol(ncounts),dimnames=list(pwn,colnames(ncounts)))
  for (i in seq_len(length(pathways))){
    
    # average ligand z-score
    zz <- z[ligands[[i]],]
    if (is.matrix(zz))
      mL <- colSums(zz)/length(ligands[[i]])
    else
      mL <- zz
    
    # average receptor z-score
    zz <- z[receptors[[i]],]
    if (is.matrix(zz))
      mR <- colSums(zz)/length(receptors[[i]])
    else
      mR <- zz
    
    # average target gene z-score
    zz <- z[t.genes[[i]],]
    if (is.matrix(zz))
      mT <- colSums(zz)/length(t.genes[[i]])
    else
      mT <- zz
    
    res[i,] <- LR.weight*0.5*(mL+mR)+(1-LR.weight)*mT
  }
  
  res
  
}) # scoreLRGeneSignatures
