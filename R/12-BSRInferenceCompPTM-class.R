library(methods)

#' BulkSignalR cluster comparison-based inference object
#'
#' An S4 class to represent ligand-receptor interactions inferred from
#' a comparison of two clusters of samples. This class inherits from
#' BSRInference.
#'
#' @slot cmp.name  The name of the BSRClusterComp object in a BSRDataModelComp
#' object comp list.
#' @slot tg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot tg.logFC  A list of target gene logFC, one
#' entry per interaction
#' @slot ptm.genes  A list of target gene P-values, one
#' entry per interaction
#' @slot ptmg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot ptmg.corr  A list of target gene P-values, one
#' entry per interaction
#' @slot ptmg.logFC  A list of target gene logFC, one
#' entry per interaction
#' @slot p.genes  A list of target gene P-values, one
#' entry per interaction
#' @slot pg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot pg.logFC  A list of target gene logFC, one
#' entry per interaction
#' @slot dp.genes  A list of target gene P-values, one
#' entry per interaction
#' @slot dpg.pval  A list of target gene P-values, one
#' entry per interaction
#' @slot dpg.corr  A list of target gene P-values, one
#' entry per interaction
#' @slot dpg.logFC  A list of target gene logFC, one
#' entry per interaction
#'
#' @details This class is contains inferred LR interactions along with
#' their statistical significance. Data representation supports subsequent
#' reductions to pathways, etc. See reduction functions
#' \code{"\link[=BSRInferenceCompPTM-class]{reduceToBestPathway}"},
#' \code{"\link[=BSRInferenceCompPTM-class]{reduceToLigand}"},
#' \code{"\link[=BSRInferenceCompPTM-class]{reduceToReceptor}"} and
#' \code{"\link[=BSRInferenceCompPTM-class]{reduceToPathway}"}.
#' @export
#' @examples
#' print('BSRInferenceCompPTM class')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1,"random.example")
#' 
setClass("BSRInferenceCompPTM",
         contains="BSRInference",
         slots=c(cmp.name="character",
                 ptm.genes="list",
                 p.genes="list",
                 dp.genes="list",
                 tg.pval="list",
                 tg.logFC="list",
                 ptmg.pval="list",
                 ptmg.logFC="list",
                 ptmg.corr="list",
                 dpg.corr="list",
                 pg.pval="list",
                 pg.logFC="list",
                 dpg.pval="list",
                 dpg.logFC="list"),
         prototype=list(
           cmp.name="happy",
           LRinter=data.frame(L="A", R="B", pw.id="123",
                              pw.name="one pw", pval=1.0, qval=1.0, L.logFC=2,
                              R.logFC=1.5, LR.pval=0.6, LR.corr=0.5,
                              rank=2, len=50, rank.pval=0.6, rank.corr=0.34,
                              stringsAsFactors=FALSE),
           ptm.genes=list(c("a","b","c")),
           p.genes=list(c("a","b","c")),
           dp.genes=list(c("a","b","c")),
           tg.pval=list(c(0.05,0.1,0.008)),
           tg.logFC=list(c(-1,0,2)),
           ptmg.pval=list(c(0.05,0.1,0.008)),
           ptmg.logFC=list(c(-1,0,2)),
           pg.pval=list(c(0.05,0.1,0.008)),
           pg.logFC=list(c(-1,0,2)),
           dpg.pval=list(c(0.05,0.1,0.008)),
           dpg.logFC=list(c(-1,0,2)),
           ptmg.corr=list(c(0.05,0.1,0.008)),
           dpg.corr=list(c(-1,0,2))
         ))

setValidity("BSRInferenceCompPTM",
            function(object) {
              if (!is.character(object@cmp.name))
                return("cmp.name is not of character type")
              if (length(object@cmp.name) == 0)
                return("cmp.name must have a length > 0")
              if (!is.list(object@tg.pval))
                return("tg.pval is not a list")
              if (!is.list(object@tg.logFC))
                return("tg.logFC is not a list")
              
              TRUE
            }
)

setMethod("show", "BSRInferenceCompPTM",
          function(object) {
            callNextMethod()
            cat("Cluster comparison name:", object@cmp.name, "\n")
          }
)


# Accessors & setters ========================================================

if (!isGeneric("cmpName")) {
  if (is.function("cmpName"))
    fun <- cmpName
  else
    fun <- function(x) standardGeneric("cmpName")
  setGeneric("cmpName", fun)
}
#' Comparison name accessor
#'
#' @name cmpName
#' @aliases cmpName,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("cmpName", "BSRInferenceCompPTM", function(x) x@cmp.name)

if (!isGeneric("cmpName<-")) {
  if (is.function("cmpName<-"))
    fun <- `cmpName<-`
  else
    fun <- function(x,value) standardGeneric("cmpName<-")
  setGeneric("cmpName<-", fun)
}
#' Comparison name setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("cmpName<-", "BSRInferenceCompPTM", function(x, value){
  x@cmp.name <- value
  methods::validObject(x)
  x
})


if (!isGeneric("tgPval")) {
  if (is.function("tgPval"))
    fun <- tgPval
  else
    fun <- function(x) standardGeneric("tgPval")
  setGeneric("tgPval", fun)
}
#' Target gene P-values accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("tgPval", "BSRInferenceCompPTM", function(x) x@tg.pval)

if (!isGeneric("tgPval<-")) {
  if (is.function("tgPval<-"))
    fun <- `tgPval<-`
  else
    fun <- function(x,value) standardGeneric("tgPval<-")
  setGeneric("tgPval<-", fun)
}
#' Target gene P-values setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("tgPval<-", "BSRInferenceCompPTM", function(x, value){
  x@tg.pval <- value
  methods::validObject(x)
  x
})


if (!isGeneric("tgLogFC")) {
  if (is.function("tgLogFC"))
    fun <- tgLogFC
  else
    fun <- function(x) standardGeneric("tgLogFC")
  setGeneric("tgLogFC", fun)
}
#' Target gene logFC accessor
#'
#' @name tgLogFC
#' @aliases tgLogFC,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("tgLogFC", "BSRInferenceCompPTM", function(x) x@tg.logFC)

if (!isGeneric("tgLogFC<-")) {
  if (is.function("tgLogFC<-"))
    fun <- `tgLogFC<-`
  else
    fun <- function(x,value) standardGeneric("tgLogFC<-")
  setGeneric("tgLogFC<-", fun)
}
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("tgLogFC<-", "BSRInferenceCompPTM", function(x, value){
  x@tg.logFC <- value
  methods::validObject(x)
  x
})


if (!isGeneric("pg.pval")) {
  if (is.function("pg.pval"))
    fun <- pg.pval
  else
    fun <- function(x) standardGeneric("pg.pval")
  setGeneric("pg.pval", fun)
}
#' Target gene P-values accessor
#'
#' @name pg.pval
#' @aliases pg.pval,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("pg.pval", "BSRInferenceCompPTM", function(x) x@pg.pval)

if (!isGeneric("pg.pval<-")) {
  if (is.function("pg.pval<-"))
    fun <- `pg.pval<-`
  else
    fun <- function(x,value) standardGeneric("pg.pval<-")
  setGeneric("pg.pval<-", fun)
}
#' Target gene P-values setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("pg.pval<-", "BSRInferenceCompPTM", function(x, value){
  x@pg.pval <- value
  methods::validObject(x)
  x
})


if (!isGeneric("pg.logFC")) {
  if (is.function("pg.logFC"))
    fun <- pg.logFC
  else
    fun <- function(x) standardGeneric("pg.logFC")
  setGeneric("pg.logFC", fun)
}
#' Target gene logFC accessor
#'
#' @name pg.logFC
#' @aliases pg.logFC,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("pg.logFC", "BSRInferenceCompPTM", function(x) x@pg.logFC)

if (!isGeneric("pg.logFC<-")) {
  if (is.function("pg.logFC<-"))
    fun <- `pg.logFC<-`
  else
    fun <- function(x,value) standardGeneric("pg.logFC<-")
  setGeneric("pg.logFC<-", fun)
}
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("pg.logFC<-", "BSRInferenceCompPTM", function(x, value){
  x@pg.logFC <- value
  methods::validObject(x)
  x
})



if (!isGeneric("ptmg.pval")) {
  if (is.function("ptmg.pval"))
    fun <- ptmg.pval
  else
    fun <- function(x) standardGeneric("ptmg.pval")
  setGeneric("ptmg.pval", fun)
}
#' Target gene P-values accessor
#'
#' @name ptmg.pval
#' @aliases ptmg.pval,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("ptmg.pval", "BSRInferenceCompPTM", function(x) x@ptmg.pval)

if (!isGeneric("ptmg.pval<-")) {
  if (is.function("ptmg.pval<-"))
    fun <- `ptmg.pval<-`
  else
    fun <- function(x,value) standardGeneric("ptmg.pval<-")
  setGeneric("ptmg.pval<-", fun)
}
#' Target gene P-values setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("ptmg.pval<-", "BSRInferenceCompPTM", function(x, value){
  x@ptmg.pval <- value
  methods::validObject(x)
  x
})


if (!isGeneric("ptmg.logFC")) {
  if (is.function("ptmg.logFC"))
    fun <- ptmg.logFC
  else
    fun <- function(x) standardGeneric("ptmg.logFC")
  setGeneric("ptmg.logFC", fun)
}
#' Target gene logFC accessor
#'
#' @name ptmg.logFC
#' @aliases ptmg.logFC,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("ptmg.logFC", "BSRInferenceCompPTM", function(x) x@ptmg.logFC)

if (!isGeneric("ptmg.logFC<-")) {
  if (is.function("ptmg.logFC<-"))
    fun <- `ptmg.logFC<-`
  else
    fun <- function(x,value) standardGeneric("ptmg.logFC<-")
  setGeneric("ptmg.logFC<-", fun)
}
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("ptmg.logFC<-", "BSRInferenceCompPTM", function(x, value){
  x@ptmg.logFC <- value
  methods::validObject(x)
  x
})



if (!isGeneric("dpg.pval")) {
  if (is.function("dpg.pval"))
    fun <- dpg.pval
  else
    fun <- function(x) standardGeneric("dpg.pval")
  setGeneric("dpg.pval", fun)
}
#' Target gene P-values accessor
#'
#' @name dpg.pval
#' @aliases dpg.pval,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("dpg.pval", "BSRInferenceCompPTM", function(x) x@dpg.pval)

if (!isGeneric("dpg.pval<-")) {
  if (is.function("pdgPval<-"))
    fun <- `dpg.pval<-`
  else
    fun <- function(x,value) standardGeneric("dpg.pval<-")
  setGeneric("dpg.pval<-", fun)
}
#' Target gene P-values setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("dpg.pval<-", "BSRInferenceCompPTM", function(x, value){
  x@dpg.pval <- value
  methods::validObject(x)
  x
})


if (!isGeneric("dpg.logFC")) {
  if (is.function("dpg.logFC"))
    fun <- dpg.logFC
  else
    fun <- function(x) standardGeneric("dpg.logFC")
  setGeneric("dpg.logFC", fun)
}
#' Target gene logFC accessor
#'
#' @name dpg.logFC
#' @aliases dpg.logFC,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("dpg.logFC", "BSRInferenceCompPTM", function(x) x@dpg.logFC)

if (!isGeneric("dpg.logFC<-")) {
  if (is.function("dpg.logFC<-"))
    fun <- `dpg.logFC<-`
  else
    fun <- function(x,value) standardGeneric("dpg.logFC<-")
  setGeneric("dpg.logFC<-", fun)
}
#' Target gene logFC setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("dpg.logFC<-", "BSRInferenceCompPTM", function(x, value){
  x@dpg.logFC <- value
  methods::validObject(x)
  x
})



if (!isGeneric("ptm.genes")) {
  if (is.function("ptm.genes"))
    fun <- ptm.genes
  else
    fun <- function(x) standardGeneric("ptm.genes")
  setGeneric("ptm.genes", fun)
}
#' PTM genes accessor
#'
#' @name ptm.genes
#' @aliases ptm.genes,BSRInferenceCompPTM-method
#' @param x BSRInferance object
#' @export
setMethod("ptm.genes", "BSRInferenceCompPTM", function(x) x@ptm.genes)

if (!isGeneric("ptm.genes<-")) {
  if (is.function("ptm.genes<-"))
    fun <- `ptm.genes<-`
  else
    fun <- function(x,value) standardGeneric("ptm.genes<-")
  setGeneric("ptm.genes<-", fun)
}
#' PTM genes setter (internal use only)
#' @param x BSRInferance object
#' @param value value to be set BSRInferenceCompPTM
#' @keywords internal 
setMethod("ptm.genes<-", "BSRInferenceCompPTM", function(x, value){
  x@ptm.genes <- value
  methods::validObject(x)
  x
})

if (!isGeneric("ptmg.corr")) {
  if (is.function("ptmg.corr"))
    fun <- ptmg.corr
  else
    fun <- function(x) standardGeneric("ptmg.corr")
  setGeneric("ptmg.corr", fun)
}
#' Target gene correlations accessor
#'
#' @name ptmg.corr
#' @aliases ptmg.corr,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("ptmg.corr", "BSRInferenceCompPTM", function(x) x@ptmg.corr)

if (!isGeneric("ptmg.corr<-")) {
  if (is.function("ptmg.corr<-"))
    fun <- `ptmg.corr<-`
  else
    fun <- function(x,value) standardGeneric("ptmg.corr<-")
  setGeneric("ptmg.corr<-", fun)
}
#' PTM gene correlations setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("ptmg.corr<-", "BSRInferenceCompPTM", function(x, value){
  x@ptmg.corr <- value
  methods::validObject(x)
  x
})


if (!isGeneric("p.genes")) {
  if (is.function("p.genes"))
    fun <- p.genes
  else
    fun <- function(x) standardGeneric("p.genes")
  setGeneric("p.genes", fun)
}
#' PTM genes accessor
#'
#' @name p.genes
#' @aliases p.genes,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("p.genes", "BSRInferenceCompPTM", function(x) x@p.genes)

if (!isGeneric("p.genes<-")) {
  if (is.function("p.genes<-"))
    fun <- `p.genes<-`
  else
    fun <- function(x,value) standardGeneric("p.genes<-")
  setGeneric("p.genes<-", fun)
}
#' PTM genes setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set BSRInferenceCompPTM
#' @keywords internal 
setMethod("p.genes<-", "BSRInferenceCompPTM", function(x, value){
  x@p.genes <- value
  methods::validObject(x)
  x
})


if (!isGeneric("pg.corr")) {
  if (is.function("pg.corr"))
    fun <- pg.corr
  else
    fun <- function(x) standardGeneric("pg.corr")
  setGeneric("pg.corr", fun)
}
#' Target gene correlations accessor
#'
#' @name pg.corr
#' @aliases pg.corr,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("pg.corr", "BSRInferenceCompPTM", function(x) x@pg.corr)

if (!isGeneric("pg.corr<-")) {
  if (is.function("pg.corr<-"))
    fun <- `pg.corr<-`
  else
    fun <- function(x,value) standardGeneric("pg.corr<-")
  setGeneric("pg.corr<-", fun)
}
#' PTM gene correlations setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("pg.corr<-", "BSRInferenceCompPTM", function(x, value){
  x@pg.corr <- value
  methods::validObject(x)
  x
})


if (!isGeneric("dp.genes")) {
  if (is.function("dp.genes"))
    fun <- dp.genes
  else
    fun <- function(x) standardGeneric("dp.genes")
  setGeneric("dp.genes", fun)
}
#' PTM genes accessor
#'
#' @name dp.genes
#' @aliases dp.genes,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("dp.genes", "BSRInferenceCompPTM", function(x) x@dp.genes)

if (!isGeneric("dp.genes<-")) {
  if (is.function("dp.genes<-"))
    fun <- `dp.genes<-`
  else
    fun <- function(x,value) standardGeneric("dp.genes<-")
  setGeneric("dp.genes<-", fun)
}
#' PTM genes setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set BSRInferenceCompPTM
#' @keywords internal 
setMethod("dp.genes<-", "BSRInferenceCompPTM", function(x, value){
  x@dp.genes <- value
  methods::validObject(x)
  x
})

if (!isGeneric("dpg.corr")) {
  if (is.function("dpg.corr"))
    fun <- dpg.corr
  else
    fun <- function(x) standardGeneric("dpg.corr")
  setGeneric("dpg.corr", fun)
}
#' Target gene correlations accessor
#'
#' @name dpg.corr
#' @aliases dpg.corr,BSRInferenceCompPTM-method
#' @param x BSRInferenceCompPTM object
#' @export
setMethod("dpg.corr", "BSRInferenceCompPTM", function(x) x@dpg.corr)

if (!isGeneric("dpg.corr<-")) {
  if (is.function("dpg.corr<-"))
    fun <- `dpg.corr<-`
  else
    fun <- function(x,value) standardGeneric("dpg.corr<-")
  setGeneric("dpg.corr<-", fun)
}
#' PTM gene correlations setter (internal use only)
#' @param x BSRInferenceCompPTM object
#' @param value value to be set for bsrinf
#' @keywords internal
setMethod("dpg.corr<-", "BSRInferenceCompPTM", function(x, value){
  x@dpg.corr <- value
  methods::validObject(x)
  x
})

# Rescoring ====================================================================

if (!isGeneric("rescoreInference")) {
  if (is.function("rescoreInference"))
    fun <- rescoreInference
  else
    fun <- function(obj, ...) standardGeneric("rescoreInference")
  setGeneric("rescoreInference", fun)
}
#' Inference re-scoring
#'
#' A method to re-score an existing BSRInferenceCompPTM object
#' (P- and Q-value estimations).
#'
#' @name rescoreInference
#' @aliases rescoreInference,BSRInferenceCompPTM-method
#'
#' @param obj BSRInferecenceComp object.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @details A BSRInferenceCompPTM object should be created by calling 
#' \code{"\link[=BSRClusterComp-class]{initialInference}"}
#'
#' @return A BSRInferenceCompPTM object.
#'
#' @export
#' @examples
#' print('rescoreInference')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # rescore
#' bsrinf.less <- rescoreInference(bsrinf, param=param(bsrdm.comp), rank.p=0.75)
#'
setMethod("rescoreInference", "BSRInferenceCompPTM", function(obj, param, rank.p=0.55,
            fdr.proc=c("BH", "Bonferroni", "Holm",
                       "Hochberg", "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {
  
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  
  # extract the necessary data from the BSRInferenceCompPTM object
  pairs <- LRinter(obj)
  t.genes <- tGenes(obj)
  tg.pval <- tgPval(obj)
  tg.corr <- tgCorr(obj)
  ptm.genes <- ptm.genes(obj)
  ptmg.pval <- ptmg.pval(obj)
  ptmg.corr <- ptmg.corr(obj)
  p.genes <- p.genes(obj)
  pg.pval <- pg.pval(obj)
  pg.corr <- pg.corr(obj)
  dp.genes <- dp.genes(obj)
  dpg.pval <- dpg.pval(obj)
  dpg.corr <- dpg.corr(obj)
  
  # recompute P-values
  for (i in 1:nrow(pairs)){
    tg <- t.genes[[i]]
    spvals <- tg.pval[[i]]
    
    # get the LR correlation P-value
    p.lr <- pairs$LR.pval[i]
    
    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    len <- pairs$len[i]
    r <- min(max(1, trunc(rank.p*len)), len)
    rank.pval <- spvals[r]
    # r-1 P-values are > rank.pval, prob to have r-1 or less
    # P-values > rank.pval is given by a binomial with success rate
    # equal to the probability to get a P-value > rank.pval, i.e.,
    # 1-rank.pval. If rank.pval is low (i.e., highly significant),
    # it becomes difficult to get as little as r-1 P-values > rank.pval by chance!
    p.rt <- stats::pbinom(r-1, len, 1-rank.pval) # cdf is punif here!
    pairs$pval[i] <- p.lr*p.rt
    pairs$rank.pval[i] <- rank.pval
    pairs$rank.corr[i] <- tg.corr[[i]][r]
  }
  
  # recompute the Q-values
  rawp <- pairs$pval
  adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
  pairs$qval <- adj$adjp[order(adj$index),fdr.proc]
  
  # update the BSRInference object
  inf.param <- infParam(obj)
  inf.param$fdr.proc <- fdr.proc
  inf.param$rank.p <- rank.p
  infParam(obj) <- inf.param
  LRinter(obj) <- pairs
  
  obj
  
}) # rescoreInference


# Reduction and pathway stat methods ===========================================


if (!isGeneric("reduceToBestPathway")) {
  if (is.function("reduceToBestPathway"))
    fun <- reduceToBestPathway
  else
    fun <- function(obj, ...) standardGeneric("reduceToBestPathway")
  setGeneric("reduceToBestPathway", fun)
}
#' Keep one pathway per ligand-receptor pair
#'
#' @name reduceToBestPathway
#' @aliases reduceToBestPathway,BSRInferenceCompPTM-method
#'
#' @param obj BSRInferenceCompPTM object
#'
#' @return A BSRInferenceCompPTM object reduced to only report one pathway per
#' ligand-receptor pair. The pathway with the
#' smallest P-value is selected.
#'
#' @details Ligand-receptor pairs
#' are evaluated in relation with pathways that allow checking receptor
#' downstream correlations. It is thus possible
#' that several pathways are reported for a same LR pair.
#'
#' @export
#' @examples
#' print('reduceToBestPathway')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redBP  <- reduceToBestPathway(bsrinf)
#'
#' @importFrom rlang .data
setMethod("reduceToBestPathway", "BSRInferenceCompPTM", function(obj) {
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual data representation
  
  # get best p-value pathway per LR pair
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  LRinter <- NULL
  pairs <- obj@LRinter
  LR <- unique(pairs[, c("L","R")])
  for (i in seq_len(nrow(LR))){
    
    L <- LR$L[i]
    R <- LR$R[i]
    
    pwr <- pairs[pairs$L==L & pairs$R==R,]
    k <- which.min(pwr$pval)
    j <- which(pairs$L==L & pairs$R==R & pairs$pw.id==pwr$pw.id[k])
    ligands <- c(ligands, obj@ligands[j])
    receptors <- c(receptors, obj@receptors[j])
    t.genes <- c(t.genes, obj@t.genes[j])
    tg.corr <- c(tg.corr, obj@tg.corr[j])
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    LRinter <- rbind(LRinter, pairs[j,])
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@inf.param$pathway.reduced <- TRUE
  
  obj
  
}) # reduceToBestPathway


if (!isGeneric("reduceToReceptor")) {
  if (is.function("reduceToReceptor"))
    fun <- reduceToReceptor
  else
    fun <- function(obj, ...) standardGeneric("reduceToReceptor")
  setGeneric("reduceToReceptor", fun)
}
#' Aggregate the ligands of a same receptor
#'
#' Simplifies a ligand-receptor table to focus on the receptors.
#'
#' @name reduceToReceptor
#' @aliases reduceToReceptor,BSRInferenceCompPTM-method
#'
#' @return BSRInferenceCompPTM object reduced to one row per receptor.
#' All the ligands are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the line with the
#' pathway featuring the smallest P-value.
#' @param obj BRSInferenceComp object
#' @export
#' @examples
#' print('reduceToReceptor')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redR  <- reduceToReceptor(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToReceptor", "BSRInferenceCompPTM", function(obj){
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual data representation
  
  if (infParam(obj)$ligand.reduced)
    stop("Already reduced to receptor") # because ligands were reduced
  
  # pool the ligands
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  LRinter <- NULL
  pairs <- obj@LRinter
  for (R in unique(pairs$R)){
    lig <- pairs[pairs$R==R,]
    k <- which.min(lig$pval)
    j <- which(pairs$R==lig$R[k] & pairs$pw.id==lig$pw.id[k])[1]
    ligands <- c(ligands, list(unique(lig$L)))
    receptors <- c(receptors, list(R))
    t.genes <- c(t.genes, obj@t.genes[j])
    tg.corr <- c(tg.corr, obj@tg.corr[j])
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    to.add <- pairs[j,]
    to.add[1, "L"] <- paste0("{", paste(unique(lig$L), collapse=";"), "}")
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@inf.param$pathway.reduced <- TRUE
  obj@inf.param$ligand.reduced <- TRUE
  
  obj
  
}) # reduceToReceptor


if (!isGeneric("reduceToLigand")) {
  if (is.function("reduceToLigand"))
    fun <- reduceToLigand
  else
    fun <- function(obj, ...) standardGeneric("reduceToLigand")
  setGeneric("reduceToLigand", fun)
}
#' Aggregate the receptors of a same ligand
#'
#' Simplifies a ligand-receptor table to focus on the ligands.
#'
#' @name reduceToLigand
#' @aliases reduceToLigand,BSRInferenceCompPTM-method
#'
#' @return A BSRInferenceCompPTM object but reduced to one row per ligand.
#' All the receptors are combined in a
#' semi-colon-separated list surrounded by curly brackets in the tabular
#' slot \code{LRinter}, and in vectors in the \code{ligands} (list) slot.
#'
#' The reported P-value and target genes are those from the pathway with
#' the smallest P-value.
#' @param obj BSRInferenceCompPTM object
#' @export
#' @examples
#' print('reduceToLigand')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redL  <- reduceToLigand(bsrinf)  
#'
#' @importFrom rlang .data
setMethod("reduceToLigand", "BSRInferenceCompPTM", function(obj){
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual representation
  
  if (infParam(obj)$receptor.reduced)
    stop("Already reduced to ligand") # because receptors were reduced
  
  # pool the receptors
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  LRinter <- NULL
  pairs <- obj@LRinter
  for (L in unique(pairs$L)){
    rec <- pairs[pairs$L==L,]
    k <- which.min(rec$pval)
    j <- which(pairs$L==L & pairs$R==rec$R[k] & pairs$pw.id==rec$pw.id[k])
    ligands <- c(ligands, list(L))
    receptors <- c(receptors, list(unique(rec$R)))
    t.genes <- c(t.genes, obj@t.genes[j])
    tg.corr <- c(tg.corr, obj@tg.corr[j])
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    to.add <- pairs[j,]
    to.add[1, "R"] <- paste0("{", paste(unique(rec$R), collapse=";"), "}")
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@inf.param$pathway.reduced <- TRUE
  obj@inf.param$receptor.reduced <- TRUE
  
  obj
  
}) # reduceToLigand


if (!isGeneric("reduceToPathway")) {
  if (is.function("reduceToPathway"))
    fun <- reduceToPathway
  else
    fun <- function(obj, ...) standardGeneric("reduceToPathway")
  setGeneric("reduceToPathway", fun)
}
#' Aggregate ligands and receptors at the pathway level
#'
#' Simplifies a ligand-receptor inference object to focus on
#' the pathways.
#'
#' @name reduceToPathway
#' @aliases reduceToPathway,BSRInferenceCompPTM-method
#'
#' @return A BSRInferenceCompPTM object reduced to only report one row per pathway.
#' The information of which ligand interacted with which receptor is lost as
#' all the ligands and all the receptors forming pairs related to a certain
#' pathway are combined.
#' For a given pathway, the reported P-values and target genes are those of
#' the best ligand-receptor pair that
#' was in this pathway.
#' Receptors and ligands are combined in two semi-colon-separated
#' lists surrounded by curly brackets in the tabular slot \code{LRinter},
#' while the list representation slots (\code{ligands} and
#' \code{receptors}) are update accordingly.
#'
#' @param obj BSRInferenceCompPTM object
#' @export
#' @examples
#' print('reduceToPathway')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reduction
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' @importFrom rlang .data
setMethod("reduceToPathway", "BSRInferenceCompPTM", function(obj){
  
  # Here we access the object slots directly as this procedure
  # is dependent of actual representation
  
  if (infParam(obj)$receptor.reduced)
    stop("Already reduced to ligand") # because receptors were reduced
  if (infParam(obj)$ligand.reduced)
    stop("Already reduced to receptor") # because ligands were reduced
  
  # reduce to unique pathways
  ligands <- list()
  receptors <- list()
  t.genes <- list()
  tg.corr <- list()
  tg.pval <- list()
  tg.logFC <- list()
  LRinter <- NULL
  pairs <- obj@LRinter
  for (p in unique(pairs$pw.id)){
    j <- which(pairs$pw.id==p)[1]
    ligands <- c(ligands, list(unique(pairs$L[pairs$pw.id==p])))
    receptors <- c(receptors, list(unique(pairs$R[pairs$pw.id==p])))
    t.genes <- c(t.genes, obj@t.genes[j])
    tg.corr <- c(tg.corr, obj@tg.corr[j])
    tg.pval <- c(tg.pval, obj@tg.pval[j])
    tg.logFC <- c(tg.logFC, obj@tg.logFC[j])
    to.add <- pairs[j,]
    to.add[1, "L"] <- paste0("{", paste(unique(pairs$L[pairs$pw.id==p]),
                                        collapse=";"), "}")
    to.add[1, "R"] <- paste0("{", paste(unique(pairs$R[pairs$pw.id==p]),
                                        collapse=";"), "}")
    LRinter <- rbind(LRinter, to.add)
  }
  
  # update the object
  obj@LRinter <- LRinter
  obj@ligands <- ligands
  obj@receptors <- receptors
  obj@t.genes <- t.genes
  obj@tg.corr <- tg.corr
  obj@tg.pval <- tg.pval
  obj@tg.logFC <- tg.logFC
  obj@inf.param$ligand.reduced <- TRUE
  obj@inf.param$receptor.reduced <- TRUE
  
  obj
  
}) # reduceToPathway


# Obtain gene signatures from a BSRInference object ============================

if (!isGeneric("getLRGeneSignatures")) {
  if (is.function("getLRGeneSignatures"))
    fun <- getLRGeneSignatures
  else
    fun <- function(obj, ...) standardGeneric("getLRGeneSignatures")
  setGeneric("getLRGeneSignatures", fun)
}
#' Extract gene signatures of LR pair activity
#'
#' Obtains gene signatures reflecting ligand-receptor as well as
#' receptor downstream activity to
#' score ligand-receptor pairs across samples subsequently with
#' \code{"\link[=BSRInferenceCompPTM-class]{scoreLRGeneSignatures}"}
#'
#' @name getLRGeneSignatures
#' @aliases getLRGeneSignatures,BSRInferenceCompPTM-method
#'
#' @param obj    BSRInferenceCompPTM object.
#' @param pval.thres    P-value threshold.
#' @param qval.thres    Q-value threshold.
#' @param with.pw.id    A logical indicating whether the ID of a pathway
#' should be concatenated to its name.
#' @return A BSRSignatureComp object containing a gene signature for each triple
#' ligand-receptor pair. A reduction to the best pathway
#' for each pair is automatically performed and the gene signature is
#' comprised of the ligand, the receptor,
#' and all the target genes with rank equal or superior to \code{pairs$rank}.
#' @export
#' @examples
#' print('getLRGeneSignatures')
#' 
#' # prepare data
#' data(sdc,package='BulkSignalR')
#' normal <- grep("^N", names(sdc))
#' bsrdm <- prepareDataset(sdc[,-normal])
#' 
#' # define the comparison
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:5)
#' colB <- as.integer(8:15)
#' n <- nrow(ncounts(bsrdm.comp))
#' stats <- data.frame(pval=runif(n), logFC=rnorm(n, 0, 2))
#' rownames(stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "random.example")
#' 
#' # infer ligand-receptor interactions from the comparison
#' bsrinf <- initialInference(bsrdm.comp,max.pval=1, "random.example")
#' 
#' # reductions
#' bsrinf.redP  <- reduceToPathway(bsrinf)  
#' bsrsig.redP <- getLRGeneSignatures(bsrinf.redP,qval.thres=0.001)
#'
#' @importFrom foreach %do% %dopar%
#' @importFrom methods new
setMethod("getLRGeneSignatures", "BSRInferenceCompPTM", function(obj,
                                                          pval.thres=NULL, qval.thres=NULL, with.pw.id=FALSE){
  
  if (is.null(pval.thres) && is.null(qval.thres))
    stop("Either a P- or a Q-value threshold must be provided")
  
  # reduce and select
  obj <- reduceToBestPathway(obj)
  pairs <- LRinter(obj)
  if (!is.null(pval.thres))
    selected <- pairs$pval <= pval.thres
  else
    selected <- pairs$qval <= qval.thres
  
  # obtain the signature object
  pairs <- pairs[selected,]
  ligands <- ligands(obj)[selected]
  receptors <- receptors(obj)[selected]
  if (with.pw.id)
    pathways <- paste0(pairs$pw.id, ": ", pairs$pw.name)
  else
    pathways <- pairs$pw.name
  t.genes <- tGenes(obj)[selected]
  t.corrs <- tgCorr(obj)[selected]
  t.pvals <- tgPval(obj)[selected]
  t.logFCs <- tgLogFC(obj)[selected]
  
  for (i in seq_len(nrow(pairs))){
    tg <- t.genes[[i]]
    t.genes[[i]] <- tg[pairs$rank[i]:length(tg)]
    tc <- t.corrs[[i]]
    t.corrs[[i]] <- tc[pairs$rank[i]:length(tc)]
    tp <- t.pvals[[i]]
    t.pvals[[i]] <- tp[pairs$rank[i]:length(tp)]
    tl <- t.logFCs[[i]]
    t.logFCs[[i]] <- tl[pairs$rank[i]:length(tl)]
  }
  
  new("BSRSignatureComp", pathways=pathways, ligands=ligands,
      receptors=receptors, t.genes=t.genes, tg.corr=t.corrs,
      tg.pval=t.pvals, tg.logFC=t.logFCs,
      cmp.name=cmpName(obj))
  
}) # getLRGeneSignatures


