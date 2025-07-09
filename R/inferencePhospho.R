#' Get correlated ligand-receptor pairs.
#'
#' Internal function to compute the Spearman correlations
#' of all the ligand-receptor
#' pairs in LRdb and return those above a minimum value.
#'
#' @param ds              A BSRDataModel object.
#' @param min.cor         The minimum correlation required.
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#'
#' @return A data frame containing putative ligand-receptor pairs along
#'   with their correlations above \code{min.cor}. This table is the first step
#'   of a ligand-receptor analysis.
#'
#'
#' @details The \code{restrict.genes} parameter is used for special cases where
#'   LRdb must be further restricted to a subset.
#'   The putative ligand-receptor pairs has 3 columns : R, L and corr.
#'
#' @importFrom methods is
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @keywords internal
.getCorrelatedLR <- function(ds, min.cor = 0.25, restrict.genes = NULL) {
  
  # local binding
  i <- NULL
  
  cat("\n .getCorrelatedLR - inferencePTM 2")
  if ((min.cor < -1) || (min.cor > 1))
    stop("min.cor must lie in [-1;+1]")
  # if (!is(ds, "BSRDataModel"))
  #   stop("ds must be an object of class BSRDataModel")
  
  #cat("\n", ds@single, "\n")
  if(ds@single == TRUE){
    rnNcountsDs <- unique(symPos(ds)[,1])
  }
  else{
    rnNcountsDs <- rownames(ncounts(ds))
  }
  #cat("\n rnC: ", head(rnNcountsDs), "\n")
  lrgenes <- unique(intersect(c(LRdb$ligand,
                                LRdb$receptor), rnNcountsDs))
  # if(!is.null(lrgenes))
  #   cat("\n lrG: ", head(lrgenes), "\n") #ACTR2
  # lrgenes <- intersect(c(LRdb$ligand,
  #                        LRdb$receptor), rownames(ncounts(ds)))
  if (!is.null(restrict.genes))
    lrgenes <- intersect(lrgenes, restrict.genes)
  
  # cat("\n ncounts la: ", ncounts(ds)[1:3,1:3], "\n")
  # cat("\n lncounts: ")
  # cat(rownames(ncounts(ds))[1:3]) #SCYL3_533
  # cat("\n")
  # cat("\n cncounts: ")
  # cat(colnames(ncounts(ds))[1:3])
  # cat("\n")
  # cat(sum(lrgenes %in% rownames(ncounts(ds)))) #pb (0)
  # compute all the correlations at once
  #cat(ncounts(ds)[lrgenes[lrgenes %in% rownames(ncounts(ds))],]) #whaaat
  
  # get the pairs
  cat("\n get the pairs \n")
  if(ds@single == TRUE && grepl("_", unlist(rownames(ncounts(ds))), fixed = TRUE)){
    ## test grepl("_", unlist(rownames(ncounts(ds))), fixed = TRUE) !!!
    cat("oups \n")
    posRNC <- which(symPos(ds)[,1] %in% lrgenes)
    lrgenes2 <- paste0(symPos(ds)[posRNC,1], "_", symPos(ds)[posRNC,2])
    corlr <- stats::cor(t(ncounts(ds)[lrgenes2[lrgenes2 %in% rownames(ncounts(ds))], ]), method = "spearman")
    cc <- strsplit(rownames(corlr),'_')
    part1 <- unlist(cc)[2*(1:length(rownames(corlr)))-1]
    cat("\n oups1 \n")
    cat(length(part1))
    cat("\n part1 \n")
    cat(part1[1:3]) #CFTR_660
    corlr2 <- corlr
    cat("\n oups4 \n")
    cat(dim(corlr2))
    #corlr2[,ncol(corlr2)+1] <- part1
    corlr2 <- cbind(corlr2, part1)
    cat("\n oups5 \n")
    cat(dim(corlr2))
    cat("\n ouuuu \n")
    dup <- which(!duplicated(corlr2[,ncol(corlr2)]))
    corlr2 <- corlr2[dup,c(dup,ncol(corlr2))]
    cat("\n ouuuu2 \n")
    tmp2 <- corlr2[,ncol(corlr2)]
    corlr2 <- corlr2[,-ncol(corlr2)]
    cat("\n oups3 \n")
    cat(dim(corlr2))
    cat("\n oupsb \n")
    rownames(corlr2) <- colnames(corlr2) <- tmp2
    cat(rownames(corlr2)[1:3])
    #cat("\n oupsc \n")
    #cat(sum(lrgenes %in% rownames(ncounts(ds))))
    
    pairs <- foreach::foreach(i = seq_len(nrow(LRdb)),
                              .combine = rbind) %do% {
                                if (LRdb$ligand[i] %in% part1 &&
                                    LRdb$receptor[i] %in% part1)
                                  data.frame(L = LRdb$ligand[i],
                                             R = LRdb$receptor[i],
                                             corr = corlr2[LRdb$ligand[i],
                                                           LRdb$receptor[i]],
                                             stringsAsFactors = FALSE)
                                else
                                  NULL
                              }
  }
  else{
    corlr <- stats::cor(t(ncounts(ds)[lrgenes[lrgenes %in% rownames(ncounts(ds))], ]), method = "spearman")
    pairs <- foreach::foreach(i = seq_len(nrow(LRdb)),
                              .combine = rbind) %do% {
                                if (LRdb$ligand[i] %in% rownames(corlr) &&
                                    LRdb$receptor[i] %in% rownames(corlr))
                                  data.frame(L = LRdb$ligand[i],
                                             R = LRdb$receptor[i],
                                             corr = corlr[LRdb$ligand[i],
                                                          LRdb$receptor[i]],
                                             stringsAsFactors = FALSE)
                                else
                                  NULL
                              }
  }
  
  # cat("pairsLR \n")
  # cat(unlist(pairs)[1:3])
  good <- pairs$corr >= min.cor
  pairs[good,]
  
}  # .getCorrelatedLR


#' Get correlated receptor-PTMprotein pairs. (sert a rien ?)
#'
#' Internal function to compute the Spearman correlations
#' of all the receptor-PTMprotein
#' pairs in LRdb and return those above a minimum value.
#'
#' @param ds              A BSRDataModel object.
#' @param min.cor         The minimum correlation required.
#' @param restrict.genes  A list of gene symbols that restricts PTMprotein and
#'   receptors.
#'
#' @return A data frame containing putative receptor-PTMprotein pairs along
#'   with their correlations above \code{min.cor}. This table is the first step
#'   of a receptor-PTMprotein analysis.
#'
#'
#' @details The \code{restrict.genes} parameter is used for special cases where
#'   LRdb must be further restricted to a subset.
#'   The putative receptor-PTMprotein pairs has 3 columns : R, P and corr.
#'
#' @importFrom methods is
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @keywords internal
.getCorrelatedRP <- function(ds, min.cor = 0.25, restrict.genes = NULL) {
  cat("\n .getCorrelatedRP - inferencePTM ")
  # local binding
  i <- NULL
  
  if ((min.cor < -1) || (min.cor > 1))
    stop("min.cor must lie in [-1;+1]")
  # if (!is(ds, "BSRDataModel"))
  #   stop("ds must be an object of class BSRDataModel")
  
  if(ds@single == TRUE){
    rnNcountsDs <- unique(symPos(ds)[,1])
  }
  else{
    rnNcountsDs <- rownames(ncounts(ds))
  }
  
  lrgenes <- intersect(c(LRdb$ligand,
                         LRdb$receptor), rnNcountsDs)
  if (!is.null(restrict.genes))
    lrgenes <- intersect(lrgenes, restrict.genes)
  
  # compute all the correlations at once
  corlr <- stats::cor(t(ncounts(ds)[lrgenes, ]), method = "spearman")
  
  # get the pairs
  pairs <- foreach::foreach(i = seq_len(nrow(LRdb)),
                            .combine = rbind) %do% {
                              if (LRdb$ligand[i] %in% rownames(corlr) &&
                                  LRdb$receptor[i] %in% rownames(corlr))
                                data.frame(L = LRdb$ligand[i],
                                           R = LRdb$receptor[i],
                                           corr = corlr[LRdb$ligand[i],
                                                        LRdb$receptor[i]],
                                           stringsAsFactors = FALSE)
                              else
                                NULL
                            }
  cat("pairsRP \n")
  cat(unlist(pairs)[1:3])
  good <- pairs$corr >= min.cor
  pairs[good,]
  
}  # .getCorrelatedRP


#' Internal function to check receptor signaling downstream
#'
#' @param lr              A data frame as returned by
#'   \code{.getCorrelatedLR()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway
#'   IDs).
#' @param rncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param id.col          Column index or name in \code{pw} for the pathway IDs.
#' @param gene.col        Column index or name in \code{pw} for the gene
#'   symbols.
#' @param pw.col          Column index or name in \code{pw} for the pathway
#'   names.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#'
#' @return A table reporting all the receptor-PTMprotein pairs provided in \code{lr}
#'   along with the pathways found and data about target gene correlations with
#'   the receptor.
#'
#' @importFrom foreach %do% %dopar%
#' @keywords internal

#on s'en fout en PTM !!
.downstreamSignaling <- function(lr, pw, pw.size, rncounts, PTM, id.col, gene.col,
                                 pw.col, min.positive=4, with.complex = TRUE) {
  if (!is.matrix(rncounts))
    stop("rncounts must be a matrix")
  
  # local binding
  r <- p <- pl <- id <- NULL
  
  cat("\n .downstreamSignaling - inferencePTM ")
  # define interaction types
  control.int <- "controls-expression-of"
  incomplex.int <- c("in-complex-with","interacts-with")
  directed.int <- c("controls-state-change-of", "catalysis-precedes",
                    "controls-expression-of", "controls-transport-of",
                    "controls-phosphorylation-of")
  PTMph.int <- "controls-phosphorylation-of"
  if (with.complex)
    correlated.int <- union(control.int, incomplex.int)
  else
    correlated.int <- control.int
  
  # compute downstream correlations
  corrg <- stats::cor(t(rncounts), method = "spearman")
  
  # the global computation above is faster than restricted to the receptors
  corrg <- corrg[unique(lr$R), ]
  #cat("\n a")
  n<-0
  # check each putative LR pair, loop over the receptors
  reg.proc <- foreach::foreach(r=unique(lr$R),.combine=rbind) %do% {
    # reg.proc <- NULL
    # for (r in unique(lr$putative.pairs$R)){
    n<-n+1
    # if(n==1 | n%%100==0)
    # cat("\n a1")
    # loop over the pathways containing the receptor r
    pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
    if (length(pa)>0){
      receptor.ligands <- unique(lr$L[lr$R==r])
      # if(n==1 | n%%100==0)
      # cat("\n a2")
      best.2nd <- foreach::foreach(p=pa,.combine=rbind) %do% {
        # best.2nd <- NULL
        # for (p in pa){
        int <- SingleCellSignalR::PwC_ReactomeKEGG[
          SingleCellSignalR::PwC_ReactomeKEGG$a.gn %in% pw[pw[[id.col]]==p,gene.col] &
            SingleCellSignalR::PwC_ReactomeKEGG$b.gn %in% pw[pw[[id.col]]==p,gene.col],
        ]
        directed <- int$type %in% directed.int
        
        # double the undirected interactions and generate a directed graph
        ret <- int[!directed,c("a.gn", "b.gn")]
        from <- ret$a.gn
        ret$a.gn <- ret$b.gn
        ret$b.gn <- from
        d.int <- unique(rbind(int[,c("a.gn", "b.gn")],ret))
        g <- igraph::graph_from_data_frame(d.int, directed=TRUE)
        # if(n==1 | n%%100==0)
        # cat("\n a3")
        # extract the target genes of receptor r
        if (r %in% d.int$a.gn || r %in% d.int$b.gn){
          # putative targets in the pathway
          target.genes <- setdiff(c(
            int[int$type %in% correlated.int & int$a.gn==r, "b.gn"],
            int[int$type %in% correlated.int & int$b.gn==r, "a.gn"],
            int[int$type %in% directed.int, "b.gn"]),
            r
          )
          
          # reduce putative to reachable from the receptor
          sp <- igraph::shortest.paths(g, r, target.genes)
          target.genes <- colnames(sp)[!is.infinite(sp[r,])]
          if(n==1 | n%%100==0)
          cat("\n b")
          # eliminate ligands of the receptor if present
          target.genes <- setdiff(target.genes, receptor.ligands)
          
          
          ### /!\ MODIF ICI POUR STEP 2 /!\ ###
          if (length(target.genes) >= min.positive){
            # if all conditions are met, list all target genes with
            # their correlations to the receptor in a data frame
            # row. Target genes are sorted wrt correlations.
            c <- corrg[r, target.genes]
            o <- order(c)
            c <- c[o]
            target.genes <- target.genes[o]
            #cat("\n c")
            data.frame(pathway=p, target.corr=paste(c,collapse=";"),
                       target.genes=paste(target.genes,collapse=";"),
                       len=length(c),
                       stringsAsFactors=FALSE)
          }
          else
            NULL
        }
        else
          NULL
      }
      if (!is.null(best.2nd))
        # one or several pathways containing the receptor r were found,
        # combine them in |-separated strings
        data.frame(R=r, pathways=paste(best.2nd$pathway, collapse="|"),
                   target.corr=paste(best.2nd$target.corr, collapse='|'),
                   target.genes=paste(best.2nd$target.genes, collapse='|'),
                   len=paste(best.2nd$len, collapse='|'),
                   stringsAsFactors=FALSE)
      else
        NULL
    }
    else
      NULL
  }
  
  # combine LR pair correlations with R-target gene correlations
  rownames(reg.proc) <- reg.proc$R
  conf.pairs <- lr[lr$R %in% reg.proc$R,]
  conf.pairs$pwid <- reg.proc[conf.pairs$R, "pathways"]
  conf.pairs$target.corr <- reg.proc[conf.pairs$R, "target.corr"]
  conf.pairs$len <- reg.proc[conf.pairs$R, "len"]
  conf.pairs$target.genes <- reg.proc[conf.pairs$R, "target.genes"]
  pw.name <- unique(pw[,c(id.col, pw.col)])
  pw2name <- stats::setNames(pw.name[[2]], pw.name[[1]])
  conf.pairs$pwname <- foreach::foreach(pl=conf.pairs$pwid,.combine=c) %do% {
    paste(foreach::foreach(id=unlist(strsplit(pl,"\\|")),
                           .combine=c) %do% {
                             pw2name[id]
                           },
          collapse="|"
    )
  }
  cat("\n conf.pairs: \n")
  cat(unlist(head(conf.pairs[,c("L", "R", "corr", "pwid", "pwname", "len", "target.genes",
                          "target.corr")])))
  cat("\n .downstreamSignaling - end \n")
  conf.pairs[,c("L", "R", "corr", "pwid", "pwname", "len", "target.genes",
                "target.corr")]
  
}  # .downstreamSignaling

#' Internal function to check receptor signaling downstream
#'
#' @param lr              A data frame as returned by
#'   \code{.getCorrelatedRP()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway
#'   IDs).
#' @param rncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param pncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param id.col          Column index or name in \code{pw} for the pathway IDs.
#' @param gene.col        Column index or name in \code{pw} for the gene
#'   symbols.
#' @param pw.col          Column index or name in \code{pw} for the pathway
#'   names.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#'
#' @return A table reporting all the receptor-PTMprotein pairs provided in \code{lr}
#'   along with the pathways found and data about target gene correlations with
#'   the receptor.
#'
#' @importFrom foreach %do% %dopar%
#' @keywords internal
.downstreamSignalingPTM <- function(lr, pw, pw.size, rncounts, PTM, id.col, gene.col,
                                 pw.col, min.positive=1, with.complex = TRUE,symPos=NULL) {
  if (!is.matrix(rncounts))
    stop("rncounts must be a matrix")
  
  # local binding
  r <- p <- pl <- id <- NULL
  cat("\n .downstreamSignalingPTM - inferencePTM ")
  # define interaction types
  control.int <- "controls-expression-of"
  incomplex.int <- c("in-complex-with","interacts-with")
  directed.int <- c("controls-state-change-of", "catalysis-precedes",
                    "controls-expression-of", "controls-transport-of",
                    "controls-phosphorylation-of")
  PTMph.int <- "controls-phosphorylation-of"
  if (with.complex)
    correlated.int <- union(control.int, incomplex.int)
  else
    correlated.int <- control.int
  
  # compute downstream correlations
  cat("\n Avant filtre: ")
  cat("\n dim ncounts: ")
  cat(dim(rncounts))
  cat("\n dim PTM: ")
  cat(dim(PTM))
  
  #rncounts <- rncounts[rownames(rncounts) %in% rownames(PTM),]
  rncounts <- rncounts[,colnames(rncounts) %in% colnames(PTM)]
  rncounts <- rncounts[,order(colnames(rncounts))]
  #rncounts <- rncounts[order(rownames(rncounts)),]
  PTM <- PTM[,colnames(PTM) %in% colnames(rncounts)]
  PTM <- PTM[,order(colnames(PTM))]
  #PTM <- PTM[order(rownames(PTM)),]
  
  cat("\n Apres filtre: ")
  cat("\n dim ncounts: ")
  cat(dim(rncounts))
  cat("\n dim PTM: ")
  cat(dim(PTM))
  
  corrg <- stats::cor(t(rncounts), t(PTM), method = "spearman")
  cat("\n dim corrg: ")
  cat(dim(corrg))
  # the global computation above is faster than restricted to the receptors
  corrg <- corrg[unique(lr$R), ]
  cat("\n 1")
  cat("\n dim corrg bis: ")
  cat(dim(corrg))
  # check each putative LR pair, loop over the receptors
  reg.proc <- foreach::foreach(r=unique(lr$R),.combine=rbind) %do% {
    # reg.proc <- NULL
    # for (r in unique(lr$putative.pairs$R)){
    
    # loop over the pathways containing the receptor r
    pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
    if (length(pa)>0){
      #cat("\n 2")
      receptor.ligands <- unique(lr$L[lr$R==r])
      m<-0
      best.2nd <- foreach::foreach(p=pa,.combine=rbind) %do% {
        # best.2nd <- NULL
        # for (p in pa){
        m <- m+1
        int <- SingleCellSignalR::PwC_ReactomeKEGG[ #get toutes interactions
          SingleCellSignalR::PwC_ReactomeKEGG$a.gn %in% pw[pw[[id.col]]==p,gene.col] &
            SingleCellSignalR::PwC_ReactomeKEGG$b.gn %in% pw[pw[[id.col]]==p,gene.col],
        ]
        directed <- int$type %in% directed.int
        
        # double the undirected interactions and generate a directed graph
        ret <- int[!directed,c("a.gn", "b.gn")]
        from <- ret$a.gn
        ret$a.gn <- ret$b.gn
        ret$b.gn <- from
        d.int <- unique(rbind(int[,c("a.gn", "b.gn")],ret))
        g <- igraph::graph_from_data_frame(d.int, directed=TRUE)
        
        # extract the target genes of receptor r
        if (r %in% d.int$a.gn || r %in% d.int$b.gn){
          # putative targets in the pathway
          target.genes <- setdiff(c( #pas forcement dans ncounts
            int[int$type %in% correlated.int & int$a.gn==r, "b.gn"],
            int[int$type %in% correlated.int & int$b.gn==r, "a.gn"],
            int[int$type %in% directed.int, "b.gn"]),
            r
          )
          #cat("\n 3")
          # reduce putative to reachable from the receptor
          sp <- igraph::shortest.paths(g, r, target.genes)
          target.genes <- colnames(sp)[!is.infinite(sp[r,])]
          if(m==1 | m%%20==0){
            cat("\n dim target.genes: ")
            cat(length(target.genes))
          }
          # eliminate ligands of the receptor if present
          target.genes <- setdiff(target.genes, receptor.ligands)
          
          #get target genes that can be PTMrylated 
          #target.genes.bis <- mergedPwCtrl["PathwayName"==p & "source" %in% c(target.genes,receptor.ligands) & "target" %in% target.genes, "genePos"]
          target.genes.bis <- subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands) & target %in% target.genes)$genePos
          #target.genes.bis.name <- mergedPwCtrl["PathwayName"==p & "source" %in% target.genes, "target"]
          #target.genes.bis.pos <- mergedPwCtrl["PathwayName"==p & "source" %in% target.genes, "pos"]
          
          target.genes <- target.genes.bis
          if(m==1 | m%%20==0){
            cat("\n dim target.genes bis: ")
            cat(length(target.genes))
          }
          if (length(target.genes) >= min.positive && any(target.genes %in% colnames(corrg))){
            # if all conditions are met, list all target genes with
            # their correlations to the receptor in a data frame
            # row. Target genes are sorted wrt correlations.
            cat("\n length(target.genes) >= min.positive \n")
            cat(target.genes)
            cat("\n corrg \n")
            cat(corrg[1:2,1:2])
            cat("\n r corrg : ")
            cat(r %in% rownames(corrg))
            cat("\n tg corrg : ")
            cat(target.genes %in% colnames(corrg))
            cat("\n rn corrg : ")
            cat(rownames(corrg)[1:2])
            cat("\n cn corrg : ")
            cat(colnames(corrg)[1:2])
            if(symPos==NULL)
              target.genes <- target.genes[target.genes %in% c(colnames(corrg),symPos[,1])]
            else
              target.genes <- target.genes[target.genes %in% c(colnames(corrg),rownames(corrg))]
            if (length(target.genes) > 0 && r %in% rownames(corrg)) {
              c <- corrg[r, target.genes]
              cat("\n c \n")
              cat(unlist(c))
              o <- order(c)
              c <- c[o]
            } else {
              c <- NULL
            }
            target.genes <- target.genes[o]
            cat("\n head: \n")
            cat(unlist(head(data.frame(pathway=p, target.corr=paste(c,collapse=";"),
                           target.genes=paste(target.genes,collapse=";"),
                           len=length(c),
                           stringsAsFactors=FALSE))))
            data.frame(pathway=p, target.corr=paste(c,collapse=";"),
                       target.genes=paste(target.genes,collapse=";"),
                       len=length(c),
                       stringsAsFactors=FALSE)
          }
          else{
            cat("\n no a11")
            data.frame(pathway=character(0), 
                       target.corr=character(0),
                       target.genes=character(0), 
                       len=integer(0),
                       stringsAsFactors=FALSE)
          }
        }
        else{
          cat("\n no a22")
          data.frame(pathway=character(0), 
                     target.corr=character(0),
                     target.genes=character(0), 
                     len=integer(0),
                     stringsAsFactors=FALSE)
        }
      }
      cat("\n before isnull best2 ")
      cat(is.null(best.2nd))
      #cat("\n")
      #cat(head(best.2nd))
      if (!is.null(best.2nd)){
        # one or several pathways containing the receptor r were found,
        # combine them in |-separated strings
        cat("\n in isnull best2")
        data.frame(R=r, pathways=paste(best.2nd$pathway, collapse="|"),
                   target.corr=paste(best.2nd$target.corr, collapse='|'),
                   target.genes=paste(best.2nd$target.genes, collapse='|'),
                   len=paste(best.2nd$len, collapse='|'),
                   stringsAsFactors=FALSE)
      }
      else{
        cat("\n no1")
        #NULL
        data.frame(pathway=character(0), 
                   target.corr=character(0),
                   target.genes=character(0), 
                   len=integer(0),
                   stringsAsFactors=FALSE)
      }
        
    }
    else{
      cat("\n no2")
      #NULL
      data.frame(pathway=character(0), 
                 target.corr=character(0),
                 target.genes=character(0), 
                 len=integer(0),
                 stringsAsFactors=FALSE)
    }
  }
  
  # combine LR pair correlations with R-target gene correlations
  cat("\n 5")
  cat("\n cn(reg): ")
  cat(colnames(reg.proc))
  #cat("\n reg: \n")
  #cat(unlist(head(reg.proc)))
  rownames(reg.proc) <- reg.proc$R
  conf.pairs <- lr[lr$R %in% reg.proc$R,]
  conf.pairs$pwid <- reg.proc[conf.pairs$R, "pathways"]
  conf.pairs$target.corr <- reg.proc[conf.pairs$R, "target.corr"]
  conf.pairs$len <- reg.proc[conf.pairs$R, "len"]
  conf.pairs$target.genes <- reg.proc[conf.pairs$R, "target.genes"]
  pw.name <- unique(pw[,c(id.col, pw.col)])
  pw2name <- stats::setNames(pw.name[[2]], pw.name[[1]])
  cat("\n 6 \n")
  cat(unlist(head(conf.pairs)))
  conf.pairs$pwname <- foreach::foreach(pl=conf.pairs$pwid,.combine=c) %do% {
    paste(foreach::foreach(id=unlist(strsplit(pl,"\\|")),
                           .combine=c) %do% {
                             pw2name[id]
                           },
          collapse="|"
    )
  }
  
  conf.pairs[,c("L", "R", "corr", "pwid", "pwname", "len", "target.genes",
                "target.corr")]
  
}  # .downstreamSignalingPTM



#' Internal function to check receptor signaling downstream
#'
#' @param lr              A data frame as returned by
#'   \code{.getCorrelatedRP()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway
#'   IDs).
#' @param rncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param pncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param id.col          Column index or name in \code{pw} for the pathway IDs.
#' @param gene.col        Column index or name in \code{pw} for the gene
#'   symbols.
#' @param pw.col          Column index or name in \code{pw} for the pathway
#'   names.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#'
#' @return A table reporting all the receptor-PTMprotein pairs provided in \code{lr}
#'   along with the pathways found and data about target gene correlations with
#'   the receptor.
#'
#' @importFrom foreach %do% %dopar%
#' @keywords internal
.downstreamSignalingPTMInf <- function(lr, pw, pw.size, rncounts, PTM, id.col, gene.col,
                                        pw.col, min.positive=1, with.complex = TRUE,symPos=NULL) {
  if (!is.matrix(rncounts))
    stop("rncounts must be a matrix")
  fichier_sortie <- "/data2/USERS/giroux/PTM/zhang/bisNA/sortie.txt"
  # local binding
  r <- p <- pl <- id <- NULL
  cat("\n .downstreamSignalingPTM - inferencePTM ")
  # define interaction types
  control.int <- "controls-expression-of"
  incomplex.int <- c("in-complex-with","interacts-with")
  directed.int <- c("controls-state-change-of", "catalysis-precedes",
                    "controls-expression-of", "controls-transport-of",
                    "controls-phosphorylation-of")
  PTMph.int <- "controls-phosphorylation-of"
  if (with.complex)
    correlated.int <- union(control.int, incomplex.int)
  else
    correlated.int <- control.int
  
  # compute downstream correlations
  cat("\n Avant filtre: ")
  cat("\n dim ncounts: ")
  cat(dim(rncounts))
  cat("\n dim PTM: ")
  cat(dim(PTM))
  
  #rncounts <- rncounts[rownames(rncounts) %in% rownames(PTM),]
  rncounts <- rncounts[,colnames(rncounts) %in% colnames(PTM)]
  rncounts <- rncounts[,order(colnames(rncounts))]
  #rncounts <- rncounts[order(rownames(rncounts)),]
  PTM <- PTM[,colnames(PTM) %in% colnames(rncounts)]
  PTM <- PTM[,order(colnames(PTM))]
  #PTM <- PTM[order(rownames(PTM)),]
  
  cat("\n Apres filtre: ")
  cat("\n dim ncounts: ")
  cat(dim(rncounts))
  cat("\n dim PTM: ")
  cat(dim(PTM))
  
  
  corrg <- stats::cor(t(rncounts), method = "spearman")
  corrgp <- stats::cor(t(rncounts), t(PTM), method = "spearman")
  cat("\n dim corrg: ")
  cat(dim(corrg))
  # the global computation above is faster than restricted to the receptors
  corrg <- corrg[unique(lr$R), ]
  cat("\n 1")
  cat("\n dim corrg bis: ")
  cat(dim(corrg))
  # check each putative LR pair, loop over the receptors
  reg.proc <- foreach::foreach(r=unique(lr$R),.combine=rbind) %do% {
    # reg.proc <- NULL
    # for (r in unique(lr$putative.pairs$R)){
    
    # loop over the pathways containing the receptor r
    pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
    if (length(pa)>0){
      #cat("\n 2")
      receptor.ligands <- unique(lr$L[lr$R==r])
      m<-0
      best.2nd <- foreach::foreach(p=pa,.combine=rbind) %do% {
        # best.2nd <- NULL
        # for (p in pa){
        m <- m+1
        int <- SingleCellSignalR::PwC_ReactomeKEGG[ #get toutes interactions
          SingleCellSignalR::PwC_ReactomeKEGG$a.gn %in% pw[pw[[id.col]]==p,gene.col] &
            SingleCellSignalR::PwC_ReactomeKEGG$b.gn %in% pw[pw[[id.col]]==p,gene.col],
        ]
        directed <- int$type %in% directed.int
        
        # double the undirected interactions and generate a directed graph
        ret <- int[!directed,c("a.gn", "b.gn")]
        from <- ret$a.gn
        ret$a.gn <- ret$b.gn
        ret$b.gn <- from
        d.int <- unique(rbind(int[,c("a.gn", "b.gn")],ret))
        g <- igraph::graph_from_data_frame(d.int, directed=TRUE)
        
        # extract the target genes of receptor r
        if (r %in% d.int$a.gn || r %in% d.int$b.gn){
          # putative targets in the pathway
          target.genes <- setdiff(c( #pas forcement dans ncounts
            int[int$type %in% correlated.int & int$a.gn==r, "b.gn"],
            int[int$type %in% correlated.int & int$b.gn==r, "a.gn"],
            int[int$type %in% directed.int, "b.gn"]),
            r
          )
          #cat("\n 3")
          # reduce putative to reachable from the receptor
          sp <- igraph::shortest.paths(g, r, target.genes)
          target.genes <- colnames(sp)[!is.infinite(sp[r,])]
          # if(m==1 | m%%20==0){
          #   cat("\n dim target.genes: ")
          #   cat(length(target.genes))
          # }
          # eliminate ligands of the receptor if present
          target.genes <- setdiff(target.genes, receptor.ligands)
          
          #get target genes that can be PTMrylated 
          #i.e. source and target detected in tg
          
          #target.genes.bis <- mergedPwCtrl["PathwayName"==p & "source" %in% c(target.genes,receptor.ligands) & "target" %in% target.genes, "genePos"]
          target.genes.PTM.pos <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands) & target %in% target.genes)$genePos)
          target.genes.PTM.name <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands) & target %in% target.genes)$target)
          #target.genes.bis.pos <- mergedPwCtrl["PathwayName"==p & "source" %in% target.genes, "pos"]
          if(p=="R-HSA-9607240" | p=="R-HSA-6811558"){
            cat(p, file = fichier_sortie, append = TRUE)
            cat("\n target.genes.PTM.pos : ", file = fichier_sortie, append = TRUE)
            cat(target.genes.PTM.pos, file = fichier_sortie, append = TRUE)
          }
          #target.genes <- target.genes.bis
          # if(m==1 | m%%20==0){
          #   cat("\n dim target.genes bis: ")
          #   cat(length(target.genes))
          # }
          if (length(target.genes) >= min.positive){
            # if all conditions are met, list all target genes with
            # their correlations to the receptor in a data frame
            # row. Target genes are sorted wrt correlations.
            #cat("\n length(target.genes) >= min.positive \n")
            # cat(target.genes)
            # cat("\n corrg \n")
            # cat(corrg[1:2,1:2])
            # cat("\n r corrg : ")
            # cat(r %in% rownames(corrg))
            # cat("\n tg corrg : ")
            # cat(target.genes %in% colnames(corrg))
            # cat("\n rn corrg : ")
            # cat(rownames(corrg)[1:2])
            # cat("\n cn corrg : ")
            # cat(colnames(corrg)[1:2])
            target.genes <- target.genes[target.genes %in% c(colnames(corrg),symPos[,1])]
            c <- corrg[r, target.genes]
            # cat("\n c \n")
            # cat(head(unlist(c)))
            # cat("\n cdim \n")
            # cat(length(c))
            
            
            PTMTheoricalPres <- unique(target.genes.PTM.pos[target.genes.PTM.pos %in% colnames(corrgp)])
            tmpNA <- grepl("_NA", target.genes.PTM.pos, fixed = TRUE)
            target.genes.PTM.posNA <- target.genes.PTM.pos[tmpNA]
            target.genes.PTM.nameNA <- target.genes.PTM.name[tmpNA]
            if(length(target.genes.PTM.posNA) > 0){
              listPos <- c()
              for(g in target.genes.PTM.nameNA){
                tmpCol <- grepl(paste0(g,"_"), colnames(corrgp))
                tmpPo <- NA
                tmpPo <- colnames(corrgp)[tmpCol]
                listPos <- c(listPos, tmpPo)
              }
              if(!is.null(tmpPo) && length(tmpPo) > 0)
                PTMTheoricalPres <- c(PTMTheoricalPres, tmpPo)
            }
            target.genes.PTM <- PTMTheoricalPres
            if(p=="R-HSA-9607240" | p=="R-HSA-6811558"){
              cat("\n PTMTheoricalPres : ", file = fichier_sortie, append = TRUE)
              cat(PTMTheoricalPres, file = fichier_sortie, append = TRUE)
            }
            if(length(PTMTheoricalPres)>0 && r %in% rownames(corrgp)){
              cp <- corrgp[r, PTMTheoricalPres]
              
              # cat("\n \n \n dim c: \n")
              # cat(length(c))
              # cat(", ")
              # cat(nrow(c))
              # cat("\n dim cp: \n")
              # cat(length(cp))
              # cat(", ")
              # cat(nrow(cp))
              # cat("\n")
              
              #c <- c(c, cp)
              
              # cat("\n dim c bis: \n")
              # cat(length(c))
              # cat("\n")
              # cat(c)
              
              op <- order(cp)
              cp <- cp[op]
              target.genes.PTM <- target.genes.PTM[op]
              len.PTM <- length(cp)
            }
            else{
              cp <- c(NA)
              len.PTM <- 0
            }
            if(p=="R-HSA-9607240" | p=="R-HSA-6811558"){
              cat("\n target.genes.PTM : ", file = fichier_sortie, append = TRUE)
              cat(target.genes.PTM, file = fichier_sortie, append = TRUE)
            }
            # cat("\n \n \n 2dim c: \n")
            # cat(length(c))
            # cat(", ")
            # cat(nrow(c))
            # cat("\n 2dim cp: \n")
            # cat(length(cp))
            # cat(", ")
            # cat(nrow(cp))
            # cat("\n")
            
            o <- order(c)
            c <- c[o]
            target.genes <- target.genes[o]
            # cat("\n head: \n")
            # cat(unlist(head(data.frame(pathway=p, target.corr=paste(c,collapse=";"), target.PTM.corr=paste(cp,collapse=";"),
            #                            target.genes=paste(target.genes,collapse=";"),
            #                            target.genes.PTM=paste(target.genes.PTM,collapse=";"),
            #                            len=length(c),
            #                            len.PTM=len.PTM,
            #                            stringsAsFactors=FALSE))))
            data.frame(pathway=p, target.corr=paste(c,collapse=";"), target.PTM.corr=paste(cp,collapse=";"),
                       target.genes=paste(target.genes,collapse=";"),
                       target.genes.PTM=paste(target.genes.PTM,collapse=";"),
                       len=length(c),
                       len.PTM=len.PTM,
                       stringsAsFactors=FALSE)
          }
          else{
            #cat("\n no a1")
            NULL
          }
        }
        else{
          #cat("\n no a2")
          NULL
        }
      }
      #cat("\n before isnull best2 ")
      #cat(is.null(best.2nd))
      #cat("\n")
      #cat(head(best.2nd))
      if (!is.null(best.2nd)){
        # one or several pathways containing the receptor r were found,
        # combine them in |-separated strings
        #cat("\n in isnull best2")
        data.frame(R=r, pathways=paste(best.2nd$pathway, collapse="|"),
                   target.corr=paste(best.2nd$target.corr, collapse='|'),
                   target.genes=paste(best.2nd$target.genes, collapse='|'),
                   target.PTM.corr=paste(best.2nd$target.PTM.corr, collapse='|'),
                   target.genes.PTM=paste(best.2nd$target.genes.PTM, collapse='|'),
                   len=paste(best.2nd$len, collapse='|'),
                   len.PTM=paste(best.2nd$len.PTM, collapse='|'),
                   stringsAsFactors=FALSE)
      }
      else{
        #cat("\n no1")
        NULL
      }
      
    }
    else{
      #cat("\n no2")
      NULL
    }
  }
  
  # combine LR pair correlations with R-target gene correlations
  #cat("\n 5")
  #cat("\n cn(reg): ")
  #cat(colnames(reg.proc))
  #cat("\n reg: \n")
  #cat(unlist(head(reg.proc)))
  rownames(reg.proc) <- reg.proc$R
  conf.pairs <- lr[lr$R %in% reg.proc$R,]
  conf.pairs$pwid <- reg.proc[conf.pairs$R, "pathways"]
  conf.pairs$target.corr <- reg.proc[conf.pairs$R, "target.corr"]
  conf.pairs$len <- reg.proc[conf.pairs$R, "len"]
  conf.pairs$target.genes <- reg.proc[conf.pairs$R, "target.genes"]
  conf.pairs$target.PTM.corr <- reg.proc[conf.pairs$R, "target.PTM.corr"]
  conf.pairs$len.PTM <- reg.proc[conf.pairs$R, "len.PTM"]
  conf.pairs$target.genes.PTM <- reg.proc[conf.pairs$R, "target.genes.PTM"]
  pw.name <- unique(pw[,c(id.col, pw.col)])
  pw2name <- stats::setNames(pw.name[[2]], pw.name[[1]])
  #cat("\n 6 \n")
  #cat(unlist(head(conf.pairs)))
  conf.pairs$pwname <- foreach::foreach(pl=conf.pairs$pwid,.combine=c) %do% {
    paste(foreach::foreach(id=unlist(strsplit(pl,"\\|")),
                           .combine=c) %do% {
                             pw2name[id]
                           },
          collapse="|"
    )
  }
  
  conf.pairs[,c("L", "R", "corr", "pwid", "pwname", "len", "target.genes",
                "target.corr", "len.PTM", "target.genes.PTM",
                "target.PTM.corr")]
  
}  # .downstreamSignalingPTMInf


#' Internal function to check receptor signaling downstream
#'
#' Assess the existence of correlations between a receptor,
#' part of a ligand-receptor pair, and
#' genes coding for proteins forming a complex with the receptor or genes
#' regulated by the receptor downstream signaling.
#'
#' @param ds              A BSRDataModel object.
#' @param lr              A table as returned by \code{.getCorrelatedLR()}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param infPTM    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @return A data frame extending \code{lr} content with the pathways found to
#' contain the receptors and data about target gene correlations with those
#' receptors. Strings in semi-colon-separated format are used to report
#' target genes and their Spearman correlations with the receptor in the
#' data frame. The target genes are sorted according to the correlation
#' coefficient.
#'
#' In a pathway of the reference, i.e., a Reactome pathway or the genes of a
#' GOBP term, the target genes are the
#' genes coding for proteins forming a complex with the receptor and the
#' genes in the pathway downstream the receptor,
#' which are given as regulated by the pathway. If \code{with.complex} is
#' set to \code{FALSE}, then only the
#' regulated genes are considered. Participation to a complex and being
#' regulated as well as the pathway directed topologies
#' are defined by Reactome and KEGG pathways as provided by PathwayCommons.
#'
#' The maximum pathway size is used to limit the redundancy inherent to GOBP
#' and Reactome. The minimum pathway size is
#' used to avoid overspecific, noninformative results.
#'
#' @importFrom methods is
#' @keywords internal
.checkReceptorSignaling <- function(ds, lr, reference=c("REACTOME"),
                                    max.pw.size=200, min.pw.size=5,
                                    min.positive=4, restrict.pw=NULL,
                                    with.complex=TRUE, infPTM=FALSE){
  
  # if (!is(ds, "BSRDataModel"))
  #     stop("ds must be a BSRDataModel object")
  
  ## dans function(... !!) ###
  # reference=c("REACTOME-GOBP",
  #             "REACTOME","GOBP")
  
  reference <- match.arg(reference)
  results <- list()
  cat("\n checkReceptorSignaling - inferencePTM ")
  
  # Reactome pathways
  if (reference %in% c("REACTOME-GOBP","REACTOME")){
    react <- reactome[reactome$`Gene name` %in% rownames(ncounts(ds)),]
    cat("\n react ok \n")
    if (!is.null(restrict.pw))
      react <- react[react$`Reactome ID` %in% restrict.pw,]
    pw.size <- table(react$`Reactome ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- react[react$`Gene name` %in% lr$R, "Reactome ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         react[react$`Reactome ID` %in% names(pw.size), "Gene name"])
    )
    corgenes <- unique(corgenes[corgenes %in% rownames(ncounts(ds))])
    cat("\n corgenes: \n", corgenes)
    cat("\n \n rn: \n", rownames(ncounts(ds)))
    cat("\n \n all: \n", length(corgenes))
    cat("\n \n in: \n", sum(corgenes %in% rownames(ncounts(ds))))
    #cat("\n \n alone: \n", corgenes[!c(corgenes %in% rownames(ncounts(ds)))])
    if(infPTM){
      cat("\n infPTM 1r")
      if(is.null(ds@symPos) || is.na(ds@symPos) || dim(ds@symPos) == c(1,1)){
        cat("\n if \n")
        corgenesP <- rownames(ds@PTM)[rownames(ds@PTM) %in% corgenes]
        cat("\n corgenesP: \n", corgenesP)
      }
      else{
        cat("\n else \n")
        corgenesposPbool <- ds@symPos[,1] %in% corgenes
        corgenesnamesP <- ds@symPos[corgenesposPbool,1]
        corgenesposP <- ds@symPos[corgenesposPbool,2]
        corgenesP <- paste0(corgenesnamesP, "_", corgenesposP)#AAS_12
        cat("\n corgenesP: \n", corgenesP)
      }
      cat("\n infPTM 12r")
      PTMCorr <- ds@PTM[corgenesP,]
      corgenes2 <- corgenes[corgenes %in% rownames(ncounts(ds))]
      cat("\n infPTM 2r")
      results$reactome.pairs <- .downstreamSignalingPTM(lr, react, pw.size,
                                                            ncounts(ds)[corgenes2,], PTMCorr, id.col="Reactome ID", gene.col="Gene name",
                                                            pw.col="Reactome name", min.positive, with.complex=with.complex,symPos=ds@symPos)
    }
    else{
      cat("\n notInfPTM 1r")
      results$reactome.pairs <- .downstreamSignaling(lr, react, pw.size,
                                                     ncounts(ds)[corgenes,], id.col="Reactome ID", gene.col="Gene name",
                                                     pw.col="Reactome name", min.positive, with.complex=with.complex)
    }
    
  }
  
  # GOBP
  if (reference %in% c("REACTOME-GOBP","GOBP")){
    go <- gobp[gobp$`Gene name` %in% rownames(ncounts(ds)),]
    if (!is.null(restrict.pw))
      go <- go[go$`GO ID` %in% restrict.pw,]
    pw.size <- table(go$`GO ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- go[go$`Gene name` %in% lr$R, "GO ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         go[go$`GO ID` %in% names(pw.size), "Gene name"])
    )
    if(infPTM){
      #cat("\n infPTM 1g")
      corgenesposPbool <- ds@symPos[,1] %in% corgenes
      corgenesnamesP <- ds@symPos[corgenesposPbool,1]
      corgenesposP <- ds@symPos[corgenesposPbool,2]
      corgenesP <- paste0(corgenesnamesP, "_", corgenesposP)#AAS_12
      PTMCorr <- ds@PTM[corgenesP,]
      corgenes2 <- corgenes[corgenes %in% rownames(ncounts(ds))]
      #cat("\n infPTM 2g")
      results$gobp.pairs <- .downstreamSignalingPTM(lr, go, pw.size,
                                                        ncounts(ds)[corgenes2,], PTMCorr, id.col="GO ID", gene.col="Gene name",
                                                        pw.col="GO name", min.positive, with.complex=with.complex,symPos=ds@symPos)
    }
    else{
      #cat("\n notInfPTM 1g")
      results$gobp.pairs <- .downstreamSignaling(lr, go, pw.size,
                                                 ncounts(ds)[corgenes,], id.col="GO ID", gene.col="Gene name",
                                                 pw.col="GO name", min.positive, with.complex=with.complex)
    }
    cat("\n ref ok")
  }
  
  # merge
  if (reference == "REACTOME-GOBP"){
    #cat("\n mergeRG")
    pairs <- unique(rbind(results$reactome.pairs[,1:2],
                          results$gobp.pairs[,1:2])
    )
    react.keys <- paste(results$reactome.pairs[[1]],
                        results$reactome.pairs[[2]], sep="|")
    gobp.keys <- paste(results$gobp.pairs[[1]],
                       results$gobp.pairs[[2]], sep="|")
    results$merged.pairs <- rbind(results$reactome.pairs,
                                  results$gobp.pairs[!(gobp.keys %in% react.keys),]
    )
  }
  else if (reference == "REACTOME"){
    cat("\n mergeR")
    results$merged.pairs <- results$reactome.pairs
  }
  else{
    #cat("\n mergeG")
    results$merged.pairs <- results$gobp.pairs
  }
  results$merged.pairs
  
} # .checkReceptorSignaling



#' Internal function to check receptor signaling downstream
#'
#' Assess the existence of correlations between a receptor,
#' part of a ligand-receptor pair, and
#' genes coding for proteins forming a complex with the receptor or genes
#' regulated by the receptor downstream signaling.
#'
#' @param ds              A BSRDataModel object.
#' @param lr              A table as returned by \code{.getCorrelatedLR()}.
#' @param reference       Which pathway reference should be used ("REACTOME"
#'   for Reactome, "GOBP" for GO Biological Process,
#'   or "REACTOME-GOBP" for both).
#' @param max.pw.size     Maximum pathway size to consider from the pathway
#'   reference.
#' @param min.pw.size     Minimum pathway size to consider from the pathway
#'   reference.
#' @param min.positive    Minimum number of target genes to be found in a given
#'   pathway.
#' @param restrict.pw     A list of pathway IDs to restrict the application of
#'   the function.
#' @param with.complex    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @param infPTM    A logical indicating whether receptor co-complex
#'   members should be included in the target genes.
#' @return A data frame extending \code{lr} content with the pathways found to
#' contain the receptors and data about target gene correlations with those
#' receptors. Strings in semi-colon-separated format are used to report
#' target genes and their Spearman correlations with the receptor in the
#' data frame. The target genes are sorted according to the correlation
#' coefficient.
#'
#' In a pathway of the reference, i.e., a Reactome pathway or the genes of a
#' GOBP term, the target genes are the
#' genes coding for proteins forming a complex with the receptor and the
#' genes in the pathway downstream the receptor,
#' which are given as regulated by the pathway. If \code{with.complex} is
#' set to \code{FALSE}, then only the
#' regulated genes are considered. Participation to a complex and being
#' regulated as well as the pathway directed topologies
#' are defined by Reactome and KEGG pathways as provided by PathwayCommons.
#'
#' The maximum pathway size is used to limit the redundancy inherent to GOBP
#' and Reactome. The minimum pathway size is
#' used to avoid overspecific, noninformative results.
#'
#' @importFrom methods is
#' @keywords internal
.checkReceptorSignalingInf <- function(ds, lr, reference=c("REACTOME"),
                                    max.pw.size=200, min.pw.size=5,
                                    min.positive=4, restrict.pw=NULL,
                                    with.complex=TRUE, infPTM=FALSE){
  
  # if (!is(ds, "BSRDataModel"))
  #     stop("ds must be a BSRDataModel object")
  
  ## dans function(... !!) ###
  # reference=c("REACTOME-GOBP",
  #             "REACTOME","GOBP")
  
  reference <- match.arg(reference)
  results <- list()
  cat("\n checkReceptorSignalingInf - inferencePTM ")
  
  # Reactome pathways
  if (reference %in% c("REACTOME-GOBP","REACTOME")){
    react <- reactome[reactome$`Gene name` %in% rownames(ncounts(ds)),]
    if (!is.null(restrict.pw))
      react <- react[react$`Reactome ID` %in% restrict.pw,]
    pw.size <- table(react$`Reactome ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- react[react$`Gene name` %in% lr$R, "Reactome ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         react[react$`Reactome ID` %in% names(pw.size), "Gene name"])
    )
    if(infPTM){
      #cat("\n infPTM 1r")
      # corgenesposPbool <- ds@symPos[,1] %in% corgenes
      # corgenesnamesP <- ds@symPos[corgenesposPbool,1]
      # corgenesposP <- ds@symPos[corgenesposPbool,2]
      # corgenesP <- paste0(corgenesnamesP, "_", corgenesposP)#AAS_12
      # PTMCorr <- ds@PTM[corgenesP,]
      # corgenes2 <- corgenes[corgenes %in% rownames(ncounts(ds))]
      # cat("\n infPTM 2r")
      # results$reactome.pairs <- .downstreamSignalingPTMInf(lr, react, pw.size,
      #                                                       ncounts(ds)[corgenes2,], PTMCorr, id.col="Reactome ID", gene.col="Gene name",
      #                                                       pw.col="Reactome name", min.positive, with.complex=with.complex,symPos=ds@symPos)
      results$reactome.pairs <- .downstreamSignalingPTMInf(lr, react, pw.size,
                                                               ncounts(ds), ds@PTM, id.col="Reactome ID", gene.col="Gene name",
                                                               pw.col="Reactome name", min.positive, with.complex=with.complex,symPos=ds@symPos)
    }
    else{
      #cat("\n notInfPTM 1r")
      results$reactome.pairs <- .downstreamSignaling(lr, react, pw.size,
                                                     ncounts(ds)[corgenes,], id.col="Reactome ID", gene.col="Gene name",
                                                     pw.col="Reactome name", min.positive, with.complex=with.complex)
    }
    
  }
  
  # GOBP
  if (reference %in% c("REACTOME-GOBP","GOBP")){
    go <- gobp[gobp$`Gene name` %in% rownames(ncounts(ds)),]
    if (!is.null(restrict.pw))
      go <- go[go$`GO ID` %in% restrict.pw,]
    pw.size <- table(go$`GO ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- go[go$`Gene name` %in% lr$R, "GO ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         go[go$`GO ID` %in% names(pw.size), "Gene name"])
    )
    if(infPTM){
      #cat("\n infPTM 1g")
      # corgenesposPbool <- ds@symPos[,1] %in% corgenes
      # corgenesnamesP <- ds@symPos[corgenesposPbool,1]
      # corgenesposP <- ds@symPos[corgenesposPbool,2]
      # corgenesP <- paste0(corgenesnamesP, "_", corgenesposP)#AAS_12
      # PTMCorr <- ds@PTM[corgenesP,]
      # corgenes2 <- corgenes[corgenes %in% rownames(ncounts(ds))]
      # cat("\n infPTM 2g")
      # results$gobp.pairs <- .downstreamSignalingPTM(lr, go, pw.size,
      #                                                   ncounts(ds)[corgenes2,], PTMCorr, id.col="GO ID", gene.col="Gene name",
      #                                                   pw.col="GO name", min.positive, with.complex=with.complex,symPos=ds@symPos)
      results$reactome.pairs <- .downstreamSignalingPTMInf(lr, go, pw.size,
                                                               ncounts(ds), ds@PTM, id.col="GO ID", gene.col="Gene name",
                                                               pw.col="GO name", min.positive, with.complex=with.complex,symPos=ds@symPos)
    }
    else{
      #cat("\n notInfPTM 1g")
      results$gobp.pairs <- .downstreamSignaling(lr, go, pw.size,
                                                 ncounts(ds)[corgenes,], id.col="GO ID", gene.col="Gene name",
                                                 pw.col="GO name", min.positive, with.complex=with.complex)
    }
    #cat("\n ref ok")
  }
  
  # merge
  if (reference == "REACTOME-GOBP"){
    #cat("\n mergeRG")
    pairs <- unique(rbind(results$reactome.pairs[,1:2],
                          results$gobp.pairs[,1:2])
    )
    react.keys <- paste(results$reactome.pairs[[1]],
                        results$reactome.pairs[[2]], sep="|")
    gobp.keys <- paste(results$gobp.pairs[[1]],
                       results$gobp.pairs[[2]], sep="|")
    results$merged.pairs <- rbind(results$reactome.pairs,
                                  results$gobp.pairs[!(gobp.keys %in% react.keys),]
    )
  }
  else if (reference == "REACTOME"){
    #cat("\n mergeR")
    results$merged.pairs <- results$reactome.pairs
  }
  else{
    #cat("\n mergeG")
    results$merged.pairs <- results$gobp.pairs
  }
  results$merged.pairs
  
} # .checkReceptorSignalingInf


#' Internal function to assign P-values to LR interactions
#'
#' Estimate the P-value of each ligand-receptor pair based
#' on the data frame output by \code{\link{.checkReceptorSignaling}}.
#'
#' @param pairs         A data frame output by \code{checkReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @return A data.frame with the data in \code{pairs} complemented with
#' P-values and adjusted P-values.
#' @keywords internal
.pValuesLR <- function(pairs, param, rank.p = 0.75,
                       fdr.proc = c("BH", "Bonferroni", "Holm", "Hochberg",
                                    "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {
  
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  
  #cat("\n pVal \n")
  # prepare the chosen model CDF
  LR.par <- param$LR.0$model
  RT.par <- param$RT.0$model
  RP.par <- param$RP.0$model
  #RP.par <- param$RT.0$model
  if (LR.par$distrib != RT.par$distrib)
    stop("Distinct statistical models for LR and RT nulls are not allowed")
  if (LR.par$distrib != RP.par$distrib)
    stop("Distinct statistical models for LR and RP nulls are not allowed")
  if (LR.par$distrib == 'censored_normal')
    cdf <- .cdfGaussian
  else if (LR.par$distrib == 'censored_mixed_normal')
    cdf <- .cdfMixedGaussian
  else if (LR.par$distrib == 'censored_stable')
    cdf <- .cdfAlphaStable
  else if (LR.par$distrib == 'empirical')
    cdf <- .cdfEmpirical
  else if (LR.par$distrib == 'kernel_empirical')
    cdf <- .cdfKernelEmpirical
  else
    stop(paste0("Unknown statistical model: ", LR.par$LR.0$model$distrib))
  
  # estimate P-values
  res <- NULL
  for (i in 1:nrow(pairs)){
    # all the data related to each pathway containing a given
    # receptor were collapsed separated by |
    # we need to split those pathways
    # cat("\n pVal \n")
    # cat(head(unlist(pairs)))
    # cat("\n")
    # cat(colnames(pairs))
    # cat("\n")
    # cat(pairs$target.genes.PTM)
    pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
    pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
    tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
    #pg <- unlist(strsplit(pairs$target.genes.PTM[i],split="\\|"))
    spear <- unlist(strsplit(pairs$target.corr[i],split="\\|"))
    len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))
    spear.PTM <- NA
    spear.PTM <- unlist(strsplit(pairs$target.PTM.corr[i],split="\\|"))
    len.PTM <- as.numeric(unlist(strsplit(pairs$len.PTM[i],split="\\|")))
    
    # estimate the LR correlation P-value
    if (pairs$corr[i] >= 0)
      # normal case
      p.lr <- 1 - cdf(pairs$corr[i], LR.par)
    else
      # to enable searching for inhibitory L-R interactions
      p.lr <- cdf(pairs$corr[i], LR.par)
    
    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    for (k in 1:length(len)){
      spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
      r <- min(max(1,trunc(rank.p*len[k])),len[k])
      rank.corr <- spears[r]
      # r-1 correlations are < rank.corr, prob to have r-1 or less
      # corr < rank.corr is given by a binomial with success rate
      # equal to the probability to get a corr < rank.corr, i.e.,
      # cdf(rank.corr, RT.par). If rank.corr is high, it becomes
      # difficult to get as little as r-1 corr < rank.corr by chance!
      p.rt <- stats::pbinom(r-1, len[k], cdf(rank.corr, RT.par))
      #p.rp <- stats::pbinom(r-1, len[k], cdf(rank.corr, RP.par))
      #create new df melant prot et PTM
      #learnparam(df)
      #p.rtp
      p.rtp <- p.rt
      #p.rtp[pg] <- p.rp[pg]
      
      res <- rbind(res,data.frame(pairs[i,c("L","R")],
                                  LR.corr=pairs[i,"corr"], pw.id=pwid[k],
                                  pw.name=pwname[k], rank=r, len=len[k], len.PTM=len.PTM[k],
                                  rank.corr=rank.corr, target.genes=tg[k],
                                  target.corr=spear[k], target.genes.PTM=pg[k],
                                  target.PTM.corr=spear.PTM[k], pval=p.lr*p.rtp,
                                  stringsAsFactors=FALSE))
    }
  }
  
  # avoid the impossible
  key <- paste(res$L, res$R, res$pw.id, sep="||")
  bad <- duplicated(key)
  res <- res[!bad,]
  
  # multiple hypothesis correction
  rawp <- res$pval
  adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
  res$qval <- adj$adjp[order(adj$index),fdr.proc]
  
  res
  
}  # .pValuesLR



#' Internal function to assign P-values to LRP interactions
#'
#' Estimate the P-value of each ligand-receptor pair based
#' on the data frame output by \code{\link{.checkReceptorSignaling}}.
#'
#' @param pairs         A data frame output by \code{checkReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @return A data.frame with the data in \code{pairs} complemented with
#' P-values and adjusted P-values.
#' @keywords internal
.pValuesLRP <- function(pairs, param, rank.p = 0.75,
                       fdr.proc = c("BH", "Bonferroni", "Holm", "Hochberg",
                                    "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {
  
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  
  cat("\n pVal LRP\n")
  # prepare the chosen model CDF
  LR.par <- param$LR.0$model
  RT.par <- param$RT.0$model
  RP.par <- param$RP.0$model
  #RP.par <- param$RT.0$model
  if (LR.par$distrib != RT.par$distrib)
    stop("Distinct statistical models for LR and RT nulls are not allowed")
  # if (LR.par$distrib != RP.par$distrib)
  #   stop("Distinct statistical models for LR and RP nulls are not allowed")
  if (LR.par$distrib == 'censored_normal')
    cdf <- .cdfGaussian
  else if (LR.par$distrib == 'censored_mixed_normal')
    cdf <- .cdfMixedGaussian
  else if (LR.par$distrib == 'censored_stable')
    cdf <- .cdfAlphaStable
  else if (LR.par$distrib == 'empirical')
    cdf <- .cdfEmpirical
  else if (LR.par$distrib == 'kernel_empirical')
    cdf <- .cdfKernelEmpirical
  else
    stop(paste0("Unknown statistical model: ", LR.par$LR.0$model$distrib))
  
  #rp
  if (RP.par$distrib == 'censored_normal')
    cdfp <- .cdfGaussian
  else if (RP.par$distrib == 'censored_mixed_normal')
    cdfp <- .cdfMixedGaussian
  else if (RP.par$distrib == 'censored_stable')
    cdfp <- .cdfAlphaStable
  else if (RP.par$distrib == 'empirical')
    cdfp <- .cdfEmpirical
  else if (RP.par$distrib == 'kernel_empirical')
    cdfp <- .cdfKernelEmpirical
  else
    stop(paste0("Unknown statistical model: ", RP.par$RP.0$model$distrib))
  
  # estimate P-values
  res <- NULL
  for (i in 1:nrow(pairs)){
    # all the data related to each pathway containing a given
    # receptor were collapsed separated by |
    # we need to split those pathways
    # cat("\n pVal \n")
    # cat(head(unlist(pairs)))
    # cat("\n")
    # cat(colnames(pairs))
    # cat("\n")
    # cat(pairs$target.genes.PTM)
    pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
    pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
    tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
    pg <- unlist(strsplit(pairs$target.genes.PTM[i],split="\\|"))
    spear <- unlist(strsplit(pairs$target.corr[i],split="\\|"))
    len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))
    spear.PTM <- NA
    spear.PTM <- unlist(strsplit(pairs$target.PTM.corr[i],split="\\|"))
    len.PTM <- as.numeric(unlist(strsplit(pairs$len.PTM[i],split="\\|")))
    
    # estimate the LR correlation P-value
    if (pairs$corr[i] >= 0)
      # normal case
      p.lr <- 1 - cdf(pairs$corr[i], LR.par)
    else
      # to enable searching for inhibitory L-R interactions
      p.lr <- cdf(pairs$corr[i], LR.par)
    
    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    for (k in 1:length(len)){
      spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
      r <- min(max(1,trunc(rank.p*len[k])),len[k])
      rank.corr <- spears[r]
      # r-1 correlations are < rank.corr, prob to have r-1 or less
      # corr < rank.corr is given by a binomial with success rate
      # equal to the probability to get a corr < rank.corr, i.e.,
      # cdf(rank.corr, RT.par). If rank.corr is high, it becomes
      # difficult to get as little as r-1 corr < rank.corr by chance!
      p.rt <- stats::pbinom(r-1, len[k], cdf(rank.corr, RT.par))
      #p.rp <- stats::pbinom(r-1, len[k], cdf(rank.corr, RP.par))
      #create new df melant prot et PTM
      #learnparam(df)
      #p.rtp
      p.rtp <- 1
      if(len.PTM[k] > 0){
        spears.PTM <- as.numeric(strsplit(spear.PTM[k],split=";")[[1]])
        r <- min(max(1,trunc(rank.p*len.PTM[k])),len.PTM[k])
        rank.corr.p <- spears.PTM[r]
        p.rtp <- stats::pbinom(r-1, len.PTM[k], cdfp(rank.corr.p, RP.par))
      }
      #p.rtp[pg] <- p.rp[pg]
      
      res <- rbind(res,data.frame(pairs[i,c("L","R")],
                                  LR.corr=pairs[i,"corr"], pw.id=pwid[k],
                                  pw.name=pwname[k], rank=r, len=len[k], len.PTM=len.PTM[k],
                                  rank.corr=rank.corr, target.genes=tg[k],
                                  target.corr=spear[k], target.genes.PTM=pg[k],
                                  target.PTM.corr=spear.PTM[k], pval=p.lr*p.rt*p.rtp,
                                  pvalLR=p.lr, pvalRT=p.rt, pvalRP=p.rtp,
                                  stringsAsFactors=FALSE))
    }
  }
  
  # avoid the impossible
  key <- paste(res$L, res$R, res$pw.id, sep="||")
  bad <- duplicated(key)
  res <- res[!bad,]
  
  # multiple hypothesis correction
  rawp <- res$pval
  adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
  res$qval <- adj$adjp[order(adj$index),fdr.proc]
  
  res
  
}  # .pValuesLRP

