#' Get regulated ligand-receptor pairs.
#'
#' Internal function to return all the pairs of ligands and receptors
#' having both a P-value below a given threshold.
#'
#' @param ds              A BSRDataModel object
#' @param cc              A BSRClusterComp object.
#' @param max.pval        The maximum P-value imposed to both the ligand
#' and the receptor.
#' @param min.logFC       The maximum log2 fold-change allowed for
#'   both the receptor and the ligand.
#' @param neg.receptors     A logical indicating whether receptors are only
#'   allowed to be upregulated (FALSE), or up- and downregulated (TRUE).
#' @param restrict.genes  A list of gene symbols that restricts ligands and
#'   receptors.
#'
#' @return A data frame containing putative ligand-receptor pairs along
#'   with the product of their respective P-values. This table is the first step
#'   of a ligand-receptor analysis.
#'
#' @details The \code{restrict.genes} parameter is used for special cases where
#'   LRdb must be further restricted to a subset.
#'   The putative ligand-receptor pairs has 6 columns : R, L, LR.pval, corr,
#'   L.logFC, and R.logFC.
#'
#' @importFrom methods is
#' @importFrom foreach %do% %dopar%
#' @import doParallel
#' @keywords internal
.getRegulatedLRPTM <- function(ds, cc, max.pval=0.01, min.logFC=1,
                            neg.receptors=FALSE, restrict.genes=NULL) {
 
  # local binding
  i <- NULL
  
  if ((max.pval <= 0) || (max.pval > 1))
    stop("max.pval must lie in ]0;1]")
  if (min.logFC <= 0)
    stop("min.logFC must be > 0")
  if (!is(ds, "BSRDataModelCompPTM"))
    stop("ds must be an object of class BSRDataModelCompPTM")
  if (!is(cc, "BSRClusterCompPTM"))
    stop("c must be an object of class BSRClusterComp")
  
  lrgenes <- intersect(c(LRdb$ligand,
                         LRdb$receptor), rownames(stats(cc)))
  if (!is.null(restrict.genes))
    lrgenes <- intersect(lrgenes, restrict.genes)
  
  # compute all the correlations at once
  corlr <- stats::cor(t(ncounts(ds)[lrgenes, c(colA(cc),colB(cc))]), method = "spearman")
  # get the pairs
  pairs <- NULL
  for (i in seq_len(nrow(LRdb))){
    L <- LRdb$ligand[i]
    R <- LRdb$receptor[i]
    if (L %in% lrgenes && R %in% lrgenes){
      pL <- stats(cc)[L, "pval"]
      pR <- stats(cc)[R, "pval"]
      if (pL <= max.pval && pR <= max.pval){
        fcL <- stats(cc)[L, "logFC"]
        fcR <- stats(cc)[R, "logFC"]
        if (!is.na(fcL) && !is.na(fcR) && fcL >= min.logFC &&
            ((neg.receptors && abs(fcR) >= min.logFC) || fcR >= min.logFC))
          pairs <- rbind(pairs,
                         data.frame(L=L, R=R,
                                    LR.pval=pL * pR, corr=corlr[L, R],
                                    L.logFC=fcL, R.logFC=fcR,
                                    stringsAsFactors=FALSE)
          )
      }
    }
  }

  if(is.null(pairs))
    stop("Dataframe `pairs` from `.getRegulatedLR` is NULL.")

  pairs
  
}  # .getRegulatedLR


#' Internal function to check receptor signaling downstream
#'
#' @param lr              A data frame as returned by
#'   \code{.getRegulatedLR()}.
#' @param pw              A table defining the reference pathways.
#' @param pw.size         A named vector with pathway sizes (names are pathway
#'   IDs).
#' @param rncounts        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param PTM        A matrix of normalized read counts with at
#'   least all the ligands, receptors, and genes in the reference pathways.
#' @param stats           A data.frame with a column 'pval' and
#' \code{rownames(stats)} assigned to gene symbols. Rows must at least
#'   include all the ligands, receptors, and genes in the reference pathways.
#' @param statsPTM           A data.frame with a column 'pval' and
#' \code{rownames(statsPTM)} assigned to gene symbols. Rows must at least
#'   include all the ligands, receptors, and genes in the reference pathways.
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
#' @return A table reporting all the ligand-receptor pairs provided in \code{lr}
#'   along with the pathways found and data about target gene regulation
#'   P-values. Target gene correlations with the receptor are computed as
#'   additional information provided \code{rncounts} is set.
#'
#' @importFrom foreach %do% %dopar%
#' @keywords internal
.downstreamRegulatedSignalingPTM <- function(lr, pw, pw.size, rncounts, stats, PTM, statsPTM,
                                          id.col, gene.col, pw.col,
                                          min.positive, with.complex = TRUE, single = FALSE) {

  if (!is.data.frame(stats))
    stop("stats must be a data.frame")
  if (!('pval' %in% names(stats)))
    stop("stats must contain a column named 'pval'")
  
  # local binding
  r <- p <- pl <- id <- NULL
  
  # define interaction types
  control.int <- "controls-expression-of"
  incomplex.int <- c("in-complex-with","interacts-with")
  directed.int <- c("controls-state-change-of", "catalysis-precedes",
                    "controls-expression-of", "controls-transport-of",
                    "controls-phosphorylation-of")
  addPTM.int <- "controls-phosphorylation-of"
  addPTM.int.pwc <- c("controls-phosphorylation-of", "regulates-phosphorylation-of")
  rmPTM.int.pwc <- c("controls-dephosphorylation-of", "regulates-dephosphorylation-of")
  
  if (with.complex)
    correlated.int <- union(control.int, incomplex.int)
  else
    correlated.int <- control.int
  
  # cat("\n rownames(rncounts) : ")
  # cat(rownames(rncounts)[1:3])
  # cat("\n rownames(PTM) : ")
  # cat(rownames(PTM)[1:3])
  #rncounts <- rncounts[rownames(rncounts) %in% rownames(PTM),]
  rncounts <- rncounts[,colnames(rncounts) %in% colnames(PTM)]
  rncounts <- rncounts[,order(as.vector(colnames(rncounts)))]
  #rncounts <- rncounts[order(rownames(rncounts)),]
  PTM <- PTM[,colnames(PTM) %in% colnames(rncounts)]
  PTM <- PTM[,order(as.vector(colnames(PTM)))]
  
  # compute downstream correlations
  corrg <- stats::cor(t(rncounts), method = "spearman")
  corrgp <- stats::cor(t(rncounts), t(PTM), method = "spearman")
  # the global computation above is faster than restricted to the receptors
  corrg <- corrg[unique(lr$R), ]
  corrgp <- corrgp[unique(lr$R), ]
  
  # check each putative LR pair, loop over the receptors
  reg.proc <- foreach::foreach(r=unique(lr$R),.combine=rbind) %do% {
    # reg.proc <- NULL
    # for (r in unique(lr$putative.pairs$R)){
    
    # loop over the pathways containing the receptor r
    pa <- intersect(pw[pw[[gene.col]]==r,id.col],names(pw.size))
    if (length(pa)>0){
      receptor.ligands <- unique(lr$L[lr$R==r])
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
          
          # eliminate ligands of the receptor if present
          target.genes <- setdiff(target.genes, receptor.ligands)
          
          
          # cat("\n colnames(corrgp)[1:3] : ")
          # cat(colnames(corrgp)[1:3])
          
          target.genes.addPTM.name <- target.genes.addPTM.pos <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands,r) & target %in% target.genes & relationType %in% addPTM.int.pwc)$target)
          target.genes.rmPTM.name <- target.genes.rmPTM.pos <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands,r) & target %in% target.genes & relationType %in% rmPTM.int.pwc)$target)
          
          if(single){
            #c("\n pas single \n")
            target.genes.addPTM.pos <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands,r) & target %in% target.genes & relationType %in% addPTM.int.pwc)$genePos)
            target.genes.rmPTM.pos <- unique(subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands,r) & target %in% target.genes & relationType %in% rmPTM.int.pwc)$genePos)
          }
          # target.genes.bis <- subset(mergedPwCtrl, PathwayName==p & source %in% c(target.genes,receptor.ligands) & target %in% target.genes)$genePos
          # target.genes <- target.genes.bis
          
          # PTM.genes.detected !!!
          
          if (length(target.genes) >= min.positive){
            # if all conditions are met, list all target genes with
            # their regulation P-values in a data frame
            # row. Target genes are sorted wrt P-values in decreasing
            # order to keep the compatibility with correlation analysis,
            # where the most significant values are at the end.
            pv <- stats[target.genes, "pval"]
            o <- order(pv, decreasing=TRUE)
            pv <- pv[o]
            lfc <- stats[target.genes, "logFC"]
            lfc <- lfc[o]

            target.genes <- target.genes[o]
            c <- corrg[r, target.genes]
            
            #######

            
            ## PTM ##
            # cat("\n target.genes.addPTM.pos : ") #OK
            # cat(target.genes.addPTM.pos)
            #cat("\n")
            PTMTheoricalPres <- unique(target.genes.addPTM.pos[target.genes.addPTM.pos %in% colnames(corrgp)])
            # if(length(PTMTheoricalPres>0)){
            #   cat("\n PTMThPres \n")
            #   cat(unlist(PTMTheoricalPres))
            # }
            #si presence de ABC_NA on prend toutes les positions de ABC
            tmpNA <- grepl("_NA", target.genes.addPTM.pos, fixed = TRUE)
            if(sum(tmpNA > 0)){
              target.genes.addPTM.posNA <- target.genes.addPTM.pos[tmpNA]
              target.genes.addPTM.nameNA <- target.genes.addPTM.name[tmpNA]
              #on prend toutes les pos si NA
              # cat("\n target.genes.addPTM.posNA \n")
              # cat(unlist(target.genes.addPTM.posNA))
              listPos <- c()
              for(gp in target.genes.addPTM.nameNA){
                tmpCol <- grepl(paste0(gp,"_"), colnames(corrgp))
                tmpPo <- NA
                tmpPo <- colnames(corrgp)[tmpCol]
                listPos <- c(listPos, tmpPo)
              }
              if(!is.null(tmpPo) && length(tmpPo) > 0)
                PTMTheoricalPres <- c(PTMTheoricalPres, tmpPo)
              
              #target.genes.addPTM.pos <- PTMTheoricalPres
            }
            target.genes.addPTM.pos <- PTMTheoricalPres
            if(length(PTMTheoricalPres>0)){
              # cat("\n target.genes.addPTM.pos \n")
              # cat(unlist(target.genes.addPTM.pos))
            }
            
            ## rmPTM ##
            
            rmPTMTheoricalPres <- unique(target.genes.rmPTM.pos[target.genes.rmPTM.pos %in% colnames(corrgp)])
            if(length(rmPTMTheoricalPres>0)){
              # cat("\n rmPTMThPres \n")
              # cat(unlist(rmPTMTheoricalPres))
            }
            #si presence de ABC_NA on prend toutes les positions de ABC
            tmpNA <- grepl("_NA", target.genes.rmPTM.pos, fixed = TRUE)
            if(sum(tmpNA > 0)){
              target.genes.rmPTM.posNA <- target.genes.rmPTM.pos[tmpNA]
              target.genes.rmPTM.nameNA <- target.genes.rmPTM.name[tmpNA]
              #on prend toutes les pos si NA
              # cat("\n target.genes.rmPTM.posNA \n")
              # cat(unlist(target.genes.rmPTM.posNA))
              listPos <- c()
              for(gp in target.genes.rmPTM.nameNA){
                tmpCol <- grepl(paste0(gp,"_"), colnames(corrgp))
                tmpdePo <- NA
                tmpdePo <- colnames(corrgp)[tmpCol]
                listPos <- c(listPos, tmpdePo)
              }
              if(!is.null(tmpdePo) && length(tmpdePo) > 0)
                rmPTMTheoricalPres <- c(rmPTMTheoricalPres, tmpdePo)
              
              
              
              #target.genes.addPTM.pos <- PTMTheoricalPres
            }
            target.genes.rmPTM.pos <- rmPTMTheoricalPres
            # if(length(rmPTMTheoricalPres>0)){
            #   cat("\n target.genes.rmPTM.pos \n")
            #   cat(unlist(target.genes.rmPTM.pos))
            # }
            # target.genes.addPTM.posNA <- target.genes.addPTM.pos[tmpNA]
            # target.genes.addPTM.nameNA <- target.genes.addPTM.name[tmpNA]
            # #on prend toutes les pos si NA
            # if(length(target.genes.addPTM.posNA) > 0){
            #   listPos <- c()
            #   for(gp in target.genes.addPTM.nameNA){
            #     tmpCol <- grepl(paste0(gp,"_"), colnames(corrgp))
            #     tmpPo <- colnames(corrgp)[tmpCol]
            #     listPos <- c(listPos, tmpPo)
            #   }
            #   PTMTheoricalPres <- c(PTMTheoricalPres, tmpPo)
            # }
            # target.genes.addPTM <- PTMTheoricalPres

            if(length(PTMTheoricalPres)>0 && r %in% rownames(corrgp)){
              pvp <- statsPTM[target.genes.addPTM.pos, "pval"]
              op <- order(pvp, decreasing=TRUE)
              pvp <- pvp[op]
              lfcp <- statsPTM[target.genes.addPTM.pos, "logFC"]
              lfcp <- lfcp[op]
              
              target.genes.addPTM.pos <- target.genes.addPTM.pos[op]
              cp <- corrgp[r, target.genes.addPTM.pos]
              len.addPTM <- length(cp)
              # cat("\n pTP \n")
            }
            else{
              cp <- pvp <- lfcp <- target.genes.addPTM.pos <- c(NA)
              len.addPTM <- 0
              #cat("\n pas pTP \n")
            }
            
            
            if(length(rmPTMTheoricalPres)>0 && r %in% rownames(corrgp)){
              pvdp <- statsPTM[target.genes.rmPTM.pos, "pval"]
              odp <- order(pvdp, decreasing=TRUE)
              pvdp <- pvdp[odp]
              lfcdp <- statsPTM[target.genes.rmPTM.pos, "logFC"]
              lfcdp <- lfcdp[odp]
              
              target.genes.rmPTM.pos <- target.genes.rmPTM.pos[odp]
              cdp <- corrgp[r, target.genes.rmPTM.pos]
              len.rmPTM <- length(cdp)
              # cat("\n dpTP \n")
            }
            else{
              cdp <- pvdp <- lfcdp <- target.genes.rmPTM.pos <- c(NA)
              len.rmPTM <- 0
              #cat("\n pas pTP \n")
            }
            pvPTM <- c(pvp, pvdp)
            lfcPTM <- c(lfcp, lfcdp)
            cPTM <- c(cp, cdp)
            target.genes.addPTM.pos <- c(target.genes.addPTM.pos, target.genes.rmPTM.pos)
            
            pvPTM <- pvPTM[!is.na(pvPTM)]
            lfcPTM <- lfcPTM[!is.na(lfcPTM)]
            cPTM <- cPTM[!is.na(cPTM)]
            target.genes.addPTM.pos <- target.genes.addPTM.pos[!is.na(target.genes.addPTM.pos)]
            
            len.PTM <- len.addPTM + len.rmPTM
            
            #########
            
            # PTM.genes <- target.genes.addPTM.pos[o]
            # cp <- corrgp[r, PTM.genes]
            
            data.frame(pathway=p, target.pval=paste(pv,collapse=";"),
                       target.genes=paste(target.genes,collapse=";"),
                       target.corr=paste(c, collapse=";"),
                       target.logFC=paste(lfc, collapse=";"),
                       len=length(c), 
                       
                       PTM.pval=paste(pvPTM,collapse=";"),
                       PTM.genes=paste(target.genes.addPTM.pos,collapse=";"),
                       PTM.corr=paste(cPTM, collapse=";"),
                       PTM.logFC=paste(lfcPTM, collapse=";"),
                       len.PTM=len.PTM,
                       
                       addPTM.pval=paste(pvp,collapse=";"),
                       addPTM.genes=paste(target.genes.addPTM.pos,collapse=";"),
                       addPTM.corr=paste(cp, collapse=";"),
                       addPTM.logFC=paste(lfcp, collapse=";"),
                       len.addPTM=len.addPTM,
                       
                       rmPTM.pval=paste(pvdp,collapse=";"),
                       rmPTM.genes=paste(target.genes.rmPTM.pos,collapse=";"),
                       rmPTM.corr=paste(cdp, collapse=";"),
                       rmPTM.logFC=paste(lfcdp, collapse=";"),
                       len.rmPTM=len.rmPTM, stringsAsFactors=FALSE)
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
                   target.pval=paste(best.2nd$target.pval, collapse='|'),
                   target.genes=paste(best.2nd$target.genes, collapse='|'),
                   target.corr=paste(best.2nd$target.corr, collapse='|'),
                   target.logFC=paste(best.2nd$target.logFC, collapse='|'),
                   len=paste(best.2nd$len, collapse='|'),
                   
                   PTM.pval=paste(best.2nd$PTM.pval, collapse='|'),
                   PTM.genes=paste(best.2nd$PTM.genes, collapse='|'),
                   PTM.corr=paste(best.2nd$PTM.corr, collapse='|'),
                   PTM.logFC=paste(best.2nd$PTM.logFC, collapse='|'),
                   len.PTM=paste(best.2nd$len.PTM, collapse='|'),
                   
                   addPTM.pval=paste(best.2nd$addPTM.pval, collapse='|'),
                   addPTM.genes=paste(best.2nd$addPTM.genes, collapse='|'),
                   addPTM.corr=paste(best.2nd$addPTM.corr, collapse='|'),
                   addPTM.logFC=paste(best.2nd$addPTM.logFC, collapse='|'),
                   len.addPTM=paste(best.2nd$len.addPTM, collapse='|'),
                   
                   rmPTM.pval=paste(best.2nd$rmPTM.pval, collapse='|'),
                   rmPTM.genes=paste(best.2nd$rmPTM.genes, collapse='|'),
                   rmPTM.corr=paste(best.2nd$rmPTM.corr, collapse='|'),
                   rmPTM.logFC=paste(best.2nd$rmPTM.logFC, collapse='|'),
                   len.rmPTM=paste(best.2nd$len.rmPTM, collapse='|'),
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
  conf.pairs$target.pval <- reg.proc[conf.pairs$R, "target.pval"]
  conf.pairs$len <- reg.proc[conf.pairs$R, "len"]
  conf.pairs$target.genes <- reg.proc[conf.pairs$R, "target.genes"]
  conf.pairs$target.corr <- reg.proc[conf.pairs$R, "target.corr"]
  conf.pairs$target.logFC <- reg.proc[conf.pairs$R, "target.logFC"]
  
  conf.pairs$PTM.pval <- reg.proc[conf.pairs$R, "PTM.pval"]
  conf.pairs$len.PTM <- reg.proc[conf.pairs$R, "len.PTM"]
  conf.pairs$PTM.genes <- reg.proc[conf.pairs$R, "PTM.genes"]
  conf.pairs$PTM.corr <- reg.proc[conf.pairs$R, "PTM.corr"]
  conf.pairs$PTM.logFC <- reg.proc[conf.pairs$R, "PTM.logFC"]
  
  conf.pairs$addPTM.pval <- reg.proc[conf.pairs$R, "addPTM.pval"]
  conf.pairs$len.addPTM <- reg.proc[conf.pairs$R, "len.addPTM"]
  conf.pairs$addPTM.genes <- reg.proc[conf.pairs$R, "addPTM.genes"]
  conf.pairs$addPTM.corr <- reg.proc[conf.pairs$R, "addPTM.corr"]
  conf.pairs$addPTM.logFC <- reg.proc[conf.pairs$R, "addPTM.logFC"]
  
  conf.pairs$rmPTM.pval <- reg.proc[conf.pairs$R, "rmPTM.pval"]
  conf.pairs$len.rmPTM <- reg.proc[conf.pairs$R, "len.rmPTM"]
  conf.pairs$rmPTM.genes <- reg.proc[conf.pairs$R, "rmPTM.genes"]
  conf.pairs$rmPTM.corr <- reg.proc[conf.pairs$R, "rmPTM.corr"]
  conf.pairs$rmPTM.logFC <- reg.proc[conf.pairs$R, "rmPTM.logFC"]
  
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
  
  conf.pairs[,c("L", "R", "LR.pval", "corr", "L.logFC", "R.logFC",
                "pwid", "pwname", "len", "target.genes",
                "target.pval", "target.logFC", "target.corr", 
                "len.PTM", "PTM.genes", "PTM.pval", "PTM.logFC", "PTM.corr",
                "len.addPTM", "addPTM.genes", "addPTM.pval", "addPTM.logFC", "addPTM.corr",
                "len.rmPTM", "rmPTM.genes", "rmPTM.pval", "rmPTM.logFC", "rmPTM.corr")]
  
}  # .downstreamRegulatedSignaling


#' Internal function to check receptor signaling downstream
#'
#' Assess the existence of concomitant regulations between a receptor,
#' part of a ligand-receptor pair, and
#' genes coding for proteins forming a complex with the receptor or genes
#' regulated by the receptor downstream signaling.
#'
#' @param ds              A BSRDataModel object
#' @param cc              A BSRClusterComp object.
#' @param lr              A table as returned by \code{.getRegulatedLR()}.
#' @param ds              An optional BSRDataModel object.
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
#' @return A data frame extending \code{lr} content with the pathways found to
#' contain the receptors and data about receptor target gene regulations
#' Strings in semi-colon-separated format are used to report
#' target genes and their regulation P-values in the
#' data frame. The target genes are sorted according to the P-values in
#' decreasing order.
#' 
#' In case \code{ds} is set, then correlations between the receptor and
#' target genes will be computed for documentation or additional use. The
#' row names of \code{stats(cc)} and \code{ncounts(ds)} must match
#' exactly (not necessarily in the same order).
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
.checkRegulatedReceptorSignalingPTM <- function(ds, cc, lr,
                                reference=c("REACTOME-GOBP","REACTOME","GOBP"),
                                max.pw.size=200, min.pw.size=5,
                                min.positive=4, restrict.pw=NULL,
                                with.complex=TRUE){
  
  if (!is(cc, "BSRClusterCompPTM"))
    stop("cc must be a BSRClusterComp object")
  if (!is(ds, "BSRDataModelCompPTM"))
    stop("ds must be a BSRDataModelCompPTM object")
  if (!all(rownames(stats(cc)) %in% rownames(ncounts(ds))))
    stop("The row names of stats(cc) and ncounts(ds) must match exactly")

  reference <- match.arg(reference)
  results <- list()
  
  # Reactome pathways
  if (reference %in% c("REACTOME-GOBP","REACTOME")){
    react <- reactome[reactome$`Gene name` %in% rownames(stats(cc)),]
    if (!is.null(restrict.pw))
      react <- react[react$`Reactome ID` %in% restrict.pw,]
    pw.size <- table(react$`Reactome ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- react[react$`Gene name` %in% lr$R, "Reactome ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         react[react$`Reactome ID` %in% names(pw.size), "Gene name"])
    )
    toMatch <- c(paste0(corgenes,"_"))
    #corgenesp <- unique(grep(paste(toMatch,collapse="|"), rownames(PTM(ds)), value=TRUE))
    if(single(ds))
      corgenesp <- unique(rownames(PTM(ds))[str_detect(rownames(PTM(ds)), paste(toMatch, collapse = "|"))])
    else
      corgenesp <- unique(rownames(PTM(ds)))
    
    results$reactome.pairs <- .downstreamRegulatedSignalingPTM(lr, react, pw.size,
                 ncounts(ds)[corgenes, c(colA(cc),colB(cc))], stats(cc)[corgenes,],
                 PTM(ds)[corgenesp, c(colA(cc),colB(cc))], statsPTM(cc)[corgenesp,], #/!\ corgenes doit etre gene.pos !!!
                 id.col="Reactome ID", gene.col="Gene name",
                 pw.col="Reactome name", min.positive=min.positive,
                 with.complex=with.complex, single = single(ds))
  }
  
  # GOBP
  if (reference %in% c("REACTOME-GOBP","GOBP")){
    go <- gobp[gobp$`Gene name` %in% rownames(stats(cc)),]
    if (!is.null(restrict.pw))
      go <- go[go$`GO ID` %in% restrict.pw,]
    pw.size <- table(go$`GO ID`)
    pw.size <- pw.size[pw.size>=min.pw.size & pw.size<=max.pw.size]
    contains.receptors <- go[go$`Gene name` %in% lr$R, "GO ID"]
    pw.size <- pw.size[names(pw.size) %in% contains.receptors]
    corgenes <- unique(c(lr$R,
                         go[go$`GO ID` %in% names(pw.size), "Gene name"])
    )
    
    toMatch <- c(paste0(corgenes,"_"))
    corgenesp <- unique(grep(paste(toMatch,collapse="|"), rownames(PTM(ds)), value=TRUE))
    
    results$gobp.pairs <- .downstreamRegulatedSignalingPTM(lr, go, pw.size,
                 ncounts(ds)[corgenes, c(colA(cc),colB(cc))], stats(cc)[corgenes,],
                 PTM(ds)[corgenesp, c(colA(cc),colB(cc))], statsPTM(cc)[corgenesp,],
                 id.col="GO ID", gene.col="Gene name", pw.col="GO name",
                 min.positive=min.positive, with.complex=with.complex, single = single(ds))
  }
  
  # merge
  if (reference == "REACTOME-GOBP"){
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
  else if (reference == "REACTOME")
    results$merged.pairs <- results$reactome.pairs
  else
    results$merged.pairs <- results$gobp.pairs
  
  results$merged.pairs
  
} # .checkRegulatedReceptorSignaling


#' Internal function to assign P-values to LR interactions
#'
#' Estimate the P-value of each ligand-receptor pair based
#' on the data frame output by \code{\link{.checkRegulatedReceptorSignaling}}.
#'
#' @param pairs         A data frame output by
#'   \code{checkRegulatedReceptorSignaling}.
#' @param param         A list containing the statistical model parameters.
#' @param rank.p        A number between 0 and 1 defining the rank of the last
#'   considered target genes.
#' @param fdr.proc      The procedure for adjusting P-values according to
#'   \code{\link[multtest]{mt.rawp2adjp}}.
#'
#' @return A data.frame with the data in \code{pairs} complemented with
#' P-values and adjusted P-values.
#' @keywords internal
.pValuesRegulatedLRPTM <- function(pairs, param, rank.p = 0.75,
                       fdr.proc = c("BH", "Bonferroni", "Holm", "Hochberg",
                                    "SidakSS", "SidakSD", "BY", "ABH", "TSBH")) {
  
  if (rank.p < 0 || rank.p > 1)
    stop("rank.p must lie in [0;1]")
  fdr.proc <- match.arg(fdr.proc)
  if(is.null(pairs))
    stop("Dataframe `pairs` from `.checkRegulatedReceptorSignaling` is NULL.")

  # estimate P-values
  res <- NULL
  resp <- NULL
  for (i in seq_len(nrow(pairs))){
    # all the data related to each pathway containing a given
    # receptor were collapsed separated by |
    # we need to split those pathways
    pwid <- unlist(strsplit(pairs$pwid[i],split="\\|"))
    pwname <- unlist(strsplit(pairs$pwname[i],split="\\|"))
    tg <- unlist(strsplit(pairs$target.genes[i],split="\\|"))
    spval <- unlist(strsplit(pairs$target.pval[i],split="\\|"))
    slfc <- unlist(strsplit(pairs$target.logFC[i],split="\\|"))
    spear <- unlist(strsplit(pairs$target.corr[i],split="\\|"))
    len <- as.numeric(unlist(strsplit(pairs$len[i],split="\\|")))
    
    PTMg <- unlist(strsplit(pairs$PTM.genes[i],split="\\|"))
    spvalPTM <- unlist(strsplit(pairs$PTM.pval[i],split="\\|"))
    slfcPTM <- unlist(strsplit(pairs$PTM.logFC[i],split="\\|"))
    spearPTM <- unlist(strsplit(pairs$PTM.corr[i],split="\\|"))
    lenPTM <- as.numeric(unlist(strsplit(pairs$len.PTM[i],split="\\|")))
    
    pg <- unlist(strsplit(pairs$addPTM.genes[i],split="\\|"))
    spvalp <- unlist(strsplit(pairs$addPTM.pval[i],split="\\|"))
    slfcp <- unlist(strsplit(pairs$addPTM.logFC[i],split="\\|"))
    spearp <- unlist(strsplit(pairs$addPTM.corr[i],split="\\|"))
    lenp <- as.numeric(unlist(strsplit(pairs$len.addPTM[i],split="\\|")))

    dpg <- unlist(strsplit(pairs$rmPTM.genes[i],split="\\|"))
    spvaldp <- unlist(strsplit(pairs$rmPTM.pval[i],split="\\|"))
    slfcdp <- unlist(strsplit(pairs$rmPTM.logFC[i],split="\\|"))
    speardp <- unlist(strsplit(pairs$rmPTM.corr[i],split="\\|"))
    lendp <- as.numeric(unlist(strsplit(pairs$len.rmPTM[i],split="\\|")))
    
    # get the LR correlation P-value
    p.lr <- pairs$LR.pval[i]

    # estimate the target gene correlation P-value based on rank statistics
    # for the individual correlation Gaussian model
    for (k in seq_len(length(len))){
      spvals <- as.numeric(strsplit(spval[k],split=";")[[1]])
      spears <- as.numeric(strsplit(spear[k],split=";")[[1]])
      slfcs <- as.numeric(strsplit(slfc[k],split=";")[[1]])
      r <- min(max(1,trunc(rank.p*len[k])),len[k])
      rank.pval <- spvals[r]
      rank.corr <- spears[r]
      # r-1 P-values are > rank.pval, prob to have r-1 or less
      # P-values > rank.pval is given by a binomial with success rate
      # equal to the probability to get a P-value > rank.pval, i.e.,
      # 1-rank.pval. If rank.pval is low (i.e., highly significant),
      # it becomes difficult to get as little as r-1 P-values > rank.pval by chance!
      p.rt <- stats::pbinom(r-1, len[k], 1-rank.pval) # cdf is punif here!
      #cat("1 \n")
      res <- rbind(res,data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
                                  pw.id=pwid[k], pw.name=pwname[k], rank=r,
                                  len=len[k], rank.pval=rank.pval,
                                  rank.corr=rank.corr,
                                  target.genes=tg[k], target.pval=spval[k],
                                  target.logFC=slfc[k], target.corr=spear[k],
                                  pvalLRT=p.lr*p.rt, stringsAsFactors=FALSE))
    }
    if(length(PTMg) > 0){
      for (k in seq_len(length(lenPTM))){
        spvalsPTM <- as.numeric(strsplit(spvalPTM[k],split=";")[[1]])
        spearsPTM <- as.numeric(strsplit(spearPTM[k],split=";")[[1]])
        slfcsPTM <- as.numeric(strsplit(slfcPTM[k],split=";")[[1]])
        r <- min(max(1,trunc(rank.p*lenPTM[k])),lenPTM[k])
        rank.pvalPTM <- spvalsPTM[r]
        rank.corrPTM <- spearsPTM[r]
        # r-1 P-values are > rank.pval, prob to have r-1 or less
        # P-values > rank.pval is given by a binomial with success rate
        # equal to the probability to get a P-value > rank.pval, i.e.,
        # 1-rank.pval. If rank.pval is low (i.e., highly significant),
        # it becomes difficult to get as little as r-1 P-values > rank.pval by chance!
        p.rPTM <- 1
        if(lenPTM[k] == 0){
          NULL
          
        }
        else{
          # cat("2b \n")
          # cat(unlist(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")]))
          # cat("\n")
          # cat(pwid[k], "\n", pwname[k], "\n", r, "\n",
          #     lenp[k], "\n", rank.pvalPTM, "\n",
          #     rank.corrp, "\n",
          #     pg[k], "\n", spvalp[k], "\n",
          #     slfcp[k], "\n", spearp[k], "\n",
          #     p.rp)
          p.rPTM <- stats::pbinom(r-1, lenPTM[k], 1-rank.pvalPTM) # cdf is punif here!
          
          posPh <- unlist(strsplit(PTMg[k],split=";")) %in% unlist(strsplit(pg[k],split=";"))
          posDePh <- unlist(strsplit(PTMg[k],split=";")) %in% unlist(strsplit(dpg[k],split=";"))
          # cat("\n posPh \n")
          # cat(unlist(posPh))
          # cat("\n posDePh \n")
          # cat(unlist(posDePh))
          
          lenp <- sum(posPh)
          if(lenp > 0){
            #rank.pvalp <- rank.pvalPTM[posPh] #a modifier
            #rank.corrp <- rank.corrPTM[posPh] #a modifier
            addPTM.genes <- unlist(strsplit(PTMg[k],split=";"))[posPh]
            addPTM.genes <- paste(addPTM.genes,collapse=";")
            addPTM.pval <- spvalsPTM[posPh]
            addPTM.pval <- paste(addPTM.pval,collapse=";")
            addPTM.logFC <- slfcsPTM[posPh]
            addPTM.logFC <- paste(addPTM.logFC,collapse=";")
            #PTM.logFC <- paste(PTM.logFC,collapse=";")
            addPTM.corr <- spearsPTM[posPh]
            addPTM.corr <- paste(addPTM.corr,collapse=";")
            #pvalRP <- p.rPTM[posPh]
          }
          else{
            #rank.pvalp <- NA
            #rank.corrp <- NA
            addPTM.genes <- NA
            addPTM.pval <- 1
            addPTM.logFC <- NA
            addPTM.corr <- NA
            #pvalRP <- 1
          }
          
          lendp <- sum(posDePh)
          if(lendp > 0){
            #rank.pvaldp <- rank.pvalPTM[posDePh]
            #rank.corrdp <- rank.corrPTM[posDePh]
            rmPTM.genes <- unlist(strsplit(PTMg[k],split=";"))[posDePh]
            rmPTM.genes <- paste(rmPTM.genes,collapse=";")
            rmPTM.pval <- spvalsPTM[posDePh]
            rmPTM.pval <- paste(rmPTM.pval,collapse=";")
            rmPTM.logFC <- slfcsPTM[posDePh]
            rmPTM.logFC <- paste(rmPTM.logFC,collapse=";")
            rmPTM.corr <- spearsPTM[posDePh]
            rmPTM.corr <- paste(rmPTM.corr,collapse=";")
            #pvalRdP <- p.rPTM[posDePh]
          }
          else{
            #rank.pvaldp <- NA
            #rank.corrdp <- NA
            rmPTM.genes <- NA
            rmPTM.pval <- 1
            rmPTM.logFC <- NA
            rmPTM.corr <- NA
            #pvalRdP <- 1
          }
          
# 
#           cat("\n 4 \n")
#           cat(unlist(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")]))
#           cat("\n \n PTM \n")
#           cat("pwid[k]-", pwid[k], "\n", 
#               "pwname[k]-",pwname[k], "\n", 
#               "r-", r, "\n",
#               "lenPTM[k]-",lenPTM[k], "\n", 
#               "rank.pvalPTM-", rank.pvalPTM, "\n",
#               "rank.corrPTM-", rank.corrPTM, "\n",
#               "PTMg[k]-", PTMg[k], " : ", length(PTMg[k]), "\n", 
#               "spvalPTM[k]-", spvalPTM[k], " : ", length(spvalPTM[k]), "\n", 
#               "spvalsPTM-", spvalsPTM, " : ", length(spvalsPTM), "\n", 
#               "slfcPTM[k]-", slfcPTM[k], " : ", length(slfcPTM[k]), "\n", 
#               "slfcsPTM-", slfcsPTM, " : ", length(slfcsPTM),  "\n", #pb la
#               "spearsPTM-", spearsPTM, "\n",
#               p.rPTM) 
#           cat("\n \n p \n")
#           cat(lenp, "\n", #pas coherent
#               # rank.pvalp, "\n",#pas coherent
#               # rank.corrp, "\n",#pas coherent
#               PTM.genes, " : ", length(PTM.genes), "\n", 
#               PTM.pval, " : ", length(PTM.pval), "\n", 
#               PTM.logFC, " : ", length(PTM.logFC), "\n", 
#               PTM.corr, "\n")#pas coherent
#           cat("\n \n dp \n")
#           cat(lendp, "\n", 
#               # rank.pvaldp, "\n",
#               # rank.corrdp, "\n",
#               rmPTM.genes, "\n", rmPTM.pval, "\n",
#               rmPTM.logFC, "\n", rmPTM.corr, "\n")
#           cat(nrow(resp),",", nrow(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")]), ",")
#           cat(nrow(data.frame(
#                                 pw.id=pwid[k], pw.name=pwname[k], rank=r,
#                                 lenPTM=lenPTM[k], rank.pvalPTM=rank.pvalPTM,
#                                 rank.corrPTM=rank.corrPTM,
#                                 PTM.genes=PTMg[k], PTM.pval=spvalPTM[k],
#                                 PTM.logFC=slfcPTM[k], PTM.corr=spearPTM[k],
#                                 pvalRPTM=p.rPTM,
#                                 
#                                 lenp=lenp, 
#                                 # rank.pvalp=rank.pvalp,
#                                 # rank.corrp=rank.corrp,
#                                 PTM.genes=PTM.genes, PTM.pval=PTM.pval,
#                                 PTM.logFC=PTM.logFC, PTM.corr=PTM.corr,
#                                 #pvalRP=pvalRP,
#                                 
#                                 lendp=lendp, 
#                                 # rank.pvaldp=rank.pvaldp,
#                                 # rank.corrdp=rank.corrdp,
#                                 rmPTM.genes=rmPTM.genes, rmPTM.pval=rmPTM.pval,
#                                 rmPTM.logFC=rmPTM.logFC, rmPTM.corr=rmPTM.corr,
#                                 #pvalRdP=pvalRdP,
#                                 stringsAsFactors=FALSE)))
          resp <- rbind(resp,data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
                                        pw.id=pwid[k], pw.name=pwname[k], rank=r,
                                        lenPTM=lenPTM[k], rank.pvalPTM=rank.pvalPTM,
                                        rank.corrPTM=rank.corrPTM,
                                        PTM.genes=PTMg[k], PTM.pval=spvalPTM[k],
                                        PTM.logFC=slfcPTM[k], PTM.corr=spearPTM[k],
                                        pvalRPTM=p.rPTM,
                                        
                                        lenp=lenp, 
                                        # rank.pvalp=rank.pvalp,
                                        # rank.corrp=rank.corrp,
                                        addPTM.genes=addPTM.genes, addPTM.pval=addPTM.pval,
                                        addPTM.logFC=addPTM.logFC, addPTM.corr=addPTM.corr,
                                        #pvalRP=pvalRP,
                                        
                                        lendp=lendp, 
                                        # rank.pvaldp=rank.pvaldp,
                                        # rank.corrdp=rank.corrdp,
                                        rmPTM.genes=rmPTM.genes, rmPTM.pval=rmPTM.pval,
                                        rmPTM.logFC=rmPTM.logFC, rmPTM.corr=rmPTM.corr,
                                        #pvalRdP=pvalRdP,
                                        stringsAsFactors=FALSE))
        }
        # cat("2 \n")
        # cat(unlist(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")]))
        # cat("\n")
        # cat(pwid[k], "\n", pwname[k], "\n", r, "\n",
        #     lenp[k], "\n", rank.pvalp, "\n",
        #     rank.corrp, "\n",
        #     pg[k], "\n", spvalp[k], "\n",
        #     slfcp[k], "\n", spearp[k], "\n",
        #     p.rp)
        # resp <- rbind(resp,data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
        #                             pw.id=pwid[k], pw.name=pwname[k], rank=r,
        #                             lenp=lenp[k], rank.pvalp=rank.pvalp,
        #                             rank.corrp=rank.corrp,
        #                             PTM.genes=pg[k], PTM.pval=spvalp[k],
        #                             PTM.logFC=slfcp[k], PTM.corr=spearp[k],
        #                             pvalRP=p.rp, stringsAsFactors=FALSE))
        
      }
    }
    else{
      # for (k in seq_len(length(lenp))){
      #   spvalsp <- as.numeric(strsplit(spvalp[k],split=";")[[1]])
      #   spearsp <- as.numeric(strsplit(spearp[k],split=";")[[1]])
      #   slfcsp <- as.numeric(strsplit(slfcp[k],split=";")[[1]])
      #   r <- min(max(1,trunc(rank.p*lenp[k])),lenp[k])
      #   rank.pvalp <- spvalsp[r]
      #   rank.corrp <- spearsp[r]
      #   # r-1 P-values are > rank.pval, prob to have r-1 or less
      #   # P-values > rank.pval is given by a binomial with success rate
      #   # equal to the probability to get a P-value > rank.pval, i.e.,
      #   # 1-rank.pval. If rank.pval is low (i.e., highly significant),
      #   # it becomes difficult to get as little as r-1 P-values > rank.pval by chance!
      #   p.rp <- 1 #stats::pbinom(r-1, lenp[k], 1-rank.pvalp) # cdf is punif here!
      #   if(length(pg[k]) == 0){
      #     p.rp <- 1
      #   }
      #   cat("3 \n")
      # cat("\n", ncol(resp), "\n", colnames(resp))
      # cat("\n", ncol(data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
      #                           pw.id=pwid[k], pw.name=pwname[k], rank=r,
      #                           
      #                           lenPTM=0, rank.pvalPTM=NA,
      #                           rank.corrPTM=NA,
      #                           PTM.genes=NA, PTM.pval=1,
      #                           PTM.logFC=NA, PTM.corr=NA,
      #                           pvalRPTM=1,
      #                           
      #                           lenp=0, 
      #                           # rank.pvalp=NA,
      #                           # rank.corrp=NA,
      #                           PTM.genes=NA, PTM.pval=1,
      #                           PTM.logFC=NA, PTM.corr=NA,
      #                           pvalRP=1,
      #                           
      #                           lendp=0, 
      #                           # rank.pvaldp=NA,
      #                           # rank.corrdp=NA,
      #                           rmPTM.genes=NA, rmPTM.pval=1,
      #                           rmPTM.logFC=NA, rmPTM.corr=NA,
      #                           pvalRdP=1,stringsAsFactors=FALSE)), "\n", colnames(data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
      #                                                                                         pw.id=pwid[k], pw.name=pwname[k], rank=r,
      #                                                                                         
      #                                                                                         lenPTM=0, rank.pvalPTM=NA,
      #                                                                                         rank.corrPTM=NA,
      #                                                                                         PTM.genes=NA, PTM.pval=1,
      #                                                                                         PTM.logFC=NA, PTM.corr=NA,
      #                                                                                         pvalRPTM=1,
      #                                                                                         
      #                                                                                         lenp=0, 
      #                                                                                         # rank.pvalp=NA,
      #                                                                                         # rank.corrp=NA,
      #                                                                                         PTM.genes=NA, PTM.pval=1,
      #                                                                                         PTM.logFC=NA, PTM.corr=NA,
      #                                                                                         pvalRP=1,
      #                                                                                         
      #                                                                                         lendp=0, 
      #                                                                                         # rank.pvaldp=NA,
      #                                                                                         # rank.corrdp=NA,
      #                                                                                         rmPTM.genes=NA, rmPTM.pval=1,
      #                                                                                         rmPTM.logFC=NA, rmPTM.corr=NA,
      #                                                                                         pvalRdP=1,stringsAsFactors=FALSE)))
        resp <- rbind(resp,data.frame(pairs[i,c("L","R","LR.pval","corr","L.logFC","R.logFC")],
                                      pw.id=pwid[k], pw.name=pwname[k], rank=r,
                                      
                                      lenPTM=0, rank.pvalPTM=NA,
                                      rank.corrPTM=NA,
                                      PTM.genes=NA, PTM.pval=1,
                                      PTM.logFC=NA, PTM.corr=NA,
                                      pvalRPTM=1,
                                      
                                      lenp=0, 
                                      # rank.pvalp=NA,
                                      # rank.corrp=NA,
                                      addPTM.genes=NA, addPTM.pval=1,
                                      addPTM.logFC=NA, addPTM.corr=NA,
                                      #pvalRP=1,
                                      
                                      lendp=0, 
                                      # rank.pvaldp=NA,
                                      # rank.corrdp=NA,
                                      rmPTM.genes=NA, rmPTM.pval=1,
                                      rmPTM.logFC=NA, rmPTM.corr=NA,
                                      stringsAsFactors=FALSE))
      # }
    }
    
  }
  names(res)[4] <- "LR.corr"
  names(resp)[4] <- "LR.corr"
  
  resMerged <- merge(res, resp, by=c("L","R","pw.id"), all=T) #maybe remove duplicated columns
  
  #cat("\n",unlist(colnames(resMerged)),"\n")
                          # L R pw.id LR.pval.x LR.corr.x L.logFC.x R.logFC.x pw.name.x rank.x 
                          # len rank.pval rank.corr target.genes target.pval 
                          # target.logFC target.corr pvalLRT 
                          # LR.pval.y LR.corr.y L.logFC.y R.logFC.y pw.name.y rank.y lenPTM rank.pvalPTM rank.corrPTM PTM.genes PTM.pval PTM.logFC PTM.corr pvalRPTM 
                          # lenp rank.pvalp rank.corrp PTM.genes PTM.pval PTM.logFC PTM.corr pvalRP 
                          # lendp rank.pvaldp rank.corrdp rmPTM.genes rmPTM.pval rmPTM.logFC rmPTM.corr pvalRdP 
  
  #resMerged <- resMerged[!duplicated(as.list(resMerged))]
  resMerged <- resMerged[,c("L", "R", "pw.id", "LR.pval.x", "LR.corr.x","L.logFC.x", "R.logFC.x", "pw.name.x", "rank.x",
                            "len", "rank.pval", "rank.corr", "target.genes", "target.pval", 
                            "target.logFC", "target.corr", "pvalLRT",
                            "lenPTM", "rank.pvalPTM", "rank.corrPTM", "PTM.genes", "PTM.pval", "PTM.logFC", "PTM.corr", "pvalRPTM",
                            "lenp", "addPTM.genes", "addPTM.pval", "addPTM.logFC", "addPTM.corr", 
                            "lendp", "rmPTM.genes", "rmPTM.pval", "rmPTM.logFC", "rmPTM.corr")]
  resMerged <- resMerged[!duplicated(resMerged), ]
  resMerged$pvalRPTM[is.na(resMerged$pvalRPTM)] <- 1
  resMerged$pval <- resMerged$pvalLRT*resMerged$pvalRPTM
  
  # avoid the impossible
  key <- paste(resMerged$L, resMerged$R, resMerged$pw.id, sep="||")
  bad <- duplicated(key)
  resMerged <- resMerged[!bad,]
  
  # multiple hypothesis correction
  rawp <- resMerged$pval
  adj <- multtest::mt.rawp2adjp(rawp,fdr.proc)
  resMerged$qval <- adj$adjp[order(adj$index),fdr.proc]
  #cat("\n",unlist(colnames(resMerged)))
                              # L R pw.id LR.pval.x LR.corr.x L.logFC.x R.logFC.x pw.name.x rank.x 
                              # len rank.pval rank.corr target.genes target.pval 
                              # target.logFC target.corr pvalLRT 
                              # LR.pval.y LR.corr.y L.logFC.y R.logFC.y pw.name.y rank.y lenPTM rank.pvalPTM rank.corrPTM PTM.genes PTM.pval PTM.logFC PTM.corr pvalRPTM 
                              # lenp rank.pvalp rank.corrp PTM.genes PTM.pval PTM.logFC PTM.corr pvalRP 
                              # rank.pvaldp rank.corrdp rmPTM.genes rmPTM.pval rmPTM.logFC pvalRdP pval qval
  
  resMerged <- resMerged[,c("L", "R", "pw.id", "LR.pval.x", "LR.corr.x","L.logFC.x", "R.logFC.x", "pw.name.x", "rank.x",
                            "len", "rank.pval", "rank.corr", "target.genes", "target.pval", 
                            "target.logFC", "target.corr", "pvalLRT",
                            "lenPTM", "rank.pvalPTM", "rank.corrPTM", "PTM.genes", "PTM.pval", "PTM.logFC", "PTM.corr", "pvalRPTM",
                            "lenp", "addPTM.genes", "addPTM.pval", "addPTM.logFC", "addPTM.corr", 
                            "lendp", "rmPTM.genes", "rmPTM.pval", "rmPTM.logFC", "rmPTM.corr", "pval", "qval")]
  
  colnames(resMerged) <- c("L", "R", "pw.id", "LR.pval", "LR.corr","L.logFC", "R.logFC", "pw.name", "rank",
                           "len", "rank.pval", "rank.corr", "target.genes", "target.pval", 
                           "target.logFC", "target.corr", "pvalLRT", 
                           "lenPTM", "rank.pvalPTM", "rank.corrPTM", "PTM.genes", "PTM.pval", "PTM.logFC", "PTM.corr", "pvalRPTM",
                           "lenp", "addPTM.genes", "addPTM.pval", "addPTM.logFC", "addPTM.corr",
                           "lendp", "rmPTM.genes", "rmPTM.pval", "rmPTM.logFC", "rmPTM.corr", "pval", "qval")

  resMerged
  
}  # .pValuesRegulatedLR
