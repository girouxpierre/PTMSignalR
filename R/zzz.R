.onLoad <- function(...) {

 # Initially I was using envir = .GlobalEnv
 # Here and also in dataPrepare::resetLRdb
 # But it was causing error when running
 # examples during R CMD check.

 myEnv <- new.env(parent = globalenv())
 attach(myEnv, name="LRdbEnv")

 # The ligand-receptor database (LRdb) and the PathwayCommons interaction
 # table (PwC_ReactomeKEGG) used to come from SingleCellSignalR as exported
 # objects. Recent SingleCellSignalR versions no longer ship them, and recent
 # BulkSignalR keeps them in a private environment ('.SignalR') under
 # package-prefixed names (BulkSignalR_LRdb, BulkSignalR_Network). Resolve each
 # object robustly across the possible providers so PTMSignalR keeps working
 # whichever packaging is installed.

 # object from a package namespace or its data sets, else NULL
 .fromPkg <- function(name, pkg) {
   if (!requireNamespace(pkg, quietly = TRUE))
     return(NULL)
   obj <- tryCatch(get(name, envir = asNamespace(pkg)),
                   error = function(e) NULL)
   if (!is.null(obj))
     return(obj)
   tmp <- new.env()
   loaded <- try(utils::data(list = name, package = pkg, envir = tmp),
                 silent = TRUE)
   if (!inherits(loaded, "try-error") && exists(name, envir = tmp))
     return(get(name, envir = tmp))
   NULL
 }

 # object from BulkSignalR's private '.SignalR' env (recent BulkSignalR),
 # given the BulkSignalR-internal alias, else NULL
 .fromBSR <- function(alias) {
   if (!requireNamespace("BulkSignalR", quietly = TRUE))
     return(NULL)
   sig <- tryCatch(get(".SignalR", envir = asNamespace("BulkSignalR")),
                   error = function(e) NULL)
   if (!is.null(sig) && exists(alias, envir = sig, inherits = FALSE))
     return(get(alias, envir = sig))
   NULL
 }

 .getRefData <- function(name, bsr.alias) {
   obj <- .fromPkg(name, "SingleCellSignalR")           # legacy provider
   if (is.null(obj)) obj <- .fromBSR(bsr.alias)         # recent BulkSignalR
   if (is.null(obj)) obj <- .fromPkg(name, "BulkSignalR")
   if (is.null(obj))                                    # already attached?
     obj <- tryCatch(get(name, envir = globalenv(), inherits = TRUE),
                     error = function(e) NULL)
   if (is.null(obj))
     stop("Could not obtain '", name, "'. Load a package that provides it ",
          "(SingleCellSignalR, or a recent BulkSignalR which ships it as ",
          bsr.alias, ").")
   obj
 }

 assign("LRdb", .getRefData("LRdb", "BulkSignalR_LRdb"),
        envir = as.environment("LRdbEnv"))
 assign("PwC_ReactomeKEGG", .getRefData("PwC_ReactomeKEGG", "BulkSignalR_Network"),
        envir = as.environment("LRdbEnv"))
 data(sysdata, envir=environment())

}
