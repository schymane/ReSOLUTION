#MetFrag Interpreter Functions
# E. Schymanski, 7/4/2017



#' Convert ExplPeaks MetFrag Entry into MS/MS with mz and int
#'
#' @description This extracts the Explained MS/MS peaks from MetFrag results files
#'
#' @usage ExplPeaks2MSMS(ExplPeaks)
#'
#' @param ExplPeaks A single \code{ExplPeaks} entry from the MetFrag results file
#'
#' @return Returns a data frame with mz and int fields
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso
#'
#' @export
#'
#' @examples
#' ExplPeaks <- "63.9622_516.7;77.9655_999.0;82.0298_373.2"
#' ExplPeakList <- ExplPeaks2MSMS(ExplPeaks)
#'
ExplPeaks2MSMS <- function(ExplPeaks) {
  ExplPeaks_mzs <- strsplit(ExplPeaks, ";")[[1]]
  ExplPeaks_msms <- matrix(unlist(strsplit(ExplPeaks_mzs, "_")),ncol=2,byrow=T)#,dimnames=c("mz","int"))
  ExplPeaks_msms <- as.data.frame(ExplPeaks_msms)
  colnames(ExplPeaks_msms) <- c("mz","int")
  return(ExplPeaks_msms)
}

#' Convert FormulasOfExplPeaks MetFrag Entry into MS/MS with mz and formulas
#'
#' @description This extracts the Explained MS/MS peak formulas from MetFrag results files
#'
#' @usage ExplFormulas2AnnoMZ(FormulasOfExplPeaks)
#'
#' @param FormulasOfExplPeaks A single \code{FormulasOfExplPeaks} entry from the MetFrag results file
#'
#' @return Returns a data frame with mz and formula fields (latter labeled \code{anno}).
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso
#'
#' @export
#'
#' @examples
#' FormulasOfExplPeaks <- "62.9638:[HO2P]-H-;78.959:[HO3P]-H-;153.0325:[C4H10O4P]-"
#' ExplFormulas2AnnoMZ(FormulasOfExplPeaks)
#'
ExplFormulas2AnnoMZ <- function(FormulasOfExplPeaks) {
  ExplPeaksFormulas <- strsplit(FormulasOfExplPeaks, ";")[[1]]
  ExplPeaksFormulas <- matrix(unlist(strsplit(ExplPeaksFormulas, ":")),ncol=2,byrow=T)#,dimnames=c("mz","int"))
  ExplPeaksFormulas <- as.data.frame(ExplPeaksFormulas)
  colnames(ExplPeaksFormulas) <- c("mz","anno")
  return(ExplPeaksFormulas)
}

# more to come ... in development ...
