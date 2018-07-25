
#' EIC Prescreening function for RMassBank
#' 
#' @description This is a wrapper function to perform prescreening of
#' EICs and MS/MS scans prior to running RMassBank, so entries can 
#' be checked visually. Output is a PDF and summary CSV; this is designed
#' to check the data and not look pretty. 
#' This is mainly an internal function made public just in case it helps anyone. 
#' It is recommended to run this immediately prior to msmsWorkflow in RMassBank
#' once all the appropriate files are set up. 
#'
#' @param archive_name Name of the archive (compatible with RMassBank naming)
#' @param RMB_mode Mode in RMassBank 
#' @param FileList The FileList to be used in RMassBank 
#' @param cmpd_list The CompoundList to be used in RMassBank
#' @param ppm_limit_fine The fine limit for extracting the MS/MS for \code{findMsMsHR.mass} in 
#' RMassBank (ppm)
#' @param EIC_limit The (absolute) limit for EIC extraction. 
#'
#' @return Function returns a PDF and CSV summary file.
#' 
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @export
#'
#' @examples
#' 
#' Currently no examples implemented.
#' 
RMB_EIC_prescreen <- function(archive_name, RMB_mode, FileList, cmpd_list, 
                              ppm_limit_fine=10, EIC_limit=0.001) {
  pdf_title <- paste(archive_name, "_trim_EICscan.pdf",sep="")
  pdf(pdf_title)
  n_spec <- 0
  cmpd_RT_maxI <- ""
  msms_found <- ""
  rts <- 0
  max_I_prec <- ""
  cmpd_RT_maxI_min <- ""
  file_list <- read.csv(FileList, stringsAsFactors = FALSE)
  cmpd_info <- read.csv(cmpd_list,stringsAsFactors = FALSE)
  #i <- 1
  f <- openMSfile(file_list$Files[1])
  
  for (i in 1:length(file_list$ID)) {
    #f <- openMSfile(file_list[i,1])
    cpdID <- file_list$ID[i]
    n_spec <- n_spec+1
    smiles <- findSmiles(cpdID)
    #mz <- cmpd_info$mz[i] 
    mz <- as.numeric(findMz(cpdID, RMB_mode)[3])
    #eic <- findEIC(f,mz$mzCenter, limit=0.001)
    eic <- findEIC(f,mz, limit=EIC_limit)
    msms_found[n_spec] <- FALSE # set to false and set to TRUE if found
    msms <- findMsMsHR.mass(f, mz, 0.5, ppm(mz,ppm_limit_fine,p=TRUE))
    max_I_prec_index <- which.max(eic$intensity)
    cmpd_RT_maxI[n_spec] <- eic[max_I_prec_index,1]
    max_I_prec[n_spec] <- eic[max_I_prec_index,2]
    cmpd_RT_maxI_min[n_spec] <- as.numeric(cmpd_RT_maxI[n_spec])/60
    # Retrieve spectrum data
    plot.new()
    plot.window(range(eic$rt),range(eic$intensity))
    box()
    lines(eic$intensity ~ eic$rt)
    for(specs in msms) {
      if(specs@found == TRUE) {
        df <- do.call(rbind, lapply(specs@children, function(sp) c(sp@rt, intensity = max(sp@intensity))))
        lines(intensity ~ retentionTime, data=df, type='h',col= "blue")
        msms_found[n_spec] <- TRUE
      }
    }
    title(main=cpdID, xlab="RT (sec)", ylab="Intensity")
    text(as.numeric(cmpd_RT_maxI[n_spec]),as.numeric(max_I_prec[n_spec]),
         labels=as.numeric(cmpd_RT_maxI_min[n_spec]),pos=4)
    axis(1)
    axis(2)
    # find out where precursor is most intense
    gc()
    rts[i] <- (cmpd_RT_maxI[n_spec])
  }
  
  dev.off()
  
  #summary <- cbind(cmpd_ids,cmpd_info$mz,cmpd_RT_maxI,cmpd_RT_maxI_min,max_I_prec)
  write.csv(cbind(file_list$ID,cmpd_info$mz,cmpd_info$Name, cmpd_RT_maxI,cmpd_RT_maxI_min,max_I_prec,msms_found),
            file=paste(archive_name, "_RTs_trim.csv",sep=""),row.names=F)
  #file=paste(basename(mzML_file), "_RTs_wI.csv",sep=""),row.names=F)
}
