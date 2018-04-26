# Collection of MS1 functions
# E. Schymanski, 16/3/2017

#' Extract Raw MS1 Scan from RMassBank Archive
#'
#' @description This extracts the "raw" MS1 data from a \code{RMassBank} archive file
#' and saves the corresponding MS1 files in the given directory.
#'
#' @usage extractMS1.RMBarchive(RMBarchive, MS1_dir)
#'
#' @param RMBarchive Archive file to process
#' @param MS1_dir Directory in which to save the MS1 files
#'
#' @return Returns nothing; files are saved externally.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{msmsWorkflow}}
#'
#' @export
#'
extractMS1.RMBarchive <- function(RMBarchive, MS1_dir) {
  # load the workspace
  w <- loadMsmsWorkspace(RMBarchive)
  #extract the information we need
  mzs <- lapply(w@spectra,function(s) s@parent@mz)
  ints <- lapply(w@spectra,function(s) s@parent@intensity)
  ids <- lapply(w@spectra,function(s) s@id)
  names <- lapply(w@spectra,function(s) s@name)
  # now save the MS1 files
  for (i in 1:length(mzs)) {
    name <- names[[i]]
    name <- sub(" ","",name,fixed=TRUE)
    id <- ids[[i]]
#    rt <- rts[[i]]
    file_name <- paste0(MS1_dir,"/",id,"_",name,"_MS1raw.txt")
    peaks <- cbind(mzs[[i]],ints[[i]])
    write.table(peaks,file=file_name, quote=FALSE, row.names=FALSE, col.names = FALSE)
  }
}

#' Extract Compound ID from File Path and File Name
#'
#' @description This extracts the Compound ID from the start of a
#' file name, given the full path.
#'
#' @usage IDfromFilePath(file_path)
#'
#' @param file_path The file path to extract the ID from.
#'
#' @return Returns the compound ID or whatever content is in front of the first \code{"_"}.
#' Designed to work off output of other functions.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{extractMS1.RMBarchive}}
#'
#' @export
#'
#' @examples
#' IDfromFilePath("C:/MY_DATA/MS1/3529_Indoxacarb_MS1raw.txt")
#'
IDfromFilePath <- function(file_path) {
  filename <- basename(file_path)
  id <- strsplit(filename,"_",fixed=TRUE)[[1]][1]
  return(id)
}


#' Select a Subset of Annotated Peaks for Plotting
#'
#' @description This uses the \code{\link{filterPeakSatellites}} of \code{\link{RMassBank}}
#' to select annotated peaks for plotting and thus prevent overlap.
#'
#' @usage trimAnnotation(annotated_dataframe,mzLimit=10,IntLimit=1,
#' column_names = c("mz","anno","intensity","relInt"), UseRmbDefaultSettings=TRUE)
#'
#' @param annotated_dataframe The data frame containing the annotated data to trim
#' @param mzLimit The \code{mz} window to exclude other annotations.
#' @param IntLimit The default value of 1 is recommended. See details.
#' @param column_names This should be used to defined the column names and can be retrieved from
#' \code{annotated_dataframe} using \code{\link{colnames}}. NOTE that the columns \code{mz} and
#' \code{intensity} MUST exist for this function to work.
#' @param UseRmbDefaultSettings The default \code{TRUE} calls \code{RmbDefaultSettings}; if
#' other settings are in use set this to \code{FALSE}. See Details.
#'
#' @return Returns a "trimmed" annotated_dataframe
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @details
#' \code{UseRmbDefaultSettings}: If \code{RMassBank} is not in use for any other purpose, the
#' default \code{TRUE} calls \code{RmbDefaultSettings} to ensure that settings are present for use.
#' If other settings have been defined using \code{RMassBank}, set this to \code{FALSE} to avoid
#' overwriting those settings. This function retrieves settings and adjusts them internally without
#' overwriting any stored settings.
#' \code{IntLimit}: This is used in \code{\link{filterPeakSatellites}} and is relative to 1.
#' For this context, this should be kept at 1. Values below 1 will let additional peaks through;
#' values above 1 lead to strange behaviour.
#'
#'
#' @seealso \code{\link{extractMS1.RMBarchive}}
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA000405.txt",package="ReSOLUTION")
#' MBrecord <- system.file("extdata","EA026206_Simazine.txt",package="ReSOLUTION")
#' anno_peaks <- getMBPeaksAnnos(MBrecord, fixColNames=TRUE)
#' anno_cols <- sub("Int","intensity",colnames(anno_peaks))
#' colnames(anno_peaks) <- anno_cols
#' trimAnnotation(anno_peaks, mzLimit=10, IntLimit=1, column_names=anno_cols)
#'
#'
trimAnnotation <- function(annotated_dataframe,mzLimit=10,IntLimit=1,
                           column_names = c("mz","anno","intensity","relInt"),
                           UseRmbDefaultSettings=TRUE) {
  # note this requires columns "mz" and "intensity" to work
  #column_names <- colnames(annotated_dataframe)
  if (UseRmbDefaultSettings) {
    RmbDefaultSettings()
  }
  filter_settings <- getOption("RMassBank")$filterSettings
  filter_settings$satelliteMzLimit <- mzLimit
  filter_settings$satelliteIntLimit <- IntLimit
  colnames(annotated_dataframe) <- column_names
  trimmed_frame <- filterPeakSatellites(annotated_dataframe,filterSettings=filter_settings)
  return(trimmed_frame)
}

#' Extract MS1 and MS2 Peaklists from RMassBank Archive
#'
#' @description This extracts the "raw" MS1 data and the processed MS2 data
#' from a \code{RMassBank} archive file
#' and saves the corresponding files in the given directories. In contrast
#' to \code{\link{extractMS1.RMBarchive}}, this only saves data where MS2 were found.
#'
#' @usage extractMS1and2.RMBarchive(RMBarchive, MS1_dir, MS2_dir)
#'
#' @param RMBarchive Archive file to process
#' @param MS1_dir Directory in which to save the MS1 files
#' @param MS2_dir Directory in which to save the MS2 files
#'
#' @return Returns nothing; files are saved externally.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{extractMS1.RMBarchive}}.
#'
#' @details Note that this function could be extended to improve the MS2 export.
#' The current state is very basic.
#'
#' @export
#'
extractMS1and2.RMBarchive <- function(RMBarchive, MS1_dir, MS2_dir) {
  # load the workspace
  w <- loadMsmsWorkspace(RMBarchive)
  #extract the information we need
  # indexes of the spectra found
  specs_found <- unlist(lapply(w@spectra,function(s) s@found))
  # the MS2 scans
  children <- lapply(w@spectra[specs_found], function(sp) sp@children)
  # the IDs and names
  ids <- lapply(w@spectra[specs_found], function(sp) sp@id)
  # at the moment, "names" aren't used in the file name
  names <- lapply(w@spectra[specs_found],function(s) s@name)
  # retrieve the data we need
  allDf <- lapply(children, function(c) lapply(c, getData))
  MS1 <- lapply(w@spectra[specs_found], function(sp) as.data.frame(sp@parent))

  # check dirs exist
  if (!file.exists(MS1_dir)) {
    dir.create(MS1_dir)
  }
  if (!file.exists(MS2_dir)) {
    dir.create(MS2_dir)
  }

  # save MS2 files
  # NOTE: could build in some extra functionality here to improve the export.
  for (i in 1:length(allDf)) {
    for (j in 1:length(allDf[[i]])) {
      id <- as.character(ids[i])
      name1 <- names(allDf[[i]][j])
      filename <- paste(id, name1,"MSMS.txt",sep="_")
      msms <- allDf[[i]][[j]]
      #msms <- msms[which(msms$intensity>50),]
      #msms_trim <- trimAnnotation(msms, mzLimit = 1, column_names=c("mz", "intensity"))

      write.table(msms, filename,quote=F, row.names=F)
      # add col.names=F to get rid of header...
    }
  }
  #save MS1 files
  for (i in seq_len(length(MS1))) {
    id <- as.character(ids[i])
    filename <- paste(id, "MS.txt",sep="_")
    ms <- MS1[[i]]
    write.table(ms, filename,quote=F, row.names=FALSE)
    # add col.names=F to get rid of header...
  }

  # # now save the MS1 files
  # for (i in 1:length(mzs)) {
  #   name <- names[[i]]
  #   name <- sub(" ","",name,fixed=TRUE)
  #   id <- ids[[i]]
  #   #    rt <- rts[[i]]
  #   file_name <- paste0(MS1_dir,"/",id,"_",name,"_MS1raw.txt")
  #   peaks <- cbind(mzs[[i]],ints[[i]])
  #   write.table(peaks,file=file_name, quote=FALSE, row.names=FALSE, col.names = FALSE)
  # }
}
