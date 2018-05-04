# MassBankRecordeR.R
# Functions to extract specific information from MassBank records
# 22/3/2017 Emma Schymanski
# Rewritten from getInfo scripts from Erik Mueller, contributions also
# from Michele Stravs and Tobias Schulze

# MBrecord <- "C:/DATA/MassSpecLibraries/MassBank/MassBankOpenData/Eawag/EA000405.txt"
# field_codes <- "default"
# colNames <- "default"
# directory <- "C:/DATA/Eawag/Benedikt/Scripts/recdata/"
# recursive=FALSE


#' Summarize MassBank Entries in Directory of MassBank Records
#'
#' @description This function returns a summary csv containing the entries in
#' all MassBank records in the given directory matching the list of \code{field_codes}.
#' Uses \code{\link{getMBRecordEntry}}.
#'
#' @usage getMBRecordInfo(directory, csv_name, field_codes="default", colNames="default",recursive=TRUE)
#'
#' @param directory The path to the directory to summarize. Subdirectories are also searched if
#' \code{recursive} is \code{TRUE}.
#' @param csv_name Name (with/without path) of the csv file to save the results.
#' @param field_codes A vector containing strings to match in the MassBank records. Setting
#' \code{"default"} will output default values, see \code{Details}.
#' @param colNames A vector containing column names to match \code{field_codes}. These must be
#' in the same order. If \code{"default"}, trimmed versions of the default \code{field_codes}
#' are used. If \code{"field_codes"}, then the entries in \code{field_codes} are used.
#' @param trimEntry If \code{TRUE}, the \code{field_code} is removed, if \code{FALSE},
#' the full line is returned.
#'
#' @details Use \code{\link{getMBRecordPeaks}}, \code{\link{getMBRecordPeakAnnotations}}
#' and/or \code{\link{getMBPeaksAnnos}} to retrieve peaks and peak annotations. This function
#' uses \code{\link{getMBRecordEntry}} to retrieve individual entries from records. If \code{NULL},
#' then \code{NA} is printed to the csv. If multiple entries match (e.g. for \code{CH$NAME}), then
#' one entry in the csv is made, separated with a pipe (\code{|}).
#' This function will exit with an error if the length of \code{field_codes} and \code{colNames}
#' do not match.
#' The default settings are:
#' \code{field_codes <- c("ACCESSION", "CH$NAME", "CH$SMILES", "CH$EXACT_MASS", "CH$FORMULA", "CH$IUPAC",
#' "CH$LINK: CAS", "CH$LINK: PUBCHEM CID", "CH$LINK: INCHIKEY", "CH$LINK: CHEMSPIDER",
#' "AC$MASS_SPECTROMETRY: COLLISION_ENERGY", "AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE",
#' "AC$CHROMATOGRAPHY: COLUMN_NAME", "AC$CHROMATOGRAPHY: RETENTION_TIME",
#' "MS$FOCUSED_ION: PRECURSOR_M/Z", "MS$FOCUSED_ION: PRECURSOR_TYPE")}
#' \code{colNames <- c("ACCESSION", "NAME", "SMILES", "NEUTRAL_EXACT_MASS", "FORMULA", "StdInChI",
#' "CAS_RN", "PUBCHEM_CID", "InChIKey", "CHEMSPIDER_ID", "COLLISION_ENERGY", "FRAGMENTATION_MODE",
#' "COLUMN_NAME", "RETENTION_TIME", "PRECURSOR_MZ", "PRECURSOR_TYPE")}
#'
#' @return Saves the results in \code{csvfile}, or \code{error}.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>, rewritten from getInfo
#' script from Erik Mueller.
#'
#' @seealso \code{\link{MBFileToVector}}, \code{\link{getMBRecordEntry}}, \code{\link{getMBPeaksAnnos}}.
#'
#' @export
#'
#' @examples
#' directory <- system.file("extdata", "recdata", package = "ReSOLUTION")
#' # note this overwrites any file of the same name ...
#' csv_name <- paste0(getwd(),"/test.csv")
#' getMBRecordInfo(directory, csv_name, field_codes="default", colNames="default",recursive=TRUE)
#'
getMBRecordInfo <- function(directory, csv_name, field_codes="default", colNames="default",
                            recursive=TRUE) {
  # set up the default values
  if (field_codes == "default") {
    field_codes <- c("ACCESSION", "CH$NAME", "CH$SMILES", "CH$EXACT_MASS", "CH$FORMULA", "CH$IUPAC",
                     "CH$LINK: CAS", "CH$LINK: PUBCHEM CID", "CH$LINK: INCHIKEY", "CH$LINK: CHEMSPIDER",
                     "AC$MASS_SPECTROMETRY: COLLISION_ENERGY", "AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE",
                     "AC$CHROMATOGRAPHY: COLUMN_NAME", "AC$CHROMATOGRAPHY: RETENTION_TIME",
                     "MS$FOCUSED_ION: PRECURSOR_M/Z", "MS$FOCUSED_ION: PRECURSOR_TYPE")
  }
  # if colNames = "field_codes", then set it to be the field_codes
  if (colNames == "field_codes") {
    colNames <- field_codes
  } else if (colNames == "default") {
    colNames <- c("ACCESSION", "NAME", "SMILES", "NEUTRAL_EXACT_MASS", "FORMULA", "StdInChI",
                  "CAS_RN", "PUBCHEM_CID", "InChIKey", "CHEMSPIDER_ID",
                  "COLLISION_ENERGY", "FRAGMENTATION_MODE",
                  "COLUMN_NAME", "RETENTION_TIME", "PRECURSOR_MZ", "PRECURSOR_TYPE")
  }

  # run a check that the length is the same, else it won't work.
  length_check <- length(colNames) == length(field_codes)
  if (!length_check) {
    error("Length of field_codes and colNames don't match; please check and try again")
  }
  n_cols <- length(colNames)
  # then extract the info
  files <- list.files(directory, pattern="*.txt", full.names=TRUE, recursive=recursive)
  # set up the empty matrix
  wantedmat <- matrix(0,length(files),(n_cols))
  colnames(wantedmat) <- colNames
  # get record data
  for(i in 1:length(files)){
    fileConnection <- file(normalizePath(files[i]))
    record <- readLines(fileConnection)
    close(fileConnection)
    for (j in 1:length(colNames)) {
      field_code <- field_codes[j]
      MBEntry <- getMBRecordEntry(field_code, record, trimEntry=TRUE)
      if (length(MBEntry)>1) {
        MBEntry <- paste0(MBEntry, collapse="|")
      } else if (length(MBEntry)<1) {
        MBEntry <- NA
      }
      wantedmat[i,j] <- MBEntry
    }
  }
  write.csv(wantedmat,csv_name, row.names=F)
  return("Successfully wrote the csv")
}



#' Get A Specific Entry from MassBank Record
#'
#' @description This function returns the entries in the MassBank record matching
#' the \code{field_code}, or \code{NULL} if absent. The record should be in vector
#' format, e.g. with \code{\link{MBFileToVector}}. Used to retrieve information in
#' \code{\link{getMBRecordInfo}}.
#'
#' @usage getMBRecordEntry(field_code, record, trimEntry=TRUE)
#'
#' @param field_code The string to match to the MassBank record.
#' @param record The MassBank record in vector format to search, converted with
#' \code{\link{MBFileToVector}}.
#' @param trimEntry If \code{TRUE}, the \code{field_code} is removed, if \code{FALSE},
#' the full line is returned.
#'
#' @details Use \code{\link{getMBRecordPeaks}}, \code{\link{getMBRecordPeakAnnotations}}
#' and/or \code{\link{getMBPeaksAnnos}} to retrieve peaks and peak annotations.
#'
#' @return Returns a vector containing the matching entries or \code{NULL} if the retrieval failed.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>.
#'
#' @seealso \code{\link{MBFileToVector}}, \code{\link{getMBRecordInfo}}.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' record <- MBFileToVector(MBrecord)
#' field_code <- "CH$NAME"
#' field_code <- "CH$LINK: CAS"
#' field_code <- "CH$LINK CS"
#' field_code <- "CH$LINK: PUBCHEM CID"
#' getMBRecordEntry(field_code, record)
#'
getMBRecordEntry <- function(field_code, record, trimEntry=TRUE) {
  MBEntry <- grep(field_code, record, value=TRUE, fixed=TRUE)
  containsColon <- grepl(":",field_code,fixed=TRUE)
  if (containsColon) {
    substringAdjust <- 2
  } else {
    substringAdjust <- 3
  }
  if (length(MBEntry)>0 && trimEntry) {
    MBEntry <- substring(MBEntry, (nchar(field_code)+substringAdjust))
    #MBEntry <- sub(field_code,"",MBEntry,fixed=TRUE)
  } else if (length(MBEntry)==0) {
    MBEntry <- NULL
  }
  return(MBEntry)
}


#' Read a MassBank Record (or any file) into a Vector
#'
#' @description This opens the MassBank record (text file) in the
#' path \code{MBrecord} and uses \code{readLines} to save the values
#' in a vector
#'
#' @usage MBFileToVector(MBrecord)
#'
#' @param MBrecord File name (and path) to parse.
#'
#' @return Returns a vector containing the lines in the record
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>, based on code from Erik Mueller.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' MBFileToVector(MBrecord)
#'
MBFileToVector <- function(MBrecord) {
  file.conn <- file(MBrecord)
  record <- readLines(file.conn)
  close(file.conn)
  return(record)
}

#' Get A Peak List from MassBank Record
#'
#' @description This function extracts the peak list from a MassBank record, returning
#' a variable containing the peak m/z, intensity and relative intensity entries. The
#' column names are either original, or adjusted to be more R-friendly. Use the function
#' \code{\link{getMBPeaksAnnos}} to get peak lists with annotations.
#'
#' @usage getMBRecordPeaks(MBrecord, fixColNames = FALSE)
#'
#' @param MBrecord File name (and path) to a valid MassBank record file.
#' @param fixColNames Default \code{FALSE} uses the column names in the MassBank record,
#' i.e. \code{m/z}, \code{int.} and \code{rel.int.}. Use \code{TRUE} to make this more
#' R-friendly. The function \code{\link{fixMzIntColNames}} is used to do this.
#'
#' @return Returns a data frame containing peak data, or \code{NULL} if the retrieval failed.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>, based on code from Erik Mueller.
#'
#' @seealso \code{\link{getMBPeaksAnnos}}, \code{\link{fixMzIntColNames}}.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' getMBRecordPeaks(MBrecord, fixColNames = FALSE)
#' getMBRecordPeaks(MBrecord, fixColNames = TRUE)
#'
getMBRecordPeaks <- function(MBrecord, fixColNames = FALSE) {
  #extract the lines from the record
  file.conn <- file(MBrecord)
  record <- readLines(file.conn)
  close(file.conn)
  # get the indexes we need
  PKStart <- grep("PK$PEAK:", record, fixed = TRUE) + 1
  endslash <- grep("//", record, fixed = TRUE)
  n_peaks <- endslash-PKStart
  peak_colnames <- strsplit(grep('PK$PEAK:',record, value = TRUE, fixed = TRUE),split=" ")[[1]]
  # remove first entry as it's just PK:PEAK
  peak_colnames <- peak_colnames[-1]
  # replace "." and "/" and "i"=> "I"
  if (fixColNames) {
    peak_colnames <- fixMzIntColNames(peak_colnames)
  }
  if (n_peaks>0) {
    splitted <- strsplit(record[PKStart:(endslash - 1)], " ")
    PKPeak <- matrix(nrow = endslash - PKStart, ncol = 3)
    for (k in 1:length(splitted)) {
      splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
      PKPeak[k, ] <- splitted[[k]]
    }
    PKPeak <- as.data.frame(PKPeak, stringsAsFactors = FALSE)
    PKPeak[] <- lapply(PKPeak, type.convert)
    colnames(PKPeak) <- peak_colnames
  } else {
    PKPeak <- NULL
  }
  return(PKPeak)
}


#' Retrieve Annotations from a MassBank Record
#'
#' @description This function extracts the PK$ANNOTATION entries from a MassBank record,
#' returning a data frame containing the corresponding entries, with column headers taken
#' from the record. Use the function
#' \code{\link{getMBPeaksAnnos}} to get merged peak lists and annotations.
#'
#' @usage getMBRecordPeakAnnotations(MBrecord)
#'
#' @param MBrecord File name (and path) to a valid MassBank record file.
#'
#' @return Returns a data frame containing annotations, or \code{NULL} if the retrieval failed.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>, based on code from Erik Mueller.
#'
#' @seealso \code{\link{getMBPeaksAnnos}}, \code{\link{getMBRecordPeaks}}.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' getMBRecordPeakAnnotations(MBrecord)
#'
getMBRecordPeakAnnotations <- function(MBrecord) {
  #extract the lines from the record
  file.conn <- file(MBrecord)
  record <- readLines(file.conn)
  close(file.conn)
  # get the indexes we need
  PKannotationStart <- grep("PK$ANNOTATION:", record, fixed = TRUE) + 1
  ReadAnnotation <- TRUE # for unknowns, we'll need to test this, they won't be annotated
  ReadAnnotation <- ifelse(length(PKannotationStart)<1, FALSE, TRUE)
  PKNumLine <- grep("PK$NUM_PEAK:", record, fixed = TRUE)
  n_anno_entries <- PKNumLine - PKannotationStart
  # get the column names for the annotated entries
  anno_names <- strsplit(grep('PK$ANNOTATION:',record, value = TRUE, fixed = TRUE),split=" ")[[1]]
  # remove the first entry, which is the field code
  anno_names <- anno_names[-1]
  n_anno_cols <- length(anno_names)
  # now retreive the annotations
  if (ReadAnnotation && n_anno_entries >= 1) {
    splitted <- strsplit(record[PKannotationStart:(PKNumLine - 1)], " ")
    PKannotation <- matrix(nrow = n_anno_entries, ncol = n_anno_cols)
    for (k in 1:length(splitted)) {
      splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
      PKannotation[k, ] <- splitted[[k]]
    }
    PKannotation <- as.data.frame(PKannotation, stringsAsFactors = FALSE)
    PKannotation[] <- lapply(PKannotation, type.convert)
    colnames(PKannotation) <- anno_names
  } else {
    PKannotation <- NULL
  }
  return(PKannotation)
}


#' Retrieve a Merged Peak List with Annotations from a MassBank Record
#'
#' @description This function extracts the PK$PEAK and PK$ANNOTATION entries from a MassBank record,
#' returning a merged data frame containing the corresponding entries. This uses the
#' \code{\link{getMBRecordPeakAnnotations}} and \code{\link{getMBRecordPeaks}} functions.
#'
#' @usage getMBPeaksAnnos(MBrecord, fixColNames=FALSE)
#'
#' @param MBrecord File name (and path) to a valid MassBank record file.
#' @param fixColNames Default \code{FALSE} uses the column names in the MassBank record,
#' i.e. \code{m/z}, \code{int.} and \code{rel.int.}. Use \code{TRUE} to make this more
#' R-friendly. The function \code{\link{fixMzIntColNames}} is used to do this. NOTE: if
#' \code{TRUE}, this is only performed after merging to ensure the merge works.
#'
#' @details NOTE: for this to work properly, the m/z column names must match exactly.
#' This has only been verified so far on records generated with \code{RMassBank}.
#'
#' @return Returns a data frame containing annotated peaks, or \code{NULL} if the retrieval failed.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getMBRecordPeakAnnotations}}, \code{\link{getMBRecordPeaks}},
#' \code{\link{fixMzIntColNames}}.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' getMBPeaksAnnos(MBrecord)
#' getMBPeaksAnnos(MBrecord, fixColNames=TRUE)
#'
getMBPeaksAnnos <- function(MBrecord, fixColNames=FALSE) {
  # get the peaks (fixColNames=FALSE to allow merging)
  PKPeak <- getMBRecordPeaks(MBrecord, fixColNames=FALSE)
  PKAnnotation <- getMBRecordPeakAnnotations(MBrecord)
  bothExist <- !is.null(PKPeak) && !is.null(PKAnnotation)
  if (bothExist) {
    mergedPeaks <- merge(PKPeak,PKAnnotation)
    mergedColNames <- colnames(mergedPeaks)
    if (fixColNames) {
      mergedColNames <- fixMzIntColNames(mergedColNames)
      colnames(mergedPeaks) <- mergedColNames
    }
  } else if (!is.null(PKPeak)) {
    mergedPeaks <- PKPeak
    mergedColNames <- colnames(mergedPeaks)
    if (fixColNames) {
      mergedColNames <- fixMzIntColNames(mergedColNames)
      colnames(mergedPeaks) <- mergedColNames
    }
  } else {
    mergedPeaks <- NULL
  }
  return(mergedPeaks)
}

#' Fix MassBank Column Headers to be more R-friendly
#'
#' @description This function renames column names internally in \code{\link{getMBPeaksAnnos}},
#' and \code{\link{getMBRecordPeaks}} to make them more R-friendly.
#'
#' @usage fixMzIntColNames(ColNames)
#'
#' @param ColNames Column names to edit
#'
#' @details This function changes \code{m/z} to \code{mz}, replaces a decimal place with nothing
#' and changes \code{int} to \code{Int}.
#'
#' @return Returns edited column names.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getMBPeaksAnnos}}, \code{\link{getMBRecordPeaks}}.
#'
#' @export
#'
#' @examples
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' MBcolNames <- colnames(getMBRecordPeaks(MBrecord))
#' fixMzIntColNames(MBcolNames)
#'
fixMzIntColNames <- function(ColNames) {
  ColNames <- gsub(".","",ColNames,fixed=TRUE)
  ColNames <- gsub("/","",ColNames,fixed=TRUE)
  ColNames <- gsub("int","Int",ColNames,fixed=TRUE)
  return(ColNames)
}
