# Suspect Screening Functions
# Compiled from various scripts, also from Jen/Michele/Martin
# See acknowledgements per function
# E. Schymanski, 15/2/2017
# ported to ReSOLUTION package 3/4/17


# #Exact mass screening - modified code from Jen/Michele
# # copied from 201510022 suspect screening with R.R then modified
# ppm2 <- function(mz, ppmlimit) (1e-6 * mz * ppmlimit)
# ppmMax <- 5
# This belongs in ReSOLUTION
#

#' Screen Suspects in Peak List by Mass
#'
#' @description This screens a peak list for given suspects by exact mass
#'
#' @usage screenSuspects(peaklist, suspects,ppmMax=5)
#'
#' @param peaklist A peak list containing columns \code{mz}, \code{RT}, \code{Int}.
#' @param suspects A list of suspects containing columns \code{Name} and \code{Suspect_mass}.
#' @param ppmMax Maximum \code{ppm} to match masses.
#'
#' @return Returns an ordered list of suspect hits.
#'
#' @author Jennifer Schollee, Michael Stravs, Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @export
#'
#' @examples
#'
screenSuspects <- function(peaklist, suspects,ppmMax=5) {
  #works on peaklist with column names mz, RT, Int
  # and suspect lists with columns "Name" and "Suspect_mass"
  ppm2 <- function(mz, ppmlimit) (1e-6 * mz * ppmlimit)
  #ppmMax <- 5
  suspect_hits <- do.call(rbind, (lapply(1:nrow(peaklist), function(n) {
    mass <- peaklist[n,"mz"]
    peaklist.hit <- suspects[abs(suspects$Suspect_mass - mass) < ppm2(mass, ppmMax),,drop=FALSE]

    peaklist.hit$dppm <- ((peaklist.hit$Suspect_mass / mass)-1)*1e6
    peaklist.hit$mz <- rep(peaklist[n, "mz"], nrow(peaklist.hit))
    peaklist.hit$RT <- rep(peaklist[n, "RT"], nrow(peaklist.hit))
    peaklist.hit$Int <- rep(peaklist[n, "Int"], nrow(peaklist.hit))
    #peaklist.hit$ID <- rep(peaklist[n, "ID"], nrow(peaklist.hit))
    #peaklist.hit$linkID <- rep(peaklist[n,"link ID"], nrow(peaklist.hit))
    #suspect_hits <- do.call(rbind, peaklist.hit)
    #ordered_hits <- suspect_hits[order(suspect_hits$Name, decreasing = FALSE), ]
    #peaklist.hits <- do.call(rbind, peaklist.hit)
    return(peaklist.hit)
  })))
  ordered_hits <- suspect_hits[order(suspect_hits$Name, decreasing = FALSE), ]
  return(ordered_hits)
}


# build.homol <- function(start_formula, building_block, series_name, n_start, n_max, rt=15) {
#   formulas <- c(start_formula, building_block)
#   checked <- check_chemform(isotopes,formulas)
#   member_formulas <- checked[1,2]
#   n <- 1
#   new_formula <- ""
#   neutral_monoiso_mass <- as.numeric(findMz.formula(member_formulas[1],"")[3])
#   #series_MmHm <- as.numeric(findMz.formula(series_formulas[1],"mH")[3])
#   member_names <- paste(series_name,as.character(n_start),sep="")
#   rt <- rt # in minutes
#   for(i in 2:(n_max-n_start)) {
#     #n <- n+1
#     formulas[1] <- member_formulas[n]
#     checked <- check_chemform(isotopes,formulas)
#     new_formula <- mergeform(checked[1,2],checked[2,2])
#     member_formulas[i] <- new_formula
#     neutral_monoiso_mass[i] <- as.numeric(findMz.formula(new_formula,"")[3])
#     #series_MmHm[i] <- as.numeric(findMz.formula(new_formula,"mH")[3])
#     member_names[i] <- paste(series_name,as.character(n_start+i-1),sep="")
#     rt[i] <- rt[1]
#     n <- n+1
#   }
#   homol_series <- cbind(member_names, member_formulas, rt, neutral_monoiso_mass)
#   return(homol_series)
# }

# Screen homologue series, written by ES, relies on series built in RChemMass with build.homol


#' Screen EICs for Homologue Series
#'
#' @description This is a wrapper for plotEICs designed for screening homologue series
#' built from \code{\link{build.homol}}.
#'
#' @usage screen.homol(file, file_desc, series, series_name, adduct, rt_window=NULL, ppm=10,
#' greyscale=FALSE, smiles="", formula=NULL, addMSMS=FALSE, isLog=FALSE,
#' isPos=TRUE, kekulise=TRUE, writePDF=FALSE, writeCSV=FALSE)
#'
#' @param file Path and filename to mz(X)ML file to screen
#' @param file_desc Description for the file (used in plot title and file names)
#' @param series The series from \code{\link{build.homol}} to screen
#' @param series_name The name of the series (used in legend)
#' @param adduct The adduct species of the homologue to screen. NOTE: the adduct must be exactly as written 
#' in the adducts table in the enviPat package.
#' @param rt_window The retention time window (in seconds) for EIC extraction. This is added
#' to \code{rt} to form the extraction window. Should not be zero! If \code{NULL},
#' the whole EIC is extracted.
#' @param ppm The ppm limit for EIC extraction
#' @param greyscale If \code{TRUE}, greyscale plot will be produced. If \code{FALSE},
#' rainbow colours are assigned automatically.
#' @param smiles SMILES code of the compound of interest. Can also be generic SMILES
#' (this needs to be tested...).
#' @param formula A text field to put the formula of the series if no SMILES.
#' @param addMSMS Indicates whether MSMS scans should be added to plot with default settings \code{TRUE}.
#' \code{FALSE} does not attempt to add MSMS scans.
#' @param isLog If \code{TRUE}, EICs are plotted with a log scale; \code{FALSE} with
#' absolute intensities. Zero intensities are lost with \code{TRUE}. TRUE may not work properly.
#' @param isPos Setting to pass on to \code{\link{plotEICs}}. \code{TRUE} indicates positive mode;
#' \code{FALSE} negative mode. Only relevant for switching data.
#' @param kekulise Controls aromaticity detection in SMILES for plotting.
#' @param writePDF Writes the plot to PDF rather than active device. PDF name is autocreated.
#' @param writeCSV Writes the screen summary to CSV - file name is autocreated.
#'
#' @return Returns a plot in the active device, or PDF and summary output.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{plotEICs}}, \code{\link{build.homol}}
#'
#' @export
#'
#' @examples
#' pegs <- build.homol("C4H10O3","C2H4O","PEGS",2,14)
#' library(RMassBank)
#' library(RMassBankData)
#' RmbDefaultSettings()
#' mzML_files <- list.files(system.file("spectra", package="RMassBankData"), ".mzML", full.names = TRUE)
#' mzML_file <- mzML_files[1]
#' screen.homol(mzML_file, "test_file", pegs, "PEGS_MpH", "M+H")
#' screen.homol(mzML_file, "test_file", pegs, "PEGS_MpNH4", "M+NH4")
#' #Note this doesn't look pretty as PEGS aren't present, but does test that it works.
#'
#'
screen.homol <- function(file, file_desc, series, series_name, adduct, rt_window=NULL, ppm=10,
                         greyscale=FALSE, smiles="", formula=NULL, addMSMS=FALSE, isLog=FALSE,
                         isPos=TRUE, kekulise=TRUE, writePDF=FALSE, writeCSV=FALSE) {
  #first, calculate the adduct masses from the series output from build.homol
  # NOTE: adduct must be exactly as written in Martin's enviPat package.
  data(isotopes, package="enviPat")
  data(adducts, package="enviPat")
  adduct_index <- which(as.character(adducts$Name)==adduct)
  if (length(adduct_index)<1) {
    stop("Incorrect adduct specification, adduct must match entry in enviPat")
  }
  adduct_name <- sub("+","p",adduct,fixed=TRUE)
  adduct_name <- sub("-","m",adduct_name,fixed=TRUE)
  #rt_window <- rt_window
  #formula <- formula
  mz <- as.numeric(series[,4]) + adducts[adduct_index,5]
  rt <- 60*as.numeric(series[,3])
#  mzrt <- cbind(mz,rt)
  formulas <- series[,2]
  names <- series[,1]
  #set up titles
  plot_title <- paste(series_name, adduct, file_desc)
  pdf_title <- paste(series_name, "_",adduct_name, "_",file_desc, ".pdf",sep="")
  csv_name <- paste(series_name, "_",adduct_name, "_",file_desc, ".csv",sep="")

  if (writePDF) {
    pdf(pdf_title)
  }

  #extract the EICs
  screen_summary <- plotEICs(mzML_filePath=file, mz=mz, rt=rt, rt_in_sec = TRUE, names=names,
                             smiles=smiles, plot_title=plot_title, rt_window=rt_window, ppm=10,
                             takePPM=TRUE, isPos=isPos, isLog=isLog, greyscale=greyscale,
                             kekulise=kekulise, addMSMS=addMSMS)
  if (!is.null(formula)) {
    mtext(formula, side=3, line=0)
  }
  if (writePDF) {
    dev.off()
  }

  if (writeCSV) {
    write.csv(screen_summary, file=csv_name,row.names=FALSE)
  }
  #

  # RmbDefaultSettings()
  # f <- openMSfile(file)
  # h <- header(f)
  # c <- makePeaksCache(f,h)
  # if (!is.null(rt_window)) {
  #   eics <- apply(mzrt, 1,
  #                 function(row)
  #                 {
  #                   findEIC.cached(f,row[[1]],ppm(row[[1]],ppm,p=TRUE), rtLimit = rt_window,
  #                                  headerCache=h, peaksCache=c)
  #                 })
  # } else {
  #   eics <- apply(mzrt, 1,
  #                 function(row)
  #                 {
  #                   findEIC.cached(f,row[[1]],ppm(row[[1]],ppm,p=TRUE), headerCache=h, peaksCache=c)
  #                 })
  # }
  # #rm(c)
  # #gc()
  # colors <- rainbow(length(eics), start=2/6, end=1)
  # if (greyscale) {
  #   colors <- gray((1:(length(eics)+1))/(length(eics)+1))
  #   # to flip the colour (starting light), pick the below instead
  #   #    colors <- gray(1-(1:(length(eics)+1))/(length(eics)+1))
  # }
  # #retrieve the summary data of the EICs
  # maxmin_rt <- range(unlist(lapply(eics, function(eic) {range(eic$rt)})), finite = T)
  # max_int <- max(unlist(lapply(eics, function(eic) {max(eic$intensity)})))
  # rt_at_maxI <- unlist(lapply(eics, function(eic) {eic[which.max(eic$intensity),1]}))
  # RT_minutes <- as.numeric(rt_at_maxI)/60
  # I_at_maxI <- unlist(lapply(eics, function(eic) {eic[which.max(eic$intensity),2]}))
  # screen_summary <- cbind(mz, rt, names, formulas, RT_minutes, I_at_maxI)
  # write outputs to CSV and PDF
  # write.csv(screen_summary, file=csv_name)
  # pdf(pdf_title)
  # plot.new()
  # if (!is.null(rt_window)) {
  #   x_coords <- rt_window
  # } else {
  #   x_coords <- maxmin_rt
  # }
  #
  # plot.window(x_coords, c(0,1.4*max(I_at_maxI)))
  # struct_coords <- c(min(x_coords), max(I_at_maxI), min(x_coords)+0.4*(max(x_coords)-min(x_coords)),
  #                    1.4*max(I_at_maxI))
  # #plot.window(rt_window, c(0,max(I_at_maxI)))
  # box()
  # for(n in 1:length(eics))
  # {
  #   if (!is.null(smiles)) {
  #     renderSMILES.parse2(smiles, kekulise = FALSE, coords = struct_coords)
  #   }
  #   if (logInt) {
  #     lines(intensity ~ rt ,data = eics[[n]], col=colors[[n]],lwd=1.75,ylog=TRUE)
  #   } else {
  #     lines(intensity ~ rt ,data = eics[[n]], col=colors[[n]],lwd=1.75)
  #   }
  #   title(main=plot_title, xlab="RT (sec)", ylab="Intensity")
  #   legend("topright", legend=names, col=colors, lty=1, lwd=3)
  #   axis(1, at=seq(0,maxmin_rt[2], by=60))
  #   #axis(1, at=seq(0,rt_window[2], by=60))
  #   axis(2)
  # }
  # dev.off()
  return(screen_summary)
}

# pegs <- build.homol("C4H10O3","C2H4O","PEGS",2,14)

#   screen.homol(mzXML_files[i],short_locations[i],pegs,"PEGS_MpH","M+H")
#   screen.homol(mzXML_files[i],short_locations[i],pegs,"PEGS_MpNH4","M+NH4")
#   screen.homol(mzXML_files[i],short_locations[i],pegs,"PEGS_MpNa","M+Na")
#   screen.homol(mzXML_files[i],short_locations[i],pegs,"PEGS_MpK","M+K")

# spacs <- build.homol("C9H10O5S","CH2","SPAC",3,15)
# las <- build.homol("C16H26O3S","CH2","LAS",10,14)
# spadcs <- build.homol("C9H8O7S","CH2","SPADC",1,15)
#   screen.homol(mzXML_files[i],short_locations[i],las,"LAS","M-H")
#   screen.homol(mzXML_files[i],short_locations[i],spacs,"SPACS","M-H")
#   screen.homol(mzXML_files[i],short_locations[i],spadcs,"SPADCS","M-H")

