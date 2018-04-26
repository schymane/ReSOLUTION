# Various code snippets for plotting MSMS spectra
# E. Schymanski, 11/1/2016, copied from MetFragResultPlotting.R
# Developed further in the mean time ..
# 31/3/2017 ported to ReSOLUTION


#' Plot Single MS/MS Spectrum
#'
#' @description This is a wrapper function for \code{\link{plotSpectra}} for backwards
#' compatability with previous workflows/scripts. \code{\link{plotSpectra}} offers
#' much more functionality.
#'
#' @usage plotSpectrum(mz, int, main, smiles="", labels=NULL, formatFormula=TRUE,
#' absInt=FALSE,ylim_factor=2.5,max_xlim=0,kekulise=TRUE,color="black")
#'
#' @param mz Vector containing mz values (x) to plot.
#' @param int Vector containing intensity values (y) to plot.
#' @param main Title for the plot.
#' @param smiles SMILES code of the structure to plot. Leave empty for no structure.
#' @param labels Vector containing labels for selected peaks. If \code{NULL}, no labels
#' are plotted.
#' @param formatFormula If \code{TRUE}, uses \code{chemistry2expression} to create subscript numbers.
#' Any \code{+} or \code{-} are removed to avoid errors. If \code{FALSE}, \code{labels} are printed as is.
#' This is passed to \code{chem2express} in \code{\link{plotSpectra}}.
#' @param absInt If \code{TRUE}, absolute intensity values are used. If \code{FALSE}, relative intensities
#' are calculated on individual spectra, scaled to 1000.
#' @param ylim_factor Scaling factor for y-axis. Default \code{2.5} leaves ample space for structure,
#' annotations and legend. Reduce for peaks only. Passed to \code{yfactor} in \code{\link{plotSpectra}}.
#' @param max_xlim Option to control the maximum value on the x axis. Passed to \code{xlim} in
#' \code{\link{plotSpectra}} as \code{c(0,max_xlim)}
#' @param kekulise Controls aromaticity detection of \code{smiles}. Default \code{TRUE}
#' should work best with rcdk>3.4.1. Try \code{FALSE} if rings and double bonds appear together
#' in aromatic rings.
#' @param color Optional vector to re-define default line colour. Passed to \code{line_colour} in
#' \code{\link{plotSpectra}}.
#'
#' @return Returns a plot in the current plotting device.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{plotSpectra}}, \code{\link{renderSMILES.rcdk}}, \code{\link{chemistry2expression}},
#' \code{\link{trimAnnotation}}.
#'
#' @export
#'
#' @examples
#' mz_1 <- c(58.0287, 111.0441, 168.0655, 210.1125)
#' int_1 <- c(23.0000, 999.0000, 843.5855, 999.0000)
#' labels_1 <- c("C2H4NO", "C6H7O2", "C8H10NO3", "C11H16NO3")
#' smiles_1 <- "O=C(Oc1ccccc1OC(C)C)NC"
#' plotSpectrum(mz_1, int_1, main="test spec",labels=labels_1, smiles=smiles_1)
#'
plotSpectrum <- function(mz, int, main, smiles="", labels=NULL, formatFormula=TRUE,
                         absInt=FALSE,ylim_factor=2.5,max_xlim=0,kekulise=TRUE,color="black") {
  # set xlim values to pass on to plotSpectra
  if (max_xlim == 0) {
    max_xlim <- max(mz)*1.2
  }
  plotSpectra(mz=mz, int=int, main=main, smiles=smiles, labels=labels,labelSpec=1,
              chem2express = formatFormula, absInt=absInt, yfactor=ylim_factor,
              xlim=c(0,max_xlim), kekulise=kekulise, line_colour=color, legend=main)
}

#' Plot Single Annotated MS/MS Spectrum
#'
#' @description This is a wrapper function for \code{\link{plotSpectra}} for backwards
#' compatability with previous workflows/scripts. \code{\link{plotSpectra}} offers
#' much more functionality.
#'
#' @usage plotAnnoSpectrum(mz, int, main, smiles="", labels=NULL,mz_a=NA,int_a=NA, formatFormula=TRUE,
#' absInt=FALSE,ylim_factor=2.5,max_xlim=0,kekulise=TRUE)
#'
#' @param mz Vector containing mz values (x) to plot.
#' @param int Vector containing intensity values (y) to plot.
#' @param main Title for the plot.
#' @param smiles SMILES code of the structure to plot. Leave empty for no structure.
#' @param labels Vector containing labels for selected peaks. If \code{NULL}, no labels
#' are plotted.
#' @param mz_a Vector containing mz values for second spectrum.
#' @param int_a Vector containing intensity values for second spectrum.
#' @param formatFormula If \code{TRUE}, uses \code{chemistry2expression} to create subscript numbers.
#' Any \code{+} or \code{-} are removed to avoid errors. If \code{FALSE}, \code{labels} are printed as is.
#' This is passed to \code{chem2express} in \code{\link{plotSpectra}}.
#' @param absInt If \code{TRUE}, absolute intensity values are used. If \code{FALSE}, relative intensities
#' are calculated on individual spectra, scaled to 1000.
#' @param ylim_factor Scaling factor for y-axis. Default \code{2.5} leaves ample space for structure,
#' annotations and legend. Reduce for peaks only. Passed to \code{yfactor} in \code{\link{plotSpectra}}.
#' @param max_xlim Option to control the maximum value on the x axis. Passed to \code{xlim} in
#' \code{\link{plotSpectra}} as \code{c(0,max_xlim)}
#' @param kekulise Controls aromaticity detection of \code{smiles}. Default \code{TRUE}
#' should work best with rcdk>3.4.1. Try \code{FALSE} if rings and double bonds appear together
#' in aromatic rings.
#' @param legend Vector containing legend entries for the two spectra.
#'
#' @return Returns a plot in the current plotting device.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{plotSpectra}}, \code{\link{renderSMILES.rcdk}}, \code{\link{chemistry2expression}},
#' \code{\link{trimAnnotation}}.
#'
#' @export
#'
#' @examples
#' mz_1 <- c(58.0287, 111.0441, 168.0655, 210.1125)
#' int_1 <- c(23.0000, 999.0000, 843.5855, 999.0000)
#' labels_1 <- c("C2H4NO", "C6H7O2", "C8H10NO3", "C11H16NO3")
#' smiles_1 <- "O=C(Oc1ccccc1OC(C)C)NC"
#' plotAnnoSpectrum(mz_1, int_1, main="test spec",labels=labels_1, smiles=smiles_1, mz_a = mz_1,
#' int_a = int_1)
#'
plotAnnoSpectrum <- function(mz, int, main, smiles="", labels=NULL, formatFormula=TRUE,mz_a=NA,int_a=NA,
                          absInt=FALSE,ylim_factor=2.5,max_xlim=0,kekulise=TRUE, legend=c("spec1","spec2")) {
  # set xlim values to pass on to plotSpectra
  if (max_xlim == 0) {
    max_xlim <- max(mz)*1.2
  }
  plotSpectra(mz=mz, int=int, main=main, smiles=smiles, labels=labels,labelSpec=1,
              chem2express = formatFormula, absInt=absInt, yfactor=ylim_factor,mz_a = mz_a, int_a=int_a,
              xlim=c(0,max_xlim), kekulise=kekulise, legend=legend)
}



#' Plot Multiple Spectra on a Single Plot
#'
#' @description This is the main function to plot multiple mass spectra
#' on one plot, with many options to control the output.
#'
#' @usage plotSpectra(mz, int, main, labels=NULL,labelSpec=1,label_mz=NA,
#' label_int=NA,mz_a=NA,int_a=NA, absInt=FALSE, logInt=FALSE, mz_b=NA,int_b=NA,
#' smiles="",kekulise=TRUE,legend=c("spec1","spec2","spec3"),
#' xlim=NULL,ylim=NULL, yfactor=2.5, line_colour=c("black","red","green"),
#' line_type=c("solid","dashed","dotdash"),line_width=c(2,2,1),
#' chem2express=FALSE, srt=60, pos=4,offset=-0.2)
#'
#' @param mz Vector containing mz values (x) to plot first.
#' @param int Vector containing intensity values (y) to plot first.
#' @param main Title for the plot.
#' @param labels Vector containing labels for selected peaks. If \code{NULL}, no labels
#' are plotted.
#' @param labelSpec Indicates which spectrum to use for label coordinates. Must be same
#' length as \code{labels} to avoid recycling labels. If \code{1}, \code{mz} and \code{int}
#' are used, if \code{2} then \code{mz_a} and \code{int_a}, if \code{3}, then \code{mz_b}
#' and \code{int_b}. Parameters \code{label_mz} and \code{label_int} take precedence.
#' @param label_mz Values defining starting x coordinates for \code{labels}.
#' @param label_int Values defining starting y coordinates for \code{labels}.
#' @param mz_a Vector containing mz values for second spectrum.
#' @param int_a Vector containing intensity values for second spectrum.
#' @param absInt If \code{TRUE}, absolute intensity values are used. Be careful that all
#' spectra match! If \code{FALSE}, relative intensities are calculated on individual
#' spectra, scaled to 1000.
#' @param logInt If \code{FALSE}, intensity values are used as is. If \code{TRUE}, log
#' values are used; this can be useful to find small peaks.
#' @param mz_b Vector containing mz values for third spectrum.
#' @param int_b Vector containing intensity values for third spectrum.
#' @param smiles SMILES code of the structure to plot. Leave empty for no structure.
#' @param kekulise Controls aromaticity detection of \code{smiles}. Default \code{TRUE}
#' should work best with rcdk>3.4.1. Try \code{FALSE} if rings and double bonds appear together
#' in aromatic rings.
#' @param legend Vector containing legend entries.
#' @param xlim Option to control x-axis limits. If \code{NULL}, \code{xlim=c(0,1.1*max(all_mz))}
#' @param ylim Option to control y-axis limits. If \code{NULL}, \code{ylim=c(0,yfactor*max(all_int))}
#' unless log scale is used, where \code{ylim=c(-1, yfactor+8)}.
#' @param yfactor Scaling factor for y-axis. Default \code{2.5} leaves ample space for structure,
#' annotations and legend. Reduce for peaks only.
#' @param line_colour Optional vector to re-define default line colours.
#' @param line_type Optional vector to re-define default line type.
#' @param line_width Optional vector to re-define default line widths.
#' @param chem2express If \code{TRUE}, uses \code{chemistry2expression} to create subscript numbers.
#' Any \code{+} or \code{-} are removed to avoid errors. If \code{FALSE}, \code{labels} are printed as is.
#' @param srt Optional parameter to adjust \code{srt} value for \code{labels}.
#' @param pos Optional parameter to adjust \code{pos} value for \code{labels}.
#' @param offset Optional parameter to adjust \code{offset} value for \code{labels}.
#'
#' @return Returns a plot in the current plotting device.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @details
#' For controlling label orientation, the following settings are useful:
#' Default (angle): \code{srt=60}, \code{pos=4}, \code{offset=-0.2};
#' Vertical above peak: \code{srt=90}, \code{pos=4}, \code{offset=0};
#' Horizontal, centred above peak: \code{srt=0}, \code{pos=3}, \code{offset=0.2}.
#'
#' To avoid overlapping labels, use \code{\link{trimAnnotation}}.
#'
#'
#' @seealso \code{\link{plot}}, \code{\link{renderSMILES.rcdk}}, \code{\link{chemistry2expression}},
#' \code{\link{text}}, \code{\link{trimAnnotation}}.
#'
#' @export
#'
#' @examples
#' mz_1 <- c(58.0287, 111.0441, 168.0655, 210.1125)
#' int_1 <- c(23.0000, 999.0000, 843.5855, 999.0000)
#' labels_1 <- c("C2H4NO", "C6H7O2", "C8H10NO3", "C11H16NO3")
#' smiles_1 <- "O=C(Oc1ccccc1OC(C)C)NC"
#' plotSpectra(mz_1, int_1, main="test spec",labels=labels_1, labelSpec=1, smiles=smiles_1)
#' plotSpectra(mz_1, int_1, main="test spec",labels=labels_1, labelSpec=1, smiles=smiles_1,
#' chem2express=TRUE, mz_a=mz_1, int_a = int_1, legend=c("spec1","spec1 again"))
#'
plotSpectra <- function(mz, int, main, labels=NULL,labelSpec=1,label_mz=NA, label_int=NA,
                        mz_a=NA,int_a=NA, absInt=FALSE, logInt=FALSE, mz_b=NA,int_b=NA,
                        smiles="",kekulise=TRUE,legend=c("spec1","spec2","spec3"),
                        xlim=NULL,ylim=NULL, yfactor=2.5, line_colour=c("black","red","green"),
                        line_type=c("solid","dashed","dotdash"),line_width=c(2,2,1),
                        chem2express=FALSE, srt=60, pos=4,offset=-0.2) {
  # set the limits to data limits
  xrange <- c(min(mz,mz_a,mz_b,na.rm = TRUE), max(mz,mz_a,mz_b,na.rm = TRUE))
  yrange <- c(min(int,int_a,int_b,na.rm = TRUE), max(int,int_a,int_b,na.rm = TRUE))
  if (yrange[1]>0) {
    yrange[1] <- 0
  }
  if (logInt) {
    yrange <- log10(yrange)
  }
  if (is.null(xlim)) {
    xlim <- c(0, 1.1*xrange[2])
  }
  if (is.null(ylim)) {
    if (absInt) {
      ylim <- c(yrange[1], yfactor*yrange[2])
      if (logInt) {
        ylim <- c(-2, yfactor+yrange[2])
      }
    } else {
      ylim <- c(yrange[1], yfactor*1000)
      if (logInt) {
        ylim <- c(-1, yfactor+8)
      }
    }
  }
  #normalise the peak lists before plotting!!!
  if (!absInt) {
    int <- 1000*int/max(int)
  }
  if (logInt) {
    plot(mz, log10(int), type="h", xlim=xlim, ylim=ylim, main=main,
         cex=1, lwd=line_width[1], col=line_colour[1], xlab="m/z", ylab="log10(Intensity)")
  } else {
    plot(mz, int, type="h", xlim=xlim, ylim=ylim, main=main,
         cex=1, lwd=line_width[1], col=line_colour[1], xlab="m/z", ylab="Intensity")
  }
  n_spec=1
  # next spectrum
  if (!is.na(mz_a[1])) {
    if (!absInt) {
      int_a <- 1000*int_a/max(int_a)
    }
    if (logInt) {
      int_a <- log10(int_a)
    }
    lines(mz_a, int_a, type="h",col= line_colour[2], lty=line_type[2],lwd=line_width[2])
    n_spec <- n_spec+1
  }
  if (!is.na(mz_b[1])) {
    if (!absInt) {
      int_b <- 1000*int_b/max(int_b)
    }
    if (logInt) {
      int_b <- log10(int_b)
    }
    lines(mz_b, int_b, type="h",col= line_colour[3], lty=line_type[3],lwd=line_width[3])
    n_spec <- n_spec+1
  }
  legend("topright", legend=legend[1:n_spec],
         col=line_colour[1:n_spec], lty=line_type[1:n_spec], lwd=line_width[1:n_spec], bty="n")
  #plot the molecule onto the spectrum
  if (nchar(smiles)>1) {
    renderSMILES.rcdk(smiles, coords=c(xlim[1],ylim[2]/yfactor,xlim[1]+0.3*(xlim[2]-xlim[1]),ylim[2]))
  }
  # then plot the annotations
  if (!is.null(labels)) {
    if (chem2express) {
      labels <- gsub("+","",labels,fixed=TRUE)
      labels <- gsub("-","",labels,fixed=TRUE)
      labels <- chemistry2expression(labels)
    }
    if (!is.na(label_mz) && !is.na(label_int)) {
      if (logInt) {
        label_int <- log10(label_int)
      }
      text(label_mz, label_int, labels=labels,srt=srt,pos=pos,offset=offset)
    } else if (labelSpec==1) {
      if (logInt) {
        int <- log10(int)
      }
      text(mz, int, labels=labels,srt=srt,pos=pos,offset=offset)
    } else if (labelSpec==2) {
      text(mz_a, int_a, labels=labels,srt=srt,pos=pos,offset=offset)
    } else if (labelSpec==3) {
      text(mz_b, int_b, labels=labels,srt=srt,pos=pos,offset=offset)
    }
  }

  invisible(NULL)
}

#' Convert Molecular Formulas to Regular Expression for Labels
#'
#' @description Creates subscripts for numerals in molecular formulas for plotting
#'
#' @usage chemistry2expression(formulas)
#'
#' @param formulas Vector containing formulas (as strings) to convert
#'
#' @return Expression for plotting labels.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> and Steffen Neumann <sneumann@@ipb-halle.de>
#'
#' @seealso \code{\link{plotSpectra}}, \code{\link{trimAnnotation}}, \code{\link{text}}.
#'
#' @export
#'
#' @examples
#' labels_1 <- c("C2H4NO", "C6H7O2", "C8H10NO3", "C11H16NO3")
#' chemistry2expression(labels_1)
#'
chemistry2expression <- function(formulas) {
  exprs <- sub("\\*$", "", gsub("([0-9]+)", "[\\1]*", formulas, fixed=FALSE))
  parse(text=exprs)
}

#' Plot MS1 Spectrum
#'
#' @description This is a wrapper function for \code{\link{plotSpectra}} for backwards
#' compatability with previous workflows/scripts. \code{\link{plotSpectra}} offers
#' much more functionality. \code{xlim} and \code{ylim} are fixed suited to MS1 spectra.
#'
#' @usage plotMS(mz, int, main, labels=NULL, formatFormula=TRUE, absInt=FALSE,ylim_factor=2)
#'
#' @param mz Vector containing mz values (x) to plot.
#' @param int Vector containing intensity values (y) to plot.
#' @param main Title for the plot.
#' @param labels Vector containing labels for selected peaks. If \code{NULL}, no labels
#' are plotted.
#' @param formatFormula If \code{TRUE}, uses \code{chemistry2expression} to create subscript numbers.
#' Any \code{+} or \code{-} are removed to avoid errors. If \code{FALSE}, \code{labels} are printed as is.
#' This is passed to \code{chem2express} in \code{\link{plotSpectra}}.
#' @param absInt If \code{TRUE}, absolute intensity values are used. If \code{FALSE}, relative intensities
#' are calculated on individual spectra, scaled to 1000.
#' @param ylim_factor Scaling factor for y-axis. Default \code{2.5} leaves ample space for structure,
#' annotations and legend. Reduce for peaks only. Passed to \code{yfactor} in \code{\link{plotSpectra}}.
#'
#' @return Returns a plot in the current plotting device.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{plotSpectra}}, \code{\link{renderSMILES.rcdk}}, \code{\link{chemistry2expression}},
#' \code{\link{trimAnnotation}}.
#'
#' @export
#'
#' @examples
#' ms1_mz <- c(210.1125, 211.1154, 212.1187)
#' ms1_int <- c(420889280, 52790856, 2265951.5)
#' ms1_labels <- c("[M+H]+", "13C-M+1", "13C-M+2")
#' smiles_1 <- "O=C(Oc1ccccc1OC(C)C)NC"
#' plotMS(ms1_mz, ms1_int, main="test MS1 spec",labels=ms1_labels, formatFormula=FALSE,
#' absInt=TRUE)
#' plotMS(ms1_mz, ms1_int, main="test MS1 spec",labels=ms1_labels, formatFormula=FALSE,
#' absInt=FALSE, smiles=smiles_1)
#'
plotMS <- function(mz, int, main, labels=NULL, formatFormula=TRUE,
                         absInt=FALSE,ylim_factor=2, smiles="") {
  # set limits
  if (!absInt) {
    int <- 1000*int/max(int)
  }
  xlim=c((min(mz)-5),(max(mz)+5))
  ylim=c(0,(max(int)*ylim_factor))
  # plot
  plotSpectra(mz=mz, int=int, main=main, labels=labels,labelSpec=1, xlim=xlim, ylim=ylim,
              chem2express = formatFormula, absInt=absInt, yfactor=ylim_factor, legend=main,smiles=smiles)
}

#outdated function, replaced with renderSMILES.rcdk
#' Add Structure to Plot (outdated)
#'
#' @description This is an outdated function to add a structure to a plot,
#' kept for backwards compatability. Use \code{\link{renderSMILES.rcdk}} instead.
#'
#' @param smiles SMILES of structure to add to plot
#' @param kekulise Control aromaticity detection. Default \code{TRUE} is best with
#' recent \code{\link{rcdk}} versions.
#' @param coords The coordinates (xmin, ymin, xmax, ymax) for the structure.
#'
#' @return Adds the structure to the plot in the active device
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{renderSMILES.rcdk}}.
#'
#' @export
#'
#' @examples
#' ms1_mz <- c(210.1125, 211.1154, 212.1187)
#' ms1_int <- c(420889280, 52790856, 2265951.5)
#' ms1_labels <- c("[M+H]+", "13C-M+1", "13C-M+2")
#' plotMS(ms1_mz, ms1_int, main="test MS1 spec",labels=ms1_labels, formatFormula=FALSE,
#' absInt=FALSE)
#' addStructToPlot("O=C(Oc1ccccc1OC(C)C)NC", kekulise=TRUE, coords=c(205,1000,209,2000))
#' addStructToPlot("O=C(Oc1ccccc1OC(C)C)NC", kekulise=FALSE, coords=c(205,1000,209,2000))
#'
addStructToPlot <- function(smiles,kekulise=TRUE,coords=c(0,1000,100,2500)) {
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise=kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not plotted: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol, width = 150, height = 150))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    # img <- view.image.2d(mol, width = 150, height = 150)
    if (length(img)>2) {
      rasterImage(img, coords[1],coords[2],coords[3],coords[4])
    }
  }
}


#' Compare Mass Spectra (basic wrapper for OrgMassSpecR function)
#'
#' @description This uses the \code{\link{SpectrumSimilarity}} function of
#' \code{OrgMassSpecR} to calculate a similarity value, output the table to
#' a separate text file and a plot to the active device
#'
#' @usage SpecComp.OrgMassSpecR(spec1, spec2, output_file_name, t=0.005, b=0.1, header=FALSE)
#'
#' @param spec1 File containing \code{mz} and \code{intensity} values of first (top) spectrum.
#' @param spec2 File containing \code{mz} and \code{intensity} values of second (bottom) spectrum.
#' @param output_file_name File name for the output table. If empty no output is written.
#' @param t Numeric value specifying the tolerance used to align the m/z values of the two spectra.
#' @param b Numeric value specifying the baseline threshold for peak identification.
#' Expressed as a percent of the maximum intensity.
#' @param header Sets whether \code{spec1} and \code{spec2} contain header row (\code{TRUE} or not (\code{FALSE}).
#'
#' @return Returns a similarity value, a file containing the comparison table and a plot in the active device.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>,
#' Nathan Dodder (\code{\link{SpectrumSimilarity}} function).
#'
#' @seealso \code{\link{SpectrumSimilarity}}.
#'
#' @export
#'
#' @examples
#' spec1 <- system.file("extdata","EA000405_peaks.txt",package="ReSOLUTION")
#' spec2 <- system.file("extdata","EA000407_peaks.txt",package="ReSOLUTION")
#' sim <- SpecComp.OrgMassSpecR(spec1,spec2, header=F)
#' sim
#'
SpecComp.OrgMassSpecR <- function(spec1, spec2, output_file_name="", t=0.005, b=0.1, header=FALSE) {
  top_spec <- read.table(spec1, header=header)
  top_spec_name <- basename(spec1)
  bottom_spec <- read.table(spec2, header=header)
  bottom_spec_name <- basename(spec2)
  max.x <- max(max(top_spec[,1]), max(bottom_spec[,1]))
  if (nchar(output_file_name)>3) {
    capture.output(sim <- SpectrumSimilarity(top_spec, bottom_spec, t=t, b=b,
                                             top.label=top_spec_name, bottom.label <- bottom_spec_name,
                                             xlim=c(0,1.5*max.x)), file=output_file_name)
  } else {
    sim <- SpectrumSimilarity(top_spec, bottom_spec, t=t, b=b,top.label=top_spec_name,
                              bottom.label <- bottom_spec_name,xlim=c(0,1.5*max.x))

  }
  mtext("Spectrum similarity", line=2)
  sim <- format(sim, digits=4)
  mtext(as.character(sim), line=1)

  return(sim)
}


