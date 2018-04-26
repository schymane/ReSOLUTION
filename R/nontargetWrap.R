#nontarget wrapper functions
# E. Schymanski, 5/1/2017
# Documentation 3/4/2017

# data(isotopes, package="enviPat")
# data(adducts, package="enviPat")


#' Wrapper Function for Package nontarget
#'
#' @description This wrapper function runs the nontarget workflow with
#' selected default settings for positive and negative mode, creating
#' many output files to interpret the results.
#'
#' @usage nontargetWrap(peaklist_file, desc, path, pos=TRUE, inclHomol=FALSE,
#' colnums=c(1,3,5))
#'
#' @param peaklist_file File containing the peak list
#' @param desc Description - used to create file names
#' @param path Path in which to save the results files
#' @param pos Define the mode. Default \code{TRUE} indicates positive data,
#' \code{FALSE} indicates negative data.
#' @param inclHomol Indicate whether to include homologue series detection. Default \code{FALSE}.
#' If \code{TRUE}, all plus certain selected series are screened. Can be very
#' time consuming and lead to memory errors.
#' @param colnums Define the location of the mz, int and RT columns in \code{peaklist_file}.
#' The default \code{c(1,3,5)} works on enviPick output. Must be redefined for other peaklists
#' to avoid strange results!
#'
#' @return Results are saved into many output files in \code{path}
#'
#' @author Emma Schymanski (<emma.schymanski@@uni.lu>, wrapper), Martin Loos
#' (\code{\link{nontarget}}).
#'
#' @seealso \code{\link{nontarget}}.
#'
#' @export
#'
#' @examples
#'
#' #Note: this produces a lot of files, caution before trying.
#' peaklist_path <- system.file("extdata","Blind_nano_IS_toPeak750.txt",package="ReSOLUTION")
#' test_path <- paste0(getwd(),"/nontargettest")
#' dir.create(test_path)
#' nontargetWrap(peaklist_path, "BlindNanoIS_testWrap_to750", test_path, pos=TRUE, inclHomol=FALSE, colnums=c(1,3,5))
#' gc() #clear memory after.
#'
#'
nontargetWrap <- function(peaklist_file, desc, path, pos=TRUE, inclHomol=FALSE,
                          colnums=c(1,3,5)) {
  data(isotopes, package="enviPat")
  data(adducts, package="enviPat")

  peaklist <- read.table(peaklist_file)
  peaklist <- peaklist[,colnums]
  # colnums 1,3,5 selects the correct columns from enviPick output.
  # for other files, have to redefine. mz, int and RT
  dir_name <- desc
  full_dir_name <- paste(path,desc,sep="/")
  dir.create(full_dir_name)
  setwd(full_dir_name)

  #remove satelite peaks
  peaklist <- rm.sat(peaklist, dmz=0.3, drt=(0.1*60), intrat=0.015, spar=0.8,
                     corcut=-1000,plotit=TRUE)
  peaklist<-peaklist[peaklist[,4],1:3]

  # isotope pattern grouping
  # define the isotopes
  if (pos) {
      iso<-make.isos(isotopes,
                 #take out 41K as this seems to only cause problems  later
                 use_isotopes=c("13C","15N","34S","37Cl","81Br","13C","15N","34S","37Cl","81Br"),
                                  use_charges=c(1,1,1,1,1,2,2,2,2,2))
                 #use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
                 #use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
  } else {
    iso<-make.isos(isotopes,
                   #take out 41K as this seems to only cause problems  later
                   use_isotopes=c("13C","15N","34S","37Cl","81Br","13C","15N","34S","37Cl","81Br"),
                   #use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
                   use_charges=c(-1,-1,-1,-1,-1,-2,-2,-2,-2,-2))
  }
  # run the grouping
  pattern<-pattern.search(peaklist, iso, cutint=1000, rttol=60*c(-0.05,0.05), mztol=2,
                          mzfrac=0.1, ppm=TRUE, inttol=0.2,
                          rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                          deter=FALSE, entry=50)
  # plot the results
  plot_name_ends <- c("_PlotIsotopes.pdf", "_PlotDefect_N.pdf", "_PlotDefect_C.pdf",
                      "_PlotDefect_S.pdf","_PlotDefect_Cl.pdf","_PlotDefect_Br.pdf","_PlotDefect_K.pdf")
  plot_names <- paste(dir_name,plot_name_ends,sep="")
  plot_titles <- sub(".pdf","", plot_names)
  iso_plot_name <- paste(dir_name,"_PlotIsotopesAndElements.pdf",sep="")
  pdf(iso_plot_name,paper="a4r")
  a <- plotisotopes(pattern)
  mtext(plot_titles[1], line=1)
  plotdefect(pattern,elements=c("N"))
  mtext(plot_titles[2], line=1)
  plotdefect(pattern,elements=c("C"))
  mtext(plot_titles[3], line=1)
  plotdefect(pattern,elements=c("S"))
  mtext(plot_titles[4], line=1)
  plotdefect(pattern,elements=c("Cl"))
  mtext(plot_titles[5], line=1)
  plotdefect(pattern,elements=c("Br"))
  mtext(plot_titles[6], line=1)
  #plotdefect(pattern,elements=c("K"))
  #mtext(plot_titles[7], line=1)
  dev.off()

  # text file output
  write.table(pattern[[1]],file=paste(dir_name,"_pattern.txt", sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(a,file=paste(dir_name,"_IsotopeSummary.txt",sep=""),
              row.names=FALSE,quote=FALSE)

  # run grouping for the adducts and plot
  if (pos) {
    adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
                          use_adducts=c("M+K","M+H","M+Na","M+NH4"), ion_mode="positive")
  } else {
    adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
                          use_adducts=c("M-H","M+FA-H"), ion_mode="negative")
    # adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
    #                       use_adducts=c("M-H","M+FA-H","M+Cl"), ion_mode="negative")
  }

  # show single pattern groups and relations
  pdf(paste(dir_name,"_PlotAdducts_IsoAdd.pdf",sep=""), paper="a4r")
  b <- plotadduct(adduct)
  mtext(paste(dir_name,"_PlotAdduct",sep=""), line=1)
  plotall(pattern,adduct)
  mtext(paste(dir_name,"_PlotAllIsoAdduct",sep=""), line=1)
  dev.off()
  # plot for one group only - but wait until components are grouped later
  #plotgroup(pattern,adduct,groupID=1,massrange=10,allmass=FALSE);
  # text file output
  write.table(adduct[[1]],file=paste(dir_name,"_adduct.txt", sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(b,file=paste(dir_name,"_AdductSummary.txt", sep=""),
              row.names=TRUE,quote=TRUE)#,sep=";")

  if (inclHomol) {
  # Screen for homologous series
  #RT is in seconds from enviPick, need much bigger tolerances than default (for min)
  homol<-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"),use_C=TRUE,
                      #charge=c(-1,-2),
                      #charge=c(1,2), use_C=TRUE,
                      minmz=5, maxmz=60, minrt=0, maxrt=120,
                      ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=FALSE, vec_size=3E6)
  #For a few distinct masses # C2H4: 28.0313 # C2H4O: 44.026215 # CH2: 14.0compI1565 # C3H6O: 58.041865
  homol_select <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"),use_C=TRUE,
                              #charge=c(-1),
                              #charge=c(1), use_C=TRUE,
                              minmz=5, maxmz=60, minrt=0, maxrt=120,
                              ppm=TRUE, mztol=3.5, rttol=30, minlength=4, vec_size=3E6,
                              mzfilter=c(14.01565,28.0313,44.026215,58.041865))
  # for mz14 only:
  homol_mz14 <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"), use_C=TRUE,
                            #charge=c(-1),
                            # charge=c(1), use_C=TRUE,
                            minmz=14, maxmz=15, minrt=0, maxrt=120,
                            ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=14.01565, vec_size=3E6)
  # for mz44 only:
  homol_mz44 <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"), use_C=TRUE,
                            # charge=c(1), use_C=TRUE,
                            #charge=c(-1),
                            minmz=44.026215, maxmz=45, minrt=0, maxrt=120,
                            ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=44.026215, vec_size=3E6)

  # do the plots and writing output
  pdf(paste(dir_name,"_HomologueSeries.pdf",sep=""), paper="a4r")
  plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
  mtext(paste(dir_name,"_HomologueSeries_mz5to60_rt0to120",sep=""),line=1)
  # mz 14,28,44,58
  plothomol(homol_select,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
  mtext(paste(dir_name,"_HomologueSeries_mz14_28_44_58_rt0to120",sep=""),line=1)
  # mz 14
  plothomol(homol_mz14,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
  mtext(paste(dir_name,"_HomologueSeries_mz14_rt0to120",sep=""),line=1)
  # mz 44
  plothomol(homol_mz44,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
  mtext(paste(dir_name,"_HomologueSeries_mz44_rt0to120",sep=""),line=1)
  dev.off()

  # all
  write.table(homol[[1]],file=paste(dir_name,"_mz5to60_homol1.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(homol[[3]],file=paste(dir_name,"_mz5to60_homol3.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  # mz 14,28,44,58
  write.table(homol_select[[1]],file=paste(dir_name,"_mz14_28_44_58_homol1.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(homol_select[[3]],file=paste(dir_name,"_mz14_28_44_58_homol3.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  # mz 14
  write.table(homol_mz14[[1]],file=paste(dir_name,"_mz14_homol1.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(homol_mz14[[3]],file=paste(dir_name,"_mz14_homol3.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  # mz 44
  write.table(homol_mz44[[1]],file=paste(dir_name,"_mz44_homol1.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  write.table(homol_mz44[[3]],file=paste(dir_name,"_mz44_homol3.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")



  #combining grouping results into components
  comp<-combine(pattern, adduct, homol,
                dont=FALSE, rules=c(TRUE,FALSE,FALSE))
  } else {
    comp<-combine(pattern, adduct, homol=FALSE,
                  dont=FALSE, rules=c(TRUE,FALSE,FALSE))
  }

  pdf(paste(dir_name,"_PlotIsotopes_Comp.pdf",sep=""), paper="a4r")
  isoPerComp <- plotisotopes(comp)
  mtext(paste(dir_name,"_PlotIsotopes_Components",sep=""),line=1)
  dev.off()

  write.table(comp[[1]],file=paste(dir_name,"_comp1.txt",sep=""),
              row.names=FALSE,quote=TRUE)#,sep=";")
  comp1 <- read.table(paste(dir_name,"_comp1.txt",sep=""), header=TRUE)
  write.table(isoPerComp,file=paste(dir_name,"_IsoPerComponent.txt",sep=""),
              row.names=FALSE,quote=TRUE)

  # create subdir for component peaklists
  subdir <- paste(full_dir_name,"/MS1_Comp_all/",sep="")
  dir.create(subdir)
  compID <- 0
  subdir2 <- paste(full_dir_name,"/MS1_Comp_trim/",sep="")
  dir.create(subdir2)

  pdf(paste(dir_name,"_Plot_Individual_Components.pdf",sep=""), paper="a4r")
  for(i in 1:length(comp[[1]][,1])) {
    # check whether isotope or adduct group exists
    ID_PatternGroup <- as.character(comp[[1]][i,2]) #stops this being a factor but allows string comparison of "-"
    ID_AdductGroup <- as.character(comp[[1]][i,4])
    compID <- comp[[1]][i,1]
    ID_test <- length(grep("-",ID_PatternGroup)) + length(grep("-",ID_AdductGroup))
    # 0 if neither are "-", 1 if one is "-", 2 if both are "-" and then we don't want peak lists
    if(ID_test<2) {
      a <- plotcomp(comp,compoID=compID,peakID=FALSE)
      mtext(paste(dir_name,"_PlotComponent_",compID,sep=""),line=1)
      #do files containing all peaks in range
      filename <- paste(subdir,dir_name,"_CompPeakListAll_",compID,".txt",sep="")
      #comp_peaklist <- cbind(a$'concerned peaks'[2],a$'concerned peaks'[3])
      comp_peaklist <- cbind(a$'all peaks within range'[2],a$'all peaks within range'[3])
      write.table(comp_peaklist,filename,row.names=FALSE,col.names=FALSE)
      # do files with only the selected peaks
      filename <- paste(subdir2,dir_name,"_CompPeakList_",compID,".txt",sep="")
      #comp_peaklist <- cbind(a$'concerned peaks'[2],a$'concerned peaks'[3])
      comp_peaklist <- cbind(a$a[2],a$a[3])
      write.table(comp_peaklist,filename,row.names=FALSE,col.names=FALSE)

    }
  }
  dev.off()
  #return()
}


# nontargetWrap.csv.FIA <- function(peaklist_file, desc, path, pos=TRUE, inclHomol=FALSE,
#                               remSat=TRUE) {
#   peaklist <- read.csv(peaklist_file)
#   # This assumes peaklist as csv has mz, int and RT in correct order
#   #peaklist <- peaklist[,c(1,3,5)]
#   dir_name <- desc
#   full_dir_name <- paste(path,desc,sep="/")
#   dir.create(full_dir_name)
#   setwd(full_dir_name)
#
#   #remove satelite peaks
#   if (remSat) {
#     peaklist <- rm.sat(peaklist, dmz=0.3, drt=(0.1*60), intrat=0.015, spar=0.8,
#                        corcut=-1000,plotit=TRUE)
#     peaklist<-peaklist[peaklist[,4],1:3]
#   }
#
#   # isotope pattern grouping
#   # define the isotopes
#   if (pos) {
#     iso<-make.isos(isotopes,
#                    #take out 41K as this seems to only cause problems  later
#                    use_isotopes=c("13C","15N","34S","37Cl","81Br","13C","15N","34S","37Cl","81Br"),
#                    use_charges=c(1,1,1,1,1,2,2,2,2,2))
#     #use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
#     #use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2))
#     #use_charges=c(-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2)
#   } else {
#     iso<-make.isos(isotopes,
#                    #take out 41K as this seems to only cause problems  later
#                    use_isotopes=c("13C","15N","34S","37Cl","81Br","13C","15N","34S","37Cl","81Br"),
#                    #use_charges=c(1,1,1,1,1,2,2,2,2,2)
#                    #use_isotopes=c("13C","15N","34S","37Cl","81Br","41K","13C","15N","34S","37Cl","81Br","41K"),
#                    #use_charges=c(1,1,1,1,1,1,2,2,2,2,2,2)
#                    use_charges=c(-1,-1,-1,-1,-1,-2,-2,-2,-2,-2))
#
#   }
#   # run the grouping
#
# #   pattern<-pattern.search(peaklist, iso, cutint=1000, rttol=60*c(-0.05,0.05), mztol=2,
# #                           mzfrac=0.1, ppm=TRUE, inttol=0.2,
# #                           rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
# #                           deter=FALSE, entry=50)
#   # try some new params for FIA / TOF.
#   # First try: kept cutint at 1000, mztol down to 1 from 2, inttol to 0.01 from 0.2
#   pattern<-pattern.search(peaklist, iso, cutint=1000, rttol=60*c(-0.05,0.05), mztol=1,
#                           mzfrac=0.1, ppm=TRUE, inttol=0.01,
#                           rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
#                           deter=FALSE, entry=50)
#   # plot the results
#   plot_name_ends <- c("_PlotIsotopes.pdf", "_PlotDefect_N.pdf", "_PlotDefect_C.pdf",
#                       "_PlotDefect_S.pdf","_PlotDefect_Cl.pdf","_PlotDefect_Br.pdf","_PlotDefect_K.pdf")
#   plot_names <- paste(dir_name,plot_name_ends,sep="")
#   plot_titles <- sub(".pdf","", plot_names)
#   iso_plot_name <- paste(dir_name,"_PlotIsotopesAndElements.pdf",sep="")
#   pdf(iso_plot_name,paper="a4r")
#   a <- plotisotopes(pattern)
#   mtext(plot_titles[1], line=1)
#   plotdefect(pattern,elements=c("N"))
#   mtext(plot_titles[2], line=1)
#   plotdefect(pattern,elements=c("C"))
#   mtext(plot_titles[3], line=1)
#   plotdefect(pattern,elements=c("S"))
#   mtext(plot_titles[4], line=1)
#   plotdefect(pattern,elements=c("Cl"))
#   mtext(plot_titles[5], line=1)
#   plotdefect(pattern,elements=c("Br"))
#   mtext(plot_titles[6], line=1)
#   #plotdefect(pattern,elements=c("K"))
#   #mtext(plot_titles[7], line=1)
#   dev.off()
#
#   # text file output
#   write.table(pattern[[1]],file=paste(dir_name,"_pattern.txt", sep=""),
#               row.names=FALSE,quote=TRUE)#,sep=";")
#   write.table(a,file=paste(dir_name,"_IsotopeSummary.txt",sep=""),
#               row.names=FALSE,quote=FALSE)
#
#   # run grouping for the adducts and plot
#   if (pos) {
#     adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
#                           use_adducts=c("M+K","M+H","M+Na","M+NH4"), ion_mode="positive")
#   } else {
#     adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
#                           use_adducts=c("M-H","M+FA-H"), ion_mode="negative")
#     # adduct<-adduct.search(peaklist, adducts, rttol=60*0.05, mztol=3, ppm=TRUE,
#     #                       use_adducts=c("M-H","M+FA-H","M+Cl"), ion_mode="negative")
#   }
#
#   # show single pattern groups and relations
#   pdf(paste(dir_name,"_PlotAdducts_IsoAdd.pdf",sep=""), paper="a4r")
#   b <- plotadduct(adduct)
#   mtext(paste(dir_name,"_PlotAdduct",sep=""), line=1)
#   plotall(pattern,adduct)
#   mtext(paste(dir_name,"_PlotAllIsoAdduct",sep=""), line=1)
#   dev.off()
#   # plot for one group only - but wait until components are grouped later
#   #plotgroup(pattern,adduct,groupID=1,massrange=10,allmass=FALSE);
#   # text file output
#   write.table(adduct[[1]],file=paste(dir_name,"_adduct.txt", sep=""),
#               row.names=FALSE,quote=TRUE)#,sep=";")
#   write.table(b,file=paste(dir_name,"_AdductSummary.txt", sep=""),
#               row.names=TRUE,quote=TRUE)#,sep=";")
#
#   if (inclHomol) {
#     # Screen for homologous series
#     #RT is in seconds from enviPick, need much bigger tolerances than default (for min)
#     homol<-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"),use_C=TRUE,
#                         #charge=c(-1,-2),
#                         #charge=c(1,2), use_C=TRUE,
#                         minmz=5, maxmz=60, minrt=0, maxrt=120,
#                         ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=FALSE, vec_size=3E6)
#     #For a few distinct masses # C2H4: 28.0313 # C2H4O: 44.026215 # CH2: 14.0compI1565 # C3H6O: 58.041865
#     homol_select <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"),use_C=TRUE,
#                                 #charge=c(-1),
#                                 #charge=c(1), use_C=TRUE,
#                                 minmz=5, maxmz=60, minrt=0, maxrt=120,
#                                 ppm=TRUE, mztol=3.5, rttol=30, minlength=4, vec_size=3E6,
#                                 mzfilter=c(14.01565,28.0313,44.026215,58.041865))
#     # for mz14 only:
#     homol_mz14 <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"), use_C=TRUE,
#                               #charge=c(-1),
#                               # charge=c(1), use_C=TRUE,
#                               minmz=14, maxmz=15, minrt=0, maxrt=120,
#                               ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=14.01565, vec_size=3E6)
#     # for mz44 only:
#     homol_mz44 <-homol.search(peaklist, isotopes, elements=c("C","H","O","S","N"), use_C=TRUE,
#                               # charge=c(1), use_C=TRUE,
#                               #charge=c(-1),
#                               minmz=44.026215, maxmz=45, minrt=0, maxrt=120,
#                               ppm=TRUE, mztol=3.5, rttol=30, minlength=4, mzfilter=44.026215, vec_size=3E6)
#
#     # do the plots and writing output
#     pdf(paste(dir_name,"_HomologueSeries.pdf",sep=""), paper="a4r")
#     plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
#     mtext(paste(dir_name,"_HomologueSeries_mz5to60_rt0to120",sep=""),line=1)
#     # mz 14,28,44,58
#     plothomol(homol_select,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
#     mtext(paste(dir_name,"_HomologueSeries_mz14_28_44_58_rt0to120",sep=""),line=1)
#     # mz 14
#     plothomol(homol_mz14,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
#     mtext(paste(dir_name,"_HomologueSeries_mz14_rt0to120",sep=""),line=1)
#     # mz 44
#     plothomol(homol_mz44,xlim=FALSE,ylim=FALSE,plotlegend=TRUE)
#     mtext(paste(dir_name,"_HomologueSeries_mz44_rt0to120",sep=""),line=1)
#     dev.off()
#
#     # all
#     write.table(homol[[1]],file=paste(dir_name,"_mz5to60_homol1.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     write.table(homol[[3]],file=paste(dir_name,"_mz5to60_homol3.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     # mz 14,28,44,58
#     write.table(homol_select[[1]],file=paste(dir_name,"_mz14_28_44_58_homol1.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     write.table(homol_select[[3]],file=paste(dir_name,"_mz14_28_44_58_homol3.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     # mz 14
#     write.table(homol_mz14[[1]],file=paste(dir_name,"_mz14_homol1.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     write.table(homol_mz14[[3]],file=paste(dir_name,"_mz14_homol3.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     # mz 44
#     write.table(homol_mz44[[1]],file=paste(dir_name,"_mz44_homol1.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#     write.table(homol_mz44[[3]],file=paste(dir_name,"_mz44_homol3.txt",sep=""),
#                 row.names=FALSE,quote=TRUE)#,sep=";")
#
#
#
#     #combining grouping results into components
#     comp<-combine(pattern, adduct, homol,
#                   dont=FALSE, rules=c(TRUE,FALSE,FALSE))
#   } else {
#     comp<-combine(pattern, adduct, homol=FALSE,
#                   dont=FALSE, rules=c(TRUE,FALSE,FALSE))
#   }
#
#   pdf(paste(dir_name,"_PlotIsotopes_Comp.pdf",sep=""), paper="a4r")
#   isoPerComp <- plotisotopes(comp)
#   mtext(paste(dir_name,"_PlotIsotopes_Components",sep=""),line=1)
#   dev.off()
#
#   write.table(comp[[1]],file=paste(dir_name,"_comp1.txt",sep=""),
#               row.names=FALSE,quote=TRUE)#,sep=";")
#   comp1 <- read.table(paste(dir_name,"_comp1.txt",sep=""), header=TRUE)
#   write.table(isoPerComp,file=paste(dir_name,"_IsoPerComponent.txt",sep=""),
#               row.names=FALSE,quote=TRUE)
#
#   # create subdir for component peaklists
#   subdir <- paste(full_dir_name,"/MS1_Comp_all/",sep="")
#   dir.create(subdir)
#   compID <- 0
#   subdir2 <- paste(full_dir_name,"/MS1_Comp_trim/",sep="")
#   dir.create(subdir2)
#
#   pdf(paste(dir_name,"_Plot_Individual_Components.pdf",sep=""), paper="a4r")
#   for(i in 1:length(comp[[1]][,1])) {
#     # check whether isotope or adduct group exists
#     ID_PatternGroup <- as.character(comp[[1]][i,2]) #stops this being a factor but allows string comparison of "-"
#     ID_AdductGroup <- as.character(comp[[1]][i,4])
#     compID <- comp[[1]][i,1]
#     ID_test <- length(grep("-",ID_PatternGroup)) + length(grep("-",ID_AdductGroup))
#     # 0 if neither are "-", 1 if one is "-", 2 if both are "-" and then we don't want peak lists
#     if(ID_test<2) {
#       a <- plotcomp(comp,compoID=compID,peakID=FALSE)
#       mtext(paste(dir_name,"_PlotComponent_",compID,sep=""),line=1)
#       #do files containing all peaks in range
#       filename <- paste(subdir,dir_name,"_CompPeakListAll_",compID,".txt",sep="")
#       #comp_peaklist <- cbind(a$'concerned peaks'[2],a$'concerned peaks'[3])
#       comp_peaklist <- cbind(a$'all peaks within range'[2],a$'all peaks within range'[3])
#       write.table(comp_peaklist,filename,row.names=FALSE,col.names=FALSE)
#       # do files with only the selected peaks
#       filename <- paste(subdir2,dir_name,"_CompPeakList_",compID,".txt",sep="")
#       #comp_peaklist <- cbind(a$'concerned peaks'[2],a$'concerned peaks'[3])
#       comp_peaklist <- cbind(a$a[2],a$a[3])
#       write.table(comp_peaklist,filename,row.names=FALSE,col.names=FALSE)
#
#     }
#   }
#   dev.off()
#   #return()
# }


