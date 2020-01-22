# setting up configuration files for MetFragBetaPlus and performing the calcs
# E. Schymanski, 9/6/2015
# Rev. 1 16/11/2015 - MetFrag2.2
# NOTE: 15/2/2017: Moving to ReSOLUTION package
# C:/DATA/R/ReSOLUTION_scripts/MetFragConfigR.R
# NOTE 4/4/2018: MetFrag2.4.3 & 2.4.4MS-ready documenting etc.
# fresh download from http://c-ruttkies.github.io/MetFrag/projects/metfragcl/

# requires readxl

# set up default MetFrag Config Files
# minimum info is mass, adduct type, msms peak list name

#' Create MetFrag Configuration Files
#'
#' @description This function provides options to set up configuration files to run MetFrag Command Line
#' in batch mode. Minimum information is mass, adduct type and MS/MS peak list. MetFrag Command Line
#' is available from \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/}
#'
#' @usage MetFragConfig(mass, adduct_type, results_filename, peaklist_path, base_dir,
#' DB=c("PubChem"),
#' localDB_path="", output="XLS", token="", neutralPrecursorMass=FALSE, 
#' ppm=5, mzabs=0.001, frag_ppm=5, IsPosMode=TRUE,
#' tree_depth=2, num_threads=1, add_refs=TRUE, minInt=0, rt_file_path="", rt_exp=0,suspect_path="",
#' suspect_filter=FALSE, UDS_Category="", UDS_Weights="", DB_IDs="", mol_form="", useFormula=FALSE,
#' useMoNAMetFusion=TRUE, useMonaIndiv=TRUE, MoNAoffline=TRUE, incl_el="",excl_el="", incl_exclusive=FALSE,
#' incl_smarts_filter="", incl_smarts_score="", excl_smarts_filter="",excl_smarts_score="", filter_isotopes=TRUE,
#' filter_by_InChIKey=TRUE)
#'
#' @param mass The mass with which to search the candidate database (\code{DB}). Use \code{neutralPrecursorMass} and
#' \code{adduct_type} to set whether this is monoisotopic mass or an adduct species.
#' @param adduct_type The adduct species used to define mass (if \code{neutralPrecursorMass=FALSE}) and fragmentation settings
#' in the config file, entered as either \code{PrecursorIonType} (text) or \code{PrecursorIonmode} (a number). The available
#' options are given in the system file \code{MetFragAdductTypes.csv} in the \code{extdata} folder. If 
#' \code{neutralPrecursorMass=TRUE}, set \code{adduct_type=0}.
#' Recommended default values (if ion state is unclear) are \code{[M+H]+} (1) for positive and \code{[M-H]-} (-1) for negative mode.
#' @param results_filename Enter a base filename for naming the results files - do not include file endings
#' @param peaklist_path Enter the full path and file name to the peak list for this config file
#' @param base_dir Enter the directory name to set up the subfolders for MetFrag batch results. If the folders don't exist,
#' subfolders \code{config}, \code{log} and \code{results} are created; the output of this function is saved in \code{config}.
#' @param DB Enter query database name. Current options \code{KEGG}, \code{PubChem}, \code{ExtendedPubChem}, \code{ChemSpider},
#' \code{FOR_IDENT}, \code{MetaCyc}, \code{LocalCSV}, \code{LocalPSV} or \code{LocalSDF}. For \code{HMDB}, \code{LipidMaps} and
#' \code{KEGG-derivatised} use the \code{LocalCSV} option with respective files downloaded from
#' \url{https://msbi.ipb-halle.de/~cruttkie/databases/}.
#' @param localDB_path Full path and file name to the local database for \code{LocalCSV, LocalPSV or localSDF}. Otherwise leave empty.
#' If the file is not found, the config file defaults to \code{DB=PubChem}.
#' @param output Select output format(s) desired. Current options include one or more of
#' \code{SDF, XLS, CSV, ExtendedXLS, ExtendedFragmentsXLS} entered as a string. Not tested; incorrect entries will lead to CL failure.
#' @param token ChemSpider token, only required for \code{DB=ChemSpider}. See \url{http://www.chemspider.com/MassSpecAPI.asmx} for
#' more details about which services require tokens and \url{http://www.chemspider.com/help-create-a-chemspider-account.aspx} for
#' information how to obtain your token. If an invalid token is provided (not length=36), \code{DB} defaults to \code{PubChem}.
#' @param neutralPrecursorMass Controls whether \code{mass} is treated as a neutral or charged mass. If \code{TRUE}, treated
#' as neutral. If \code{FALSE} (default), this is entered as a charged mass, adjusted in MetFragCL with the \code{adduct_type} setting.
#' @param mol_form A string containing the molecular formula (used in candidate retreival)
#' @param useFormula Default \code{FALSE} means an exact mass search is performed. If \code{TRUE}, \code{mol_form} must be given and
#' candidate retreival is based on this formula. Note some databases are sensitive to the order of elements in the formula.
#' @param DB_IDs Use this to select only certain candidates using (comma-separated) database identifiers consistent with \code{DB}.
#' @param ppm The ppm error to perform the exact mass search for candidate retrieval (default 5 ppm)
#' @param mzabs The absolute error (in Da/Th) used to match fragments to observed MS/MS peaks. Additive with \code{frag_ppm}.
#' Default 0.001 Da (Th).
#' @param frag_ppm The relative error (in ppm) used to match fragments to observed MS/MS peaks. Additive with \code{mzabs}.
#' Default 5 ppm.
#' @param IsPosMode Controls the mode for both candidate retrieval and fragmentation consistently. Default \code{TRUE} sets
#' positive mode, switch to \code{FALSE} for negative mode data.
#' @param tree_depth Sets the number of fragmentation steps. Default=2 is recommended. Higher values lead to long calculation times.
#' @param num_threads Sets the number of threads used to run calculations. Default=1; set higher for faster results.
#' @param add_refs If set to (default) \code{TRUE}, reference scoring terms will be added for \code{DB=PubChem} and \code{DB=ChemSpider}.
#' Two terms (references, patents) are added for \code{PubChem}, weighted 0.5; four terms weighted 0.25 for \code{ChemSpider}.
#' These setttings can be overwritten by setting \code{add_refs=FALSE} and adding the desired terms to \code{UDS_Category} and
#' \code{UDS_weights}.
#' @param minInt Minimum intensity value to consider peaks in the MS/MS file. Default 0, this is merely a convenience function to
#' allow users to do a bare minimum noise reduction if required.
#' @param rt_file_path Full path to the CSV file containing InChIs and retention times (RTs) of standards to build the RT model.
#' The file should contain two column separated columns with a header row with the column names \code{InChI} and \code{RetentionTime}.
#' The example system file \code{Eawag_rt_inchi.csv} in the \code{extdata} folder is the correct dataset for Eawag MassBank
#' records measured on the XBridge C18 column.
#' @param rt_exp The experimental retention time. The chromatography and RT unit must match with the file in \code{rt_file_path}.
#' @param suspect_path Path to the suspect lists to be used as a filter or scoring term.
#' @param suspect_filter Default \code{FALSE} means suspect lists in \code{suspect_path} are used to increase the score of
#' candidates present in the suspect lists given (added as a scoring term).
#' If \code{TRUE}, suspect lists are used as a filter instead (only candidates present in the suspect lists are processed).
#' @param UDS_Category A string containing the exact column headers of additional User Defined Scores (UDS) to use, separated by
#' a comma. These column headers must match exactly, cannot be repeated and
#' must be present in the default database chosen or in the LocalCSV, PSV or SDF files used as a local database.
#' This can also be used to overwrite the default reference information in \code{add_refs}.
#' @param UDS_Weights A string containing comma-separated weight values corresponding to \code{UDS_Category}. This must
#' match exactly or an exception is thrown during processing.
#' @param useMoNAMetFusion Default \code{TRUE} means that the MoNA MetFusion Score is added by default. Use \code{FALSE} to exclude.
#' @param useMonaIndiv Default \code{TRUE} means that the MoNA Individual Score is added by default. Use \code{FALSE} to exclude.
#' This performs a direct lookup by InChIKey and returns the highest similarity value over all matches. A good match is a very good sign;
#' a poor match means there is a spectrum in MoNA for that compound but this may have been recorded with vastly different settings, so
#' a poor match does not necessarily indicate that the candidate is wrong.
#' @param MoNAoffline Default \code{TRUE} means the local MoNA instance (in the jar file) is used to avoid server issues. Use
#' \code{FALSE} to perform this live, however this may not work.
#' @param incl_el A string containing comma-separated elements that must be present in candidates. This allows coupling of an
#' exact mass search with the presence of elements containing distinct isotope patterns.
#' @param excl_el A string containing comma-separated elements that must not be present in candidates. This allows coupling of an
#' exact mass search with the absence of elements containing distinct isotope patterns.
#' @param incl_exclusive Default \code{FALSE} indicates that the elements in \code{incl_el} must be present, but other elements
#' could still be present. If \code{TRUE}, only these elements are considered (use this option with caution!)
#' @param incl_smarts_filter A string containing SMARTS codes (comma-separated) used to define substructures present (candidates that
#' do not contain these SMARTS are filtered out).
#' @param incl_smarts_score A string containing SMARTS codes (comma-separated) used to increase the score of candidates with
#' certain substructures present.
#' @param excl_smarts_filter A string containing SMARTS codes to exclude candidates with these substructures present.
#' @param excl_smarts_score A string containing SMARTS codes to penalize candidate scores with these substructures present.
#' @param filter_isotopes Default \code{TRUE} removes all candidates containing non-standard isotopes.
#' @param filter_by_InChIKey Default \code{TRUE} collapses the candidate result lists by the first block of the InChIKey, presenting
#' only the candidate with the best score across all categories. If \code{FALSE}, all candiates are included in the results.
#'
#' @return Creates a MetFrag config file matching the given parameters and returns the file name.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL author).
#'
#' @seealso \code{\link{runMetFrag}} to run the config files.
#'
#' @export
#'
#' @examples
#' # Do not run unless you adjusted test_dir to an existing file location
#' peaklist_path <- system.file("extdata","EA026206_Simazine_peaks.txt",package="ReSOLUTION")
#' # change this directory to an existing one, or this example won't work
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' testCSV <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
#'
#' config_file <- MetFragConfig(201.0776,"[M+H]+","Simazine_neutralMass_PubChem",peaklist_path, test_dir, DB="PubChem",neutralPrecursorMass=TRUE)
#' config_file2 <- MetFragConfig(202.0854,1,"Simazine_precMass_localCSV",peaklist_path,test_dir,DB="LocalCSV",localDB_path=testCSV)
#' config_file2 <- MetFragConfig(202.0854,1,"Simazine_precMass_10ppm",peaklist_path,test_dir,DB="LocalCSV",localDB_path=testCSV,ppm=10)
#' config_file2 <- MetFragConfig(202.0854,1,"Simazine_precMass_10ppm_InChIFilterOff",peaklist_path,test_dir,DB="LocalCSV",
#'                               localDB_path=testCSV,ppm=10,filter_by_InChIKey = FALSE)
#'
#' #to find out the adduct states:
#' MetFragAdductTypes <- read.csv(system.file("extdata","MetFrag_AdductTypes.csv",package="ReSOLUTION"))
#'
#' # to run the config files
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' # warning: this first query takes a while, for quick testing run config_file2
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#' runMetFrag(config_file2, metfrag_dir, MetFragCL_name)
#'
MetFragConfig <- function(mass, adduct_type, results_filename, peaklist_path, base_dir,
                          DB=c("PubChem"),
                          localDB_path="",output="XLS", token="",
                          neutralPrecursorMass=FALSE, mol_form="",useFormula=FALSE,DB_IDs="",
                          ppm=5, mzabs=0.001, frag_ppm=5, IsPosMode=TRUE,
                          tree_depth=2, num_threads=1, add_refs=TRUE, minInt=0,
                          rt_file_path="",rt_exp=0,suspect_path="",suspect_filter=FALSE,
                          UDS_Category="",UDS_Weights="",
                          useMoNAMetFusion=TRUE,useMonaIndiv=TRUE,MoNAoffline=TRUE,
                          incl_el="",excl_el="",incl_exclusive=FALSE,
                          incl_smarts_filter="",incl_smarts_score="",
                          excl_smarts_filter="",excl_smarts_score="",
                          filter_isotopes=TRUE,filter_by_InChIKey=TRUE) {
  # check database definition
  # @param DB Enter query database name. Current options \code{KEGG}, \code{PubChem}, \code{ExtendedPubChem}, \code{ChemSpider},
  # \code{FOR_IDENT}, \code{MetaCyc}, \code{LocalCSV}, \code{LocalPSV} or \code{LocalSDF}. For \code{HMDB}, \code{LipidMaps} and
  # \code{KEGG-derivatised} use the \code{LocalCSV} option with respective files downloaded from
  # \url{https://msbi.ipb-halle.de/~cruttkie/databases/}. Coming soon: CompTox Dashboard.

  if(length(DB) > 1) {
    stop("Please define only one database")
  }
  if(!(DB %in% c("KEGG","PubChem","ExtendedPubChem","ChemSpider", "FOR-IDENT", "MetaCyc",
                 "LocalCSV", "LocalPSV", "LocalSDF"))) {
    stop("Incorrect database: select one of KEGG, PubChem, ExtendedPubChem, ChemSpider, FOR-IDENT, MetaCyc,
         LocalCSV, LocalPSV or LocalSDF")
  }
  # check adduct type definition
  MetFragAdductTypes <- read.csv(system.file("extdata","MetFrag_AdductTypes.csv",package="ReSOLUTION"))
  adduct_num_list <- MetFragAdductTypes$PrecursorIonMode
  adduct_name_list <- MetFragAdductTypes$PrecursorIonType
  isPosIonMode_list <- MetFragAdductTypes$IsPositiveIonMode
  adduct_type_name <- grep(adduct_type,MetFragAdductTypes$PrecursorIonType,fixed=TRUE)
  if (length(adduct_type_name)<1) {
    adduct_type_name <- FALSE
  }
  adduct_type_num <- adduct_type %in%  MetFragAdductTypes$PrecursorIonMode
  if (is.na(adduct_type_num)) {
    adduct_type_num <- FALSE
  }
  if(!(adduct_type_name || adduct_type_num)) {
    stop(paste0("Incorrect adduct type: supported adducts are listed in ",
                system.file("extdata","MetFrag_AdductTypes.csv",package="ReSOLUTION")))
  }
  # adduct settings are compulsory, as they define fragmentation and (potentially) mass settings

  # set up files and paths
  peaklist <- peaklist_path
  config_dir <- paste(base_dir,"/config",sep="")
  if (!file.exists(config_dir)) {
    dir.create(config_dir)
  }
  # add test for whether dir exists
  config_name <- paste(config_dir,"/",paste(results_filename,"_config.txt", sep=""),sep="")
  SampleName <- results_filename
  ResultsPath <- paste(base_dir,"/results",sep="")
  LocalDatabasePath <- localDB_path
  if (!file.exists(ResultsPath)) {
    dir.create(ResultsPath)
  }

  # Open file and write parameters after defining them
  file.create(config_name)
  file.conn <- file(config_name)
  open(file.conn,open="at")
  writeLines(paste("SampleName = ",as.character(SampleName),sep=""),con=file.conn)
  writeLines(paste("ResultsPath = ",as.character(ResultsPath),sep=""),con=file.conn)
  writeLines(paste("PeakListPath = ",as.character(peaklist),sep=""),con=file.conn)

  ## Setting up databases and parameters
  # adduct setting compulsory:
  # if(!(adduct_type %in% (adduct_num_list || adduct_name_list))) {
  #   stop(paste0("Incorrect adduct type: supported adducts are listed in ",
  #               system.file("extdata","MetFrag_AdductTypes.csv",package="ReSOLUTION")))
  # }

  if (adduct_type_name) {
    # retrieve the adduct mode to test this
    adduct_mode <- isPosIonMode_list[grep(adduct_type,adduct_name_list,fixed=TRUE)]
    # if we have a name match and if the modes match, print PrecursorIonType, else stop.
    if (adduct_mode == IsPosMode) {
      writeLines(paste("PrecursorIonType = ",as.character(adduct_type),sep=""),con=file.conn)
    } else {
      stop("The adduct_type and IsPosMode settings mismatch. Please check and retry")
    }
  } else if (adduct_type_num) {
    # retrieve adduct mode to test
    adduct_mode <- isPosIonMode_list[match(adduct_type,adduct_num_list)]
    if (adduct_type %in%  MetFragAdductTypes$PrecursorIonMode[1]) {
      # if adduct type is zero, don't check, only write out if not neutralPrecursorMass
      if (!neutralPrecursorMass) {
        writeLines(paste("PrecursorIonMode = ",as.character(adduct_type),sep=""),con=file.conn)
      }
    } else if (adduct_mode == IsPosMode) {
      # if we have a number match and the modes match, print PrecursorIonMode, else stop
      writeLines(paste("PrecursorIonMode = ",as.character(adduct_type),sep=""),con=file.conn)
    } else {
      stop("The adduct_type and IsPosMode settings mismatch. Please check and retry")
    }
  } else {
    # this case should never happen, but just in case ...
    stop(paste0("Error in adduct settings, please check and try again. Supported adducts are listed in ",
                system.file("extdata","MetFrag_AdductTypes.csv",package="ReSOLUTION")))
  }

  # now set up exact masses and corrections
  #IonizedPrecursorMass <- mass

  if (useFormula) {
    writeLines(paste("NeutralPrecursorMolecularFormula = ",as.character(mol_form),sep=""),con=file.conn)
    if (neutralPrecursorMass) {
      writeLines(paste("NeutralPrecursorMass = ",as.character(mass),sep=""),con=file.conn)
    } else {
      writeLines(paste("IonizedPrecursorMass = ",as.character(mass),sep=""),con=file.conn)
    }
  } else if (neutralPrecursorMass) {
    writeLines(paste("NeutralPrecursorMass = ",as.character(mass),sep=""),con=file.conn)
  } else {
    writeLines(paste("IonizedPrecursorMass = ",as.character(mass),sep=""),con=file.conn)
  }
  if (IsPosMode) {
    IsPositiveIonMode <- "True"
  } else {
    IsPositiveIonMode <- "False"
  }
  # add test for adduct type here - checking the mode too
  DatabaseSearchRelativeMassDeviation <- ppm
  FragmentPeakMatchAbsoluteMassDeviation <- mzabs
  FragmentPeakMatchRelativeMassDeviation <- frag_ppm
  MinimumAbsolutePeakIntensity <- minInt # optional: to filter noise peaks
  MaximumTreeDepth <- tree_depth
  NumberThreads <- num_threads
  # SDF, XLS, CSV, ExtendedXLS, ExtendedFragmentsXLS
  MetFragCandidateWriter <- output

#
  writeLines(paste("IsPositiveIonMode = ",as.character(IsPositiveIonMode),sep=""),con=file.conn)
  writeLines(paste("DatabaseSearchRelativeMassDeviation = ",
                   as.character(DatabaseSearchRelativeMassDeviation),sep=""),con=file.conn)
  if (MinimumAbsolutePeakIntensity > 0) {
    writeLines(paste("MinimumAbsolutePeakIntensity = ",
                     as.character(MinimumAbsolutePeakIntensity),sep=""),con=file.conn)
  }
  writeLines(paste("FragmentPeakMatchAbsoluteMassDeviation = ",
                   as.character(FragmentPeakMatchAbsoluteMassDeviation),sep=""),con=file.conn)
  writeLines(paste("FragmentPeakMatchRelativeMassDeviation = ",
                   as.character(FragmentPeakMatchRelativeMassDeviation),sep=""),con=file.conn)
  writeLines(paste("MaximumTreeDepth = ",as.character(MaximumTreeDepth),sep=""),con=file.conn)
  # preprocessing, postprocessing and score information written later
  writeLines(paste("NumberThreads = ",as.character(NumberThreads),sep=""),con=file.conn)
  writeLines(paste("MetFragCandidateWriter = ",as.character(MetFragCandidateWriter),sep=""),con=file.conn)

  # Processing Options
  MetFragPreProcessingCandidateFilter <- "UnconnectedCompoundFilter"
  if (filter_isotopes) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",IsotopeFilter",sep="")
  } # this filters compounds with non-standard isotopes
  MetFragDatabaseType <- DB
  ChemSpiderToken <- token
  if (grepl("ChemSpider",DB) && nchar(ChemSpiderToken) < 36) {
    warning("Invalid ChemSpider Token, switching to PubChem")
    DB <- "PubChem"
  }
  if (add_refs && grepl("ChemSpider",DB)) {
    MetFragScoreTypes <- "FragmenterScore,ChemSpiderReferenceCount,ChemSpiderDataSourceCount,ChemSpiderNumberPubMedReferences,ChemSpiderRSCCount"
    MetFragScoreWeights <- "1.0,0.25,0.25,0.25,0.25"
  } else if (add_refs && grepl("PubChem",DB)) {
    MetFragDatabaseType <- "ExtendedPubChem"
    MetFragScoreTypes <- "FragmenterScore,PubChemNumberPubMedReferences,PubChemNumberPatents"
    MetFragScoreWeights <- "1.0,0.5,0.5"
  } else if (grepl("ChemSpider",DB) || grepl("PubChem",DB) || grepl("KEGG",DB)) {
    MetFragScoreTypes <- "FragmenterScore"
    MetFragScoreWeights <- "1.0"
  } else if (grepl("LocalCSV",DB) || grepl("LocalSDF",DB)) {
    MetFragScoreTypes <- "FragmenterScore"
    MetFragScoreWeights <- "1.0"
    if (file.exists(LocalDatabasePath)) {
      writeLines(paste("LocalDatabasePath = ",LocalDatabasePath,sep=""),con=file.conn)
    } else {
      warning("Local database file does not exist, defaulting to PubChem without references")
      MetFragDatabaseType <- "PubChem"
    }
  } else {
    warning("Database type incorrectly defined, defaulting to PubChem without references")
    MetFragDatabaseType <- "PubChem"
    MetFragScoreTypes <- "FragmenterScore"
    MetFragScoreWeights <- "1.0"
  }
  writeLines(paste("MetFragDatabaseType = ",MetFragDatabaseType,sep=""),con=file.conn)
  if ((nchar(ChemSpiderToken)>0)&& grepl("ChemSpider",MetFragDatabaseType,fixed=TRUE)) {
    writeLines(paste("ChemSpiderToken = ",ChemSpiderToken,sep=""),con=file.conn)
  }
 if (nchar(DB_IDs)>0) {
   writeLines(paste("PrecursorCompoundIDs = ",DB_IDs,sep=""),con=file.conn)
 }
 # if (nchar(mol_form)>0) {
 #   writeLines(paste("NeutralPrecursorMolecularFormula = ",mol_form,sep=""),con=file.conn)
 # }

  # check if retention time and training file is given, if so, add
  if (nchar(rt_file_path)>1 && rt_exp > 0) {
    if (file.exists(rt_file_path)) {
      MetFragScoreTypes <- paste(MetFragScoreTypes,",RetentionTimeScore",sep="")
      MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
      RetentionTimeTrainingFile <- rt_file_path
      ExperimentalRetentionTimeValue <- rt_exp
      writeLines(paste("RetentionTimeTrainingFile = ",
                       as.character(RetentionTimeTrainingFile),sep=""),con=file.conn)
      writeLines(paste("ExperimentalRetentionTimeValue = ",
                       as.character(ExperimentalRetentionTimeValue),sep=""),con=file.conn)
    } else {
      RetentionTimeTrainingFile <- NA
      ExperimentalRetentionTimeValue <- NA
    }
  }
  # check if suspect file is given, if so, add
  if (nchar(suspect_path)>1 && !suspect_filter) {
    MetFragScoreTypes <- paste(MetFragScoreTypes,",SuspectListScore",sep="")
    MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
    ScoreSuspectLists <- suspect_path
    FilterSuspectLists <- NA
    writeLines(paste("ScoreSuspectLists = ",as.character(ScoreSuspectLists),sep=""),con=file.conn)
  } else if (nchar(suspect_path)>1 && suspect_filter) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",SuspectListFilter",sep="")
    FilterSuspectLists <- suspect_path
    ScoreSuspectLists <- NA
    writeLines(paste("FilterSuspectLists = ",as.character(FilterSuspectLists),sep=""),con=file.conn)
  } else {
    FilterSuspectLists <- NA
    ScoreSuspectLists <- NA
  }

  # check if useMoNAMetFusion option is true?
  if (MoNAoffline) {
    if (useMoNAMetFusion) {
      MetFragScoreTypes <- paste(MetFragScoreTypes,",OfflineMetFusionScore",sep="")
      MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
    }
    # check if useMoNAIndiv option is true?
    if (useMonaIndiv) {
      MetFragScoreTypes <- paste(MetFragScoreTypes,",OfflineIndividualMoNAScore",sep="")
      MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
    }
  } else {
    if (useMoNAMetFusion) {
      MetFragScoreTypes <- paste(MetFragScoreTypes,",MetFusionMoNAScore",sep="")
      MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
    }
    # check if useMoNAIndiv option is true?
    if (useMonaIndiv) {
      MetFragScoreTypes <- paste(MetFragScoreTypes,",IndividualMoNAScore",sep="")
      MetFragScoreWeights <- paste(MetFragScoreWeights,"1.0",sep=",")
    }
  }


  ## Pre-processing Inclusion/Exclusion of elements and substructures
  FilterIncludedElements <- incl_el
  FilterExcludedElements <- excl_el
  ElementInclusionExclusiveFilter <- incl_exclusive # default case false as this is very restrictive
  FilterSmartsInclusionList <- incl_smarts_filter
  FilterSmartsExclusionList <- excl_smarts_filter
  SmartsSubstructureInclusionScore <- incl_smarts_score
  SmartsSubstructureExclusionScore <- excl_smarts_score
  #filter by included elements
  if ((nchar(FilterIncludedElements)>0)&&(!ElementInclusionExclusiveFilter)) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",ElementInclusionFilter",sep="")
    writeLines(paste("FilterIncludedElements = ",
                     as.character(FilterIncludedElements),sep=""),con=file.conn)
  } else if ((nchar(FilterIncludedElements)>0)&&(ElementInclusionExclusiveFilter)) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",ElementInclusionExclusiveFilter",sep="")
    writeLines(paste("FilterIncludedElements = ",
                     as.character(FilterIncludedElements),sep=""),con=file.conn)
  }
  # filter by excluded elements
  if (nchar(FilterExcludedElements)>0) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",ElementExclusionFilter",sep="")
    writeLines(paste("FilterExcludedElements = ",
                     as.character(FilterExcludedElements),sep=""),con=file.conn)
  }
  #filter by included SMARTS
  if (nchar(FilterSmartsInclusionList)>0) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",SmartsSubstructureInclusionFilter",sep="")
    writeLines(paste("FilterSmartsInclusionList = ",
                     as.character(FilterSmartsInclusionList),sep=""),con=file.conn)
  }
  #filter by excluded SMARTS
  if (nchar(FilterSmartsExclusionList)>0) {
    MetFragPreProcessingCandidateFilter <- paste(MetFragPreProcessingCandidateFilter,
                                                 ",SmartsSubstructureExclusionFilter",sep="")
    writeLines(paste("FilterSmartsExclusionList = ",
                     as.character(FilterSmartsExclusionList),sep=""),con=file.conn)
  }
  #score by included SMARTS
  if (nchar(SmartsSubstructureInclusionScore)>0) {
    MetFragScoreTypes <- paste(MetFragScoreTypes,",SmartsSubstructureInclusionScore",sep="")
    MetFragScoreWeights <- paste(MetFragScoreWeights,",1.0",sep="")
    writeLines(paste("ScoreSmartsInclusionList = ",as.character(ScoreSmartsInclusionList),sep=""),con=file.conn)
  }
  #score by excluded SMARTS
  if (nchar(SmartsSubstructureExclusionScore)>0) {
    MetFragScoreTypes <- paste(MetFragScoreTypes,",SmartsSubstructureExclusionScore",sep="")
    MetFragScoreWeights <- paste(MetFragScoreWeights,",1",sep="")
    writeLines(paste("ScoreSmartsExclusionList = ",as.character(FilterSmartsExclusionList),sep=""),con=file.conn)
  }
  #
  #user defined scores ...
  #UDS_Category <- ""
  if ((nchar(UDS_Category)>1)) {
#    if ((nchar(UDS_Category)>1)&&(grepl("Local",DB))) {
    MetFragScoreTypes <- paste(MetFragScoreTypes,UDS_Category,sep=",")
    MetFragScoreWeights <- paste(MetFragScoreWeights,UDS_Weights,sep="")
  } #don't need to write anything extra to file as this must be in the localDB


#   #References (only for ChemSpider or ExtendedPubChem)
#   if (MetFragDatabaseType=="ChemSpider" && (nchar(ChemSpiderToken)>0)) {
#     MetFragScoreTypes <- paste(MetFragScoreTypes,
#                                ",ChemSpiderReferenceCount,ChemSpiderDataSourceCount,ChemSpiderNumberPubMedReferences,ChemSpiderRSCCount",
#                                sep="")
#     MetFragScoreWeights <- paste(MetFragScoreWeights,",0.25,0.25,0.25,0.25",sep="")
#   } else if (MetFragDatabaseType=="ExtendedPubChem") {
#     MetFragScoreTypes <- paste(MetFragScoreTypes,
#                                ",PubChemNumberPubMedReferences,PubChemNumberPatents",
#                                sep="")
#     MetFragScoreWeights <- paste(MetFragScoreWeights,",1,1",sep="")
#   }
  writeLines(paste("MetFragPreProcessingCandidateFilter = ",
                   as.character(MetFragPreProcessingCandidateFilter),sep=""),con=file.conn)
  if (filter_by_InChIKey) {
    MetFragPostProcessingCandidateFilter <- "InChIKeyFilter"
    writeLines(paste("MetFragPostProcessingCandidateFilter = ",
                     as.character(MetFragPostProcessingCandidateFilter),sep=""),con=file.conn)
  }
  writeLines(paste("MetFragScoreWeights = ",
                   as.character(MetFragScoreWeights),sep=""),con=file.conn)
  writeLines(paste("MetFragScoreTypes = ",as.character(MetFragScoreTypes),sep=""),con=file.conn)
  close(file.conn)

  return(config_name)

}


####run MetFrag ####

#' Run MetFrag Command Line from config files
#'
#' @description This function runs MetFrag Command Line for the given configuration file and directories.
#' Note that MetFragCL must be available locally to run this. MetFrag Command Line
#' is available from \url{http://c-ruttkies.github.io/MetFrag/projects/metfragcl/}
#'
#' @usage runMetFrag(config_file, MetFrag_dir, CL_name, config_dir=dirname(config_file))
#'
#' @param config_file Full path and file name to the configuration file (as returned by \code{\link{MetFragConfig}})
#' @param MetFrag_dir Full path to the cirectory containing the MetFragCL jar file
#' @param CL_name Name of the exact MetFragCL jar file to use (e.g. \code{MetFrag2.4.4-msready-CL.jar})
#' @param config_dir Full path to the directory in which config files are located. Note that parallel directories
#' \code{log} and \code{results} are created during this process for the first config file run. Defaults to the
#' directory where \code{config_file} is located.
#'
#' @return Runs MetFragCL and creates a log file and, where successful, results files as encoded in the config file.
#' If unsuccessful a status message is printed from the jar; details are saved in the log file.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL author).
#'
#' @seealso \code{\link{MetFragConfig}}
#'
#' @export
#'
#' @examples
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' peaklist_path <- system.file("extdata","EA026206_Simazine_peaks.txt",package="ReSOLUTION")
#' # change this directory to an existing one, or this example won't work
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' testCSV <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
#'
#' config_file <- MetFragConfig(201.0776,"[M+H]+","Simazine_neutralMass",peaklist_path, test_dir, DB="PubChem",neutralPrecursorMass=TRUE)
#' config_file2 <- MetFragConfig(202.0854,1,"Simazine_precMass",peaklist_path,test_dir,DB="LocalCSV",localDB_path=testCSV)
#'
#' #note this first query takes a while to run, try config_file2 for a quicker test.
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#' runMetFrag(config_file2, metfrag_dir, MetFragCL_name)
#'
runMetFrag <- function(config_file, MetFrag_dir, CL_name, config_dir=dirname(config_file)) {
  config_exists <- file.exists(config_file) && file.exists(config_dir)
  current_dir <- getwd()
  if (config_exists) {
    setwd(MetFrag_dir)
    log_dir <- gsub("config","log",config_dir)
    if (!file.exists(log_dir)) {
      dir.create(log_dir)
    }
  } else {
    warning(paste("Configuration file ",config_file,
                  " or directory not found, please try a new file"))
    stop()
  }
  MetFragCommand <- paste("java -Duser.home=",MetFrag_dir, " -jar ", CL_name, " ",config_file,sep="")
  MetFrag_out <- system(command=MetFragCommand,intern=TRUE,show.output.on.console=FALSE)
  log_file <- gsub("config","log",config_file)
  write(MetFrag_out, log_file)
  setwd(current_dir)
}


##### Prepare CompTox Dashboard XLS files for MetFrag #####

#' Prepare CompTox Dashboard MetFrag Export XLS files for MetFragCL
#'
#' @description This function prepares CompTox XLS Export files for use in MetFragCL. It converts the
#' XLS to a CSV and extracts numeric metadata field names for inclusion as scoring terms.
#'
#' @usage CompToxXLStoLocalCSVterms(xls_file, start_index=15, csv_file= "")
#'
#' @param xls_file Full path and file name to the Dashboard Export file to process
#' @param start_index The column number where the metadata columns start (default \code{15} is appropriate
#' for CompTox Dashboard MetFragBeta Export file format released March 2018)
#' @param csv_file If empty, the CSV file will have the same name as \code{xls_file} except the file ending.
#' Enter a file name here if you wish to have a different named file. If a CSV file
#' with the same name exists in the directory it is overwritten.
#'
#' @return Returns a list containing the CSV file name, a list of scoring terms to add to
#' the MetFrag config file and a corresponding score weights entry.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}
#'
#' @export
#'
#' @examples
#' CompToxXLS <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2.xls",package="ReSOLUTION")
#' LocalCSVterms <- CompToxXLStoLocalCSVterms(CompToxXLS)
#'
CompToxXLStoLocalCSVterms <- function(xls_file, start_index=15, csv_file= "") {
  # check csv file, if not defined, make default:
  if (nchar(csv_file)<1) {
    csv_file <- sub(".xls",".csv",xls_file)
  }
  # convert xls to csv
  write.csv(read_excel(xls_file),csv_file,row.names = F)
  cols <- colnames(read_excel(xls_file))
  # read in the content, get colnames and test content
  csv_content <- read.csv(csv_file)
  #cols <- colnames(csv_content)
  # now take a look at the content and run tests
  include_col_i <- vector(mode="numeric",length=0)
  include_col_n <- 1
  # take start_index as the start ... for MS-ready this is 15 (current default)
  for (i in start_index:length(cols)) {
    # this tests if there are ANY numeric values in the column, to avoid errors when running MetFragCL jar
    num_test <- suppressWarnings(length(grep("FALSE",is.na(as.numeric(as.character(csv_content[,i])))))>0)
    if (num_test) {
      include_col_i[include_col_n] <- i
      include_col_n <- include_col_n + 1
    }
  }
  # calculate the score terms and weights
  ScoreTerms <- paste(cols[include_col_i],collapse=",")
  ScoreWeights <- paste0(",", paste(rep(1,length(cols[include_col_i])),collapse=","))
  # generate the output
  LocalCSVterms <- list()
  LocalCSVterms[['CSV']] <- csv_file
  LocalCSVterms[['ScoreTerms']] <- ScoreTerms
  LocalCSVterms[['ScoreWeights']] <- ScoreWeights

  return(LocalCSVterms)

}

##### Prepare CompTox Dashboard CSV Batch files for MetFrag #####

#' Prepare CompTox Dashboard MetFrag Export CSV files for MetFragCL
#'
#' @description This function prepares CompTox CSV Export files for use in MetFragCL. It
#' extracts numeric metadata field names for inclusion as scoring terms.
#'
#' @usage CompToxCSVtoLocalCSVterms(csv_file, start_index=15)
#'
#' @param csv_file Full path and file name to the Dashboard Export file to process
#' @param start_index The column number where the metadata columns start (default \code{15} is appropriate
#' for CompTox Dashboard MetFragBeta Export file format released March 2018)
#'
#' @return Returns a list containing the CSV file name, a list of scoring terms to add to
#' the MetFrag config file and a corresponding score weights entry.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}
#'
#' @export
#'
#' @examples
#' CompToxCSV <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2_wSelectMetaData.csv",package="ReSOLUTION")
#' LocalCSVterms <- CompToxCSVtoLocalCSVterms(CompToxCSV)
#'
CompToxCSVtoLocalCSVterms <- function(csv_file, start_index=15) {
  # read in the content, get colnames and test content
  csv_content <- read.csv(csv_file,check.names=FALSE)
  cols <- colnames(csv_content)
  # now take a look at the content and run tests
  include_col_i <- vector(mode="numeric",length=0)
  include_col_n <- 1
  # take start_index as the start ... for MS-ready this is 15 (current default)
  for (i in start_index:length(cols)) {
    # this tests if there are ANY numeric values in the column, to avoid errors when running MetFragCL jar
    num_test <- suppressWarnings(length(grep("FALSE",is.na(as.numeric(as.character(csv_content[,i])))))>0)
    if (num_test) {
      include_col_i[include_col_n] <- i
      include_col_n <- include_col_n + 1
    }
  }
  # calculate the score terms and weights
  ScoreTerms <- paste(cols[include_col_i],collapse=",")
  ScoreWeights <- paste0(",", paste(rep(1,length(cols[include_col_i])),collapse=","))
  # generate the output
  LocalCSVterms <- list()
  LocalCSVterms[['CSV']] <- csv_file
  LocalCSVterms[['ScoreTerms']] <- ScoreTerms
  LocalCSVterms[['ScoreWeights']] <- ScoreWeights

  return(LocalCSVterms)

}


###### Prepare CompTox Dashboard SDF files for MetFrag #####

#' Prepare CompTox Dashboard MetFrag Export SDF files for MetFragCL
#'
#' @description This function prepares CompTox SDF MetFrag Export files for use in MetFragCL. It reads the
#' SDF and extracts numeric metadata field names for inclusion as scoring terms.
#'
#' @usage CompToxSDFtoLocalSDFterms(SDF_file)
#'
#' @param SDF_file Full path and file name to the Dashboard Export SDF file
#'
#' @return Returns a list containing the SDF file name, a list of scoring terms to add to
#' the MetFrag config file and a corresponding score weights entry.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}, \code{\link{process.sdf.file}}
#'
#' @export
#'
#' @examples
#' CompToxSDF <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2.sdf",package="ReSOLUTION")
#' LocalSDFterms <- CompToxSDFtoLocalSDFterms(CompToxSDF)
#'
#'
CompToxSDFtoLocalSDFterms <- function(SDF_file) {
  # first, process the SDF using Christoph's processSDF function
  SDF_properties <- process.sdf.file(SDF_file)
  SDF_tags <- names(SDF_properties)
  numeric_properties_index <- grep(TRUE, SDF_properties)
  SDF_tags[numeric_properties_index]
  ScoreTerms <- paste(SDF_tags[numeric_properties_index],collapse=",")
  ScoreWeights <- paste0(",",paste(rep(1,length(SDF_tags[numeric_properties_index])),collapse=","))
  LocalSDFterms <- list()
  LocalSDFterms[['SDF']] <- SDF_file
  LocalSDFterms[['ScoreTerms']] <- ScoreTerms
  LocalSDFterms[['ScoreWeights']] <- ScoreWeights

  return(LocalSDFterms)

}


##### Prepare CompTox Dashboard Full CSV file for MetFrag #####

#' Prepare CompTox Dashboard Full CSV file for MetFragCL
#'
#' @description This function extracts metadata headers from the CompTox CSV download file
#' for use in MetFragCL. It reads metadata field names for inclusion as scoring terms.
#'
#' @usage CompToxFullCSVtoLocalCSVterms(csv_file, start_index=13, TermsToRemove="default")
#'
#' @param csv_file Full path and file name to the Dashboard CSV Download file to process
#' @param start_index The column number where the metadata columns start (default \code{13} is appropriate
#' for CompTox Dashboard Download file designed specifically for MetFrag)
#' @param TermsToRemove Define column headers to remove from scoring terms. Select default values 
#' using \code{"default"} (see details) or add a vector of strings. \code{c()}
#'
#' @return Returns a list containing the CSV file name, a list of scoring terms to add to
#' the MetFrag config file and a corresponding score weights entry.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#' 
#' @details The current default \code{TermsToRemove} are 
#' \code{c("TOXCAST_NUMBER_OF_ASSAYS/TOTAL","TOXVAL_Link", "PPRTV_Link", "IRIS_Link")}.
#' This option removes the terms from MetFrag scoring to avoid processing errors, 
#' but these columns are retained in the results file, for downstream use if desired.
#' If MetFrag exits with a status=5, check the log file for terms to add to this list. 
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}
#'
#' @export
#'
#' @examples
#' # note this is an example off a limited file
#' CompToxFullCSVFile <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
#' LocalCSVterms <- CompToxCSVtoLocalCSVterms(CompToxFullCSVFile)
#'
CompToxFullCSVtoLocalCSVterms <- function(csv_file, start_index=13, TermsToRemove="default") {
  # extract column names
  cols <- colnames(read.csv(csv_file,nrows=1,check.names = FALSE))
  # this ensures that only the first row is read and not the entire file, as this is huge
  # as this is designed to work on a specific download file, no numeric test is performed.
  # instead, terms to remove are defined or can be overwritten by the user

  # calculate the score terms and weights
  ScoreTerms <- paste(cols[start_index:length(cols)],collapse=",")
  ScoreWeights <- paste0(",", paste(rep(1,(length(cols)-start_index+1)),collapse=","))
  
  # have to remove columns that do not contain numeric values.
  # use default definition, OR input.
  DefaultTermsTest <- grepl("default",TermsToRemove,fixed=T) 
  DefaultTerms <- DefaultTermsTest[1] && (length(DefaultTermsTest)==1)
  if (DefaultTerms) {
    TermsToRemove <- c("TOXCAST_NUMBER_OF_ASSAYS/TOTAL","TOXVAL_Link", "PPRTV_Link", "IRIS_Link")
  }
  # split out the score terms again
  ScoreTermsSplit <- strsplit(ScoreTerms,",")[[1]]
  ScoreWeightsSplit <- strsplit(ScoreWeights,",")[[1]]
  # calculate indices
  IndicesToRemove <- as.vector(sapply(TermsToRemove,function(string) {grep(string,ScoreTermsSplit,value=F, fixed=T)})) 
  #grep(TermsToRemove,ScoreTerms,value=F,fixed=T)
  ScoreTermsSplit <- ScoreTermsSplit[-(IndicesToRemove)]
  ScoreWeightsSplit <- ScoreWeightsSplit[-(IndicesToRemove)]
  ScoreTerms <- paste(ScoreTermsSplit,collapse=",")
  ScoreWeights <- paste(ScoreWeightsSplit,collapse=",")
  
  # generate the output
  LocalCSVterms <- list()
  LocalCSVterms[['CSV']] <- csv_file
  LocalCSVterms[['ScoreTerms']] <- ScoreTerms
  LocalCSVterms[['ScoreWeights']] <- ScoreWeights

  return(LocalCSVterms)

}




##### Create MetFrag Config Files with CompTox LocalCSV Scoring Terms from MetFrag Export #####


#' Create MetFrag Config Files with LocalCSV and Scoring Terms from CompTox MetFrag XLS Export
#'
#' @description This is a CompTox XLS or CSV specific wrapper function for \code{\link{MetFragConfig}}.
#'
#' @usage MetFragConfig.CompToxCSV(mass, adduct_type, results_filename, peaklist_path, base_dir,
#' CompToxLocalCSVterms, ...)
#'
#' @param mass The mass with which to search the candidate database (\code{DB}). Use \code{neutralPrecursorMass} and
#' \code{adduct_type} to set whether this is monoisotopic mass or an adduct species. Defaults to \code{adduct_type}.
#' @param adduct_type The adduct species used to define mass (if \code{neutralPrecursorMass=FALSE}) and fragmentation settings
#' in the config file, entered as either \code{PrecursorIonType} (text) or \code{PrecursorIonmode} (a number). The available
#' options are given in the system file \code{MetFragAdductTypes.csv} in the \code{extdata} folder.
#' Recommended default values (if ion state is unclear) are \code{[M+H]+} (1) for positive and \code{[M-H]-} (-1) for negative mode.
#' @param results_filename Enter a base filename for naming the results files - do not include file endings
#' @param peaklist_path Enter the full path and file name to the peak list for this config file
#' @param base_dir Enter the directory name to set up the subfolders for MetFrag batch results. If the folders don't exist,
#' subfolders \code{config}, \code{log} and \code{results} are created; the output of this function is saved in \code{config}.
#' @param CompToxLocalCSVterms The output of \code{\link{CompToxXLStoLocalCSVterms}}, \code{\link{CompToxCSVtoLocalCSVterms}}
#' or \code{\link{CompToxFullCSVtoLocalCSVterms}}, used to set
#' \code{DB, localDB_path, UDS_Category and UDS_Weights} in \code{\link{MetFragConfig}}
#'
#' @return Returns a MetFrag config file name
#'
#' @details Remaining parameters are described in \code{\link{MetFragConfig}}
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}, \code{\link{CompToxXLStoLocalCSVterms}},
#' \code{\link{CompToxCSVtoLocalCSVterms}}, \code{\link{CompToxFullCSVtoLocalCSVterms}}
#'
#' @export
#'
#' @examples
#' # Example from DOI: 10.1021/acs.est.7b01908
#' # Note that this scores automatically with all metadata fields in the example file, which
#' # is not necessarily ideal as not all predicted values are relevant for ranking the best candidate.
#' CompToxXLS <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2.xls",package="ReSOLUTION")
#' LocalCSVterms <- CompToxXLStoLocalCSVterms(CompToxXLS)
#' peaklist <- system.file("extdata","EQ300804_Nicotine_peaks.txt",package="ReSOLUTION")
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' config_file <- MetFragConfig.CompToxCSV(163.1230, "[M+H]+","Nicotine_PrecMass_MpHp_XLS",peaklist, test_dir, LocalCSVterms)
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#'
#' # Example of Simazine
#' # Note that this uses a CSV file with fewer scoring terms that are more relevant for candidate selection.
#' CompToxFullCSV_test <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
#' LocalCSVterms <- CompToxFullCSVtoLocalCSVterms(CompToxFullCSV_test)
#' peaklist <- system.file("extdata","EA026206_Simazine_peaks.txt",package="ReSOLUTION")
#' rt_file_path <- system.file("extdata","Eawag_rt_inchi.csv",package="ReSOLUTION")
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' MBrecord <- system.file("extdata","EA026206_Simazine.txt",package="ReSOLUTION")
#' MBinfo <- getMBRecordInfo.MetFragConfig(MBrecord,peaklist, writePeaklist=FALSE)
#' IsPosMode <- grepl(MBinfo$ion_mode,"POSITIVE")
#' adduct_type <- MBinfo$prec_type
#' config_file <- MetFragConfig.CompToxCSV(mass=MBinfo$exact_mass,adduct_type = adduct_type,
#' results_filename = paste0("Simazine", "_byExactMass_5ppm"),peaklist_path = MBinfo$peaklist,neutralPrecursorMass=TRUE,
#' base_dir = test_dir, CompToxLocalCSVterms = LocalCSVterms, IsPosMode = IsPosMode,rt_file_path = rt_file_path,
#' rt_exp = MBinfo$ret_time,filter_by_InChIKey = FALSE)
#'
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#'
#'
#' # Example of diclofenac
#' peaklist <- system.file("extdata","EA020161_Diclofenac_peaks.txt",package="ReSOLUTION")
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' MBinfo <- getMBRecordInfo.MetFragConfig(MBrecord,peaklist, writePeaklist=FALSE)
#' IsPosMode <- grepl(MBinfo$ion_mode,"POSITIVE")
#' adduct_type <- MBinfo$prec_type
#' config_file <- MetFragConfig.CompToxCSV(mass=MBinfo$prec_mass,adduct_type = adduct_type,
#' results_filename = paste0("Diclofenac", "_byPrecMass_5ppm"),peaklist_path = MBinfo$peaklist,neutralPrecursorMass=FALSE,
#' base_dir = test_dir, CompToxLocalCSVterms = LocalCSVterms, IsPosMode = IsPosMode,rt_file_path = rt_file_path,
#' rt_exp = MBinfo$ret_time,filter_by_InChIKey = FALSE)
#'
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#'
MetFragConfig.CompToxCSV <- function(mass, adduct_type, results_filename, peaklist_path, base_dir, CompToxLocalCSVterms,
                                           output="XLS", ppm=5, mzabs=0.001, frag_ppm=5, IsPosMode=TRUE, tree_depth=2, num_threads=1,
                                           add_refs=FALSE, minInt=0, rt_file_path="",rt_exp=0, suspect_path="", suspect_filter=FALSE,
                                           #DB=c("LocalCSV","LocalSDF"), localDB_path="", UDS_Category="",UDS_Weights=""
                                           token="", DB_IDs="",mol_form="",useFormula=FALSE,neutralPrecursorMass=FALSE,
                                           useMoNAMetFusion=TRUE,useMonaIndiv=TRUE,MoNAoffline=TRUE,
                                           incl_el="", excl_el="", incl_exclusive=FALSE,
                                           incl_smarts_filter="",incl_smarts_score="", excl_smarts_filter="",excl_smarts_score="",
                                           filter_isotopes=TRUE,filter_by_InChIKey=TRUE) {
  #define parameters that are missing
  #DB="localCSV"
  #localDB_path=CompToxLocalCSVterms$CSV
  #UDS_Category=CompToxLocalCSVterms$ScoreTerms
  #UDS_Weights=CompToxLocalCSVterms$ScoreWeights
  config_file <- MetFragConfig(mass=mass, adduct_type=adduct_type, results_filename=results_filename, peaklist_path=peaklist_path,
                               base_dir=base_dir, DB="LocalCSV", localDB_path = CompToxLocalCSVterms$CSV,
                               UDS_Category=CompToxLocalCSVterms$ScoreTerms,
                               UDS_Weights=CompToxLocalCSVterms$ScoreWeights,
                               output=output, ppm=ppm, mzabs=mzabs, frag_ppm=frag_ppm, IsPosMode=IsPosMode,
                               tree_depth=tree_depth, num_threads=num_threads,
                               add_refs=add_refs, minInt=minInt, rt_file_path=rt_file_path, rt_exp=rt_exp, suspect_path=suspect_path,
                               suspect_filter=suspect_filter, token=token, DB_IDs=DB_IDs,mol_form=mol_form,useFormula=useFormula,
                               neutralPrecursorMass=neutralPrecursorMass,
                               useMoNAMetFusion=useMoNAMetFusion, useMonaIndiv=useMonaIndiv, MoNAoffline=MoNAoffline,
                               incl_el=incl_el, excl_el=excl_el, incl_exclusive=incl_exclusive,
                               incl_smarts_filter=incl_smarts_filter, incl_smarts_score=incl_smarts_score,
                               excl_smarts_filter=excl_smarts_filter, excl_smarts_score=excl_smarts_score,
                               filter_isotopes=filter_isotopes,filter_by_InChIKey=filter_by_InChIKey)

  return(config_file)
}



##### Create MetFrag Config Files with CompTox LocalSDF Scoring Terms from MetFrag Export #####


#' Create MetFrag Config Files with LocalSDF and Scoring Terms from CompTox MetFrag SDF Export
#'
#' @description This is a CompTox SDF Export-specific wrapper function for \code{\link{MetFragConfig}}.
#'
#' @usage MetFragConfig.CompToxExportSDF(mass, adduct_type, results_filename, peaklist_path, base_dir,
#' CompToxSDFtoLocalSDFterms, ...)
#'
#' @param mass The mass with which to search the candidate database (\code{DB}). Use \code{neutralPrecursorMass} and
#' \code{adduct_type} to set whether this is monoisotopic mass or an adduct species. Defaults to \code{adduct_type}.
#' @param adduct_type The adduct species used to define mass (if \code{neutralPrecursorMass=FALSE}) and fragmentation settings
#' in the config file, entered as either \code{PrecursorIonType} (text) or \code{PrecursorIonmode} (a number). The available
#' options are given in the system file \code{MetFragAdductTypes.csv} in the \code{extdata} folder.
#' Recommended default values (if ion state is unclear) are \code{[M+H]+} (1) for positive and \code{[M-H]-} (-1) for negative mode.
#' @param results_filename Enter a base filename for naming the results files - do not include file endings
#' @param peaklist_path Enter the full path and file name to the peak list for this config file
#' @param base_dir Enter the directory name to set up the subfolders for MetFrag batch results. If the folders don't exist,
#' subfolders \code{config}, \code{log} and \code{results} are created; the output of this function is saved in \code{config}.
#' @param CompToxSDFtoLocalSDFterms The output of \code{\link{CompToxSDFtoLocalSDFterms}}, used to set
#' \code{DB, localDB_path, UDS_Category and UDS_Weights} in \code{\link{MetFragConfig}}
#'
#' @return Returns a MetFrag config file name
#'
#' @details Remaining parameters are described in \code{\link{MetFragConfig}}
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}, \code{\link{CompToxSDFtoLocalSDFterms}}
#'
#' @export
#'
#' @examples
#' # Example from DOI: 10.1021/acs.est.7b01908
#' CompToxSDF <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2.sdf",package="ReSOLUTION")
#' LocalSDFterms <- CompToxSDFtoLocalSDFterms(CompToxSDF)
#' peaklist <- system.file("extdata","EQ300804_Nicotine_peaks.txt",package="ReSOLUTION")
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' config_file <- MetFragConfig.CompToxExportSDF(163.1230, "[M+H]+","Nicotine_PrecMass_MpHp_SDF",peaklist, test_dir, LocalSDFterms)
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#'
MetFragConfig.CompToxExportSDF <- function(mass, adduct_type, results_filename, peaklist_path, base_dir, CompToxSDFtoLocalSDFterms,
                                           output="XLS", ppm=5, mzabs=0.001, frag_ppm=5, IsPosMode=TRUE, tree_depth=2, num_threads=1,
                                           add_refs=FALSE, minInt=0, rt_file_path="",rt_exp=0, suspect_path="", suspect_filter=FALSE,
                                           #DB=c("LocalCSV","LocalSDF"), localDB_path="", UDS_Category="",UDS_Weights=""
                                           token="", DB_IDs="",mol_form="",useFormula=FALSE,neutralPrecursorMass=FALSE,
                                           useMoNAMetFusion=TRUE,useMonaIndiv=TRUE,MoNAoffline=TRUE,
                                           incl_el="", excl_el="", incl_exclusive=FALSE,
                                           incl_smarts_filter="",incl_smarts_score="", excl_smarts_filter="",excl_smarts_score="",
                                           filter_isotopes=TRUE,filter_by_InChIKey=TRUE) {
  #define parameters that are missing
  #DB="localSDF"
  #localDB_path=CompToxSDFtoLocalSDFterms$SDF
  #UDS_Category=CompToxSDFtoLocalSDFterms$ScoreTerms
  #UDS_Weights=CompToxSDFtoLocalSDFterms$ScoreWeights
  config_file <- MetFragConfig(mass=mass, adduct_type=adduct_type, results_filename=results_filename, peaklist_path=peaklist_path,
                               base_dir=base_dir, DB="LocalSDF", localDB_path = CompToxSDFtoLocalSDFterms$SDF,
                               UDS_Category=CompToxSDFtoLocalSDFterms$ScoreTerms,
                               UDS_Weights=CompToxSDFtoLocalSDFterms$ScoreWeights,
                               output=output, ppm=ppm, mzabs=mzabs, frag_ppm=frag_ppm, IsPosMode=IsPosMode,
                               tree_depth=tree_depth, num_threads=num_threads,
                               add_refs=add_refs, minInt=minInt, rt_file_path=rt_file_path, rt_exp=rt_exp, suspect_path=suspect_path,
                               suspect_filter=suspect_filter, token=token, DB_IDs=DB_IDs,mol_form=mol_form,useFormula=useFormula,
                               neutralPrecursorMass=neutralPrecursorMass,
                               useMoNAMetFusion=useMoNAMetFusion, useMonaIndiv=useMonaIndiv, MoNAoffline=MoNAoffline,
                               incl_el=incl_el, excl_el=excl_el, incl_exclusive=incl_exclusive,
                               incl_smarts_filter=incl_smarts_filter, incl_smarts_score=incl_smarts_score,
                               excl_smarts_filter=excl_smarts_filter, excl_smarts_score=excl_smarts_score,
                               filter_isotopes=filter_isotopes,filter_by_InChIKey=filter_by_InChIKey)

  return(config_file)
}



# ##### Create MetFrag Config Files with CompTox LocalCSV Scoring Terms from Full CSV Download File #####
# ## Superceded, generalised the CSV function ... no need for specific one... as they use the same localCSVterms
#
# #' Create MetFrag Config Files with LocalCSV and Scoring Terms from CompTox Full CSV Download File
# #'
# #' @description This is a CompTox CSV Download-specific wrapper function for \code{\link{MetFragConfig}}.
# #'
# #' @usage MetFragConfig.CompToxFullCSV(mass, adduct_type, results_filename, peaklist_path, base_dir,
# #' CompToxFullCSVtoLocalCSVterms, ...)
# #'
# #' @param mass The mass with which to search the candidate database (\code{DB}). Use \code{neutralPrecursorMass} and
# #' \code{adduct_type} to set whether this is monoisotopic mass or an adduct species. Defaults to \code{adduct_type}.
# #' @param adduct_type The adduct species used to define mass (if \code{neutralPrecursorMass=FALSE}) and fragmentation settings
# #' in the config file, entered as either \code{PrecursorIonType} (text) or \code{PrecursorIonmode} (a number). The available
# #' options are given in the system file \code{MetFragAdductTypes.csv} in the \code{extdata} folder.
# #' Recommended default values (if ion state is unclear) are \code{[M+H]+} (1) for positive and \code{[M-H]-} (-1) for negative mode.
# #' @param results_filename Enter a base filename for naming the results files - do not include file endings
# #' @param peaklist_path Enter the full path and file name to the peak list for this config file
# #' @param base_dir Enter the directory name to set up the subfolders for MetFrag batch results. If the folders don't exist,
# #' subfolders \code{config}, \code{log} and \code{results} are created; the output of this function is saved in \code{config}.
# #' @param CompToxFullCSVtoLocalCSVterms The output of \code{\link{CompToxFullCSVtoLocalCSVterms}}, used to set
# #' \code{DB, localDB_path, UDS_Category and UDS_Weights} in \code{\link{MetFragConfig}}
# #'
# #' @return Returns a MetFrag config file name
# #'
# #' @details Remaining parameters are described in \code{\link{MetFragConfig}}. PROTOTYPE AT THIS STAGE
# #'
# #' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
# #' Antony J. Williams and team (CompTox Dashboard)
# #'
# #' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}, \code{\link{CompToxFullCSVtoLocalCSVterms}}
# #'
# #' @export
# #'
# #' @examples
# #' # Example of Simazine
# #' CompToxFullCSV_test <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
# #' LocalCSVterms <- CompToxFullCSVtoLocalCSVterms(CompToxFullCSV_test)
# #' peaklist <- system.file("extdata","EA026206_Simazine_peaks.txt",package="ReSOLUTION")
# #' rt_file_path <- system.file("extdata","Eawag_rt_inchi.csv",package="ReSOLUTION")
# #' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
# #' MBrecord <- system.file("extdata","EA026206_Simazine.txt",package="ReSOLUTION")
# #' MBinfo <- getMBRecordInfo.MetFragConfig(MBrecord,peaklist, writePeaklist=FALSE)
# #' IsPosMode <- grepl(MBinfo$ion_mode,"POSITIVE")
# #' adduct_type <- MBinfo$prec_type
# #' config_file <- MetFragConfig.CompToxFullCSV(mass=MBinfo$exact_mass,adduct_type = adduct_type,
# #' results_filename = paste0("Simazine", "_byExactMass_5ppm"),peaklist_path = MBinfo$peaklist,neutralPrecursorMass=TRUE,
# #' base_dir = test_dir, CompToxFullCSVtoLocalCSVterms = LocalCSVterms, IsPosMode = IsPosMode,rt_file_path = rt_file_path,
# #' rt_exp = MBinfo$ret_time,filter_by_InChIKey = FALSE)
# #'
# #' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
# #' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
# #' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
# #'
# #'
# #' # Example of diclofenac
# #' peaklist <- system.file("extdata","EA020161_Diclofenac_peaks.txt",package="ReSOLUTION")
# #' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
# #' MBinfo <- getMBRecordInfo.MetFragConfig(MBrecord,peaklist, writePeaklist=FALSE)
# #' IsPosMode <- grepl(MBinfo$ion_mode,"POSITIVE")
# #' adduct_type <- MBinfo$prec_type
# #' config_file <- MetFragConfig.CompToxFullCSV(mass=MBinfo$prec_mass,adduct_type = adduct_type,
# #' results_filename = paste0("Diclofenac", "_byPrecMass_5ppm"),peaklist_path = MBinfo$peaklist,neutralPrecursorMass=FALSE,
# #' base_dir = test_dir, CompToxFullCSVtoLocalCSVterms = LocalCSVterms, IsPosMode = IsPosMode,rt_file_path = rt_file_path,
# #' rt_exp = MBinfo$ret_time,filter_by_InChIKey = FALSE)
# #'
# #' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
# #' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
# #' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
# #'
# #'
# #'
# MetFragConfig.CompToxFullCSV <- function(mass, adduct_type, results_filename, peaklist_path, base_dir, CompToxFullCSVtoLocalCSVterms,
#                                            output="XLS", ppm=5, mzabs=0.001, frag_ppm=5, IsPosMode=TRUE, tree_depth=2, num_threads=1,
#                                            add_refs=FALSE, minInt=0, rt_file_path="",rt_exp=0, suspect_path="", suspect_filter=FALSE,
#                                            #DB=c("LocalCSV","LocalSDF"), localDB_path="", UDS_Category="",UDS_Weights=""
#                                            token="", DB_IDs="",mol_form="",useFormula=FALSE,neutralPrecursorMass=FALSE,
#                                            useMoNAMetFusion=TRUE,useMonaIndiv=TRUE,MoNAoffline=TRUE,
#                                            incl_el="", excl_el="", incl_exclusive=FALSE,
#                                            incl_smarts_filter="",incl_smarts_score="", excl_smarts_filter="",excl_smarts_score="",
#                                            filter_isotopes=TRUE,filter_by_InChIKey=TRUE) {
#   #define parameters that are missing
#   #DB="localCSV"
#   #localDB_path=CompToxXLStoLocalCSVterms$CSV
#   #UDS_Category=CompToxXLStoLocalCSVterms$ScoreTerms
#   #UDS_Weights=CompToxXLStoLocalCSVterms$ScoreWeights
#   config_file <- MetFragConfig(mass=mass, adduct_type=adduct_type, results_filename=results_filename, peaklist_path=peaklist_path,
#                                base_dir=base_dir, DB="LocalCSV", localDB_path = CompToxFullCSVtoLocalCSVterms$CSV,
#                                UDS_Category=CompToxFullCSVtoLocalCSVterms$ScoreTerms,
#                                UDS_Weights=CompToxFullCSVtoLocalCSVterms$ScoreWeights,
#                                output=output, ppm=ppm, mzabs=mzabs, frag_ppm=frag_ppm, IsPosMode=IsPosMode,
#                                tree_depth=tree_depth, num_threads=num_threads,
#                                add_refs=add_refs, minInt=minInt, rt_file_path=rt_file_path, rt_exp=rt_exp, suspect_path=suspect_path,
#                                suspect_filter=suspect_filter, token=token, DB_IDs=DB_IDs,mol_form=mol_form,useFormula=useFormula,
#                                neutralPrecursorMass=neutralPrecursorMass,
#                                useMoNAMetFusion=useMoNAMetFusion, useMonaIndiv=useMonaIndiv, MoNAoffline=MoNAoffline,
#                                incl_el=incl_el, excl_el=excl_el, incl_exclusive=incl_exclusive,
#                                incl_smarts_filter=incl_smarts_filter, incl_smarts_score=incl_smarts_score,
#                                excl_smarts_filter=excl_smarts_filter, excl_smarts_score=excl_smarts_score,
#                                filter_isotopes=filter_isotopes,filter_by_InChIKey=filter_by_InChIKey)
#
#   return(config_file)
# }


##### Create MetFrag Config Files with Formula not Mass #####


#' Create MetFrag Config Files with Molecular Formula Candidate Searching
#'
#' @description This is a molecular formula-specific wrapper function for \code{\link{MetFragConfig}}.
#'
#' @usage MetFragConfig.formula(mol_form, adduct_type, IsPosMode, results_filename, peaklist_path, base_dir, ...)
#'
#' @param mol_form The molecular formula (as a string) with which to search the candidate database (\code{DB}).
#' \code{mass} is calculated automatically.
#' @param adduct_type The adduct species used to define mass (if \code{neutralPrecursorMass=FALSE}) and fragmentation settings
#' in the config file, entered as either \code{PrecursorIonType} (text) or \code{PrecursorIonmode} (a number). The available
#' options are given in the system file \code{MetFragAdductTypes.csv} in the \code{extdata} folder.
#' Recommended default values (if ion state is unclear) are \code{[M+H]+} (1) for positive and \code{[M-H]-} (-1) for negative mode.
#' @param isPosMode Use \code{TRUE} or \code{FALSE} to set positive or negative ionization as appropriate.
#' @param results_filename Enter a base filename for naming the results files - do not include file endings
#' @param peaklist_path Enter the full path and file name to the peak list for this config file
#' @param base_dir Enter the directory name to set up the subfolders for MetFrag batch results. If the folders don't exist,
#' subfolders \code{config}, \code{log} and \code{results} are created; the output of this function is saved in \code{config}.
#'
#' @return Returns a MetFrag config file name
#'
#' @details Remaining parameters are described in \code{\link{MetFragConfig}}
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu> in partnership with Christoph Ruttkies (MetFragCL),
#' Antony J. Williams and team (CompTox Dashboard)
#'
#' @seealso \code{\link{MetFragConfig}}, \code{\link{runMetFrag}}, \code{\link{getAdductMassesFromFormula}}
#'
#' @export
#'
#' @examples
#' # Example of Simazine
#' peaklist_path <- system.file("extdata","EA026206_Simazine_peaks.txt",package="ReSOLUTION")
#' # change this directory to an existing one, or this example won't work
#' test_dir <- "C:/DATA/Workflow/MetFrag22/metfrag_test_results"
#' testCSV <- system.file("extdata","dsstox_MS_Ready_MetFragTestCSV5.csv",package="ReSOLUTION")
#'
#' config_file <- MetFragConfig.formula("C7H12ClN5","[M+H]+",IsPosMode=TRUE, "Simazine_formula_PubChem",peaklist_path, test_dir,
#' DB="PubChem")
#' config_file2 <- MetFragConfig.formula("C7H12ClN5",1,IsPosMode=TRUE, "Simazine_formula_CompTox",peaklist_path,test_dir,
#' DB="LocalCSV",localDB_path=testCSV)
#' config_file3 <- MetFragConfig.formula("C7H12ClN5",1,IsPosMode=TRUE, "Simazine_formula_CompTox_noInChIFilter",peaklist_path,test_dir,
#' DB="LocalCSV",localDB_path=testCSV,filter_by_InChIKey = FALSE)
#'
#' metfrag_dir <- "C:/DATA/Workflow/MetFrag22/"
#' MetFragCL_name <- "MetFrag2.4.4-msready-CL.jar"
#' #note this first query is a longer query using PubChem
#' runMetFrag(config_file, metfrag_dir, MetFragCL_name)
#' runMetFrag(config_file2, metfrag_dir, MetFragCL_name)
#' runMetFrag(config_file3, metfrag_dir, MetFragCL_name)
#'
MetFragConfig.formula <- function(mol_form, adduct_type, IsPosMode, results_filename, peaklist_path, base_dir,
                                  DB=c("KEGG","PubChem","ExtendedPubChem","ChemSpider","FOR-IDENT","MetaCyc",
                                       "LocalCSV","LocalPSV","LocalSDF"), localDB_path="",
                                  output="XLS", ppm=5, mzabs=0.001, frag_ppm=5, tree_depth=2, num_threads=1,
                                  add_refs=FALSE, minInt=0, rt_file_path="",rt_exp=0, suspect_path="", suspect_filter=FALSE,
                                  UDS_Category="",UDS_Weights="", token="", DB_IDs="",neutralPrecursorMass=TRUE,
                                  useMoNAMetFusion=TRUE,useMonaIndiv=TRUE,MoNAoffline=TRUE,
                                  incl_el="", excl_el="", incl_exclusive=FALSE,
                                  incl_smarts_filter="",incl_smarts_score="", excl_smarts_filter="",excl_smarts_score="",
                                  filter_isotopes=TRUE,filter_by_InChIKey=TRUE) {
  #define parameters that are missing
  #useFormula=TRUE
  #mass=getAdductMassesFromFormula(mol_form)$Monoiso_mass
  mass <- getAdductMassesFromFormula(mol_form)$Monoiso_mass
  useFormula <- TRUE
  # hard program this as we fix it in the wrapper
  neutralPrecursorMass=TRUE
  config_file <- MetFragConfig(mass=mass, adduct_type=adduct_type, results_filename=results_filename, peaklist_path=peaklist_path,
                               base_dir=base_dir, DB=DB, localDB_path = localDB_path, mol_form=mol_form, useFormula=useFormula,
                               output=output, ppm=ppm, mzabs=mzabs, frag_ppm=frag_ppm, IsPosMode=IsPosMode,
                               tree_depth=tree_depth, num_threads=num_threads,
                               add_refs=add_refs, minInt=minInt, rt_file_path=rt_file_path, rt_exp=rt_exp, suspect_path=suspect_path,
                               suspect_filter=suspect_filter, token=token, DB_IDs=DB_IDs,
                               neutralPrecursorMass=neutralPrecursorMass,
                               UDS_Category=UDS_Category,UDS_Weights=UDS_Weights,
                               useMoNAMetFusion=useMoNAMetFusion, useMonaIndiv=useMonaIndiv, MoNAoffline=MoNAoffline,
                               incl_el=incl_el, excl_el=excl_el, incl_exclusive=incl_exclusive,
                               incl_smarts_filter=incl_smarts_filter, incl_smarts_score=incl_smarts_score,
                               excl_smarts_filter=excl_smarts_filter, excl_smarts_score=excl_smarts_score,
                               filter_isotopes=filter_isotopes,filter_by_InChIKey=filter_by_InChIKey)

  return(config_file)
}


###### Extract Info from MassBank records for MetFrag Config files #####

#' Extract Info from MassBank records for MetFrag Config files
#'
#' @description This is a convenience wrapper to extract info for MetFrag Config files if source mass spectra
#' are already in MassBank format.
#'
#' @usage getMBRecordInfo.MetFragConfig(MBrecord,peaklist_filepath,writePeaklist=TRUE)
#'
#' @param MBrecord Full path to the MassBank record containing peaks of interest
#' @param peaklist_filepath Full path to the peaklist that will be used in the MetFrag Config file. It is
#' recommeded to use a temporary file.
#' @param writePeaklist Default \code{TRUE} to extract and write peaklists to a temporary file. To use an existing
#' file, set to \code{FALSE}.
#'
#' @return Returns a list containing the peaklist filename, molecular formula, exact mass, precursor type,
#' ion mode and retention time
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getMBRecordPeaks}}, \code{\link{getMBRecordEntry}}, \code{\link{MetFragConfig}}
#'
#' @export
#'
#' @examples
#'
#' peaklist <- system.file("extdata","EA020161_Diclofenac_peaks.txt",package="ReSOLUTION")
#' MBrecord <- system.file("extdata","EA020161_Diclofenac.txt",package="ReSOLUTION")
#' MBinfo <- getMBRecordInfo.MetFragConfig(MBrecord,peaklist, writePeaklist=FALSE)
#' # use writePeaklist=TRUE to create a peaklist from the MassBank record
#'
getMBRecordInfo.MetFragConfig <- function(MBrecord,peaklist_filepath,writePeaklist=TRUE) {
  # get the data
  if (writePeaklist) {
    write.table(x = getMBRecordPeaks(MBrecord), file = peaklist_filepath,quote = F,row.names=F,col.names = F)
  }
  recVec <- MBFileToVector(MBrecord)
  cmpd_Name <- getMBRecordEntry("CH$NAME:", recVec)[1] # returns the first name only
  cmpd_SMILES <- getMBRecordEntry("CH$SMILES:", recVec)
  cmpd_InChIKey <- getMBRecordEntry("CH$LINK: INCHIKEY", recVec)
  mol_form <- getMBRecordEntry("CH$FORMULA:", recVec)
  exact_mass <- getMBRecordEntry("CH$EXACT_MASS", recVec)
  prec_mass <- getMBRecordEntry("MS$FOCUSED_ION: PRECURSOR_M/Z", recVec)
  prec_type <- getMBRecordEntry("MS$FOCUSED_ION: PRECURSOR_TYPE",recVec)
  ion_mode <- getMBRecordEntry("AC$MASS_SPECTROMETRY: ION_MODE",recVec)
  ret_time <- strsplit(getMBRecordEntry("AC$CHROMATOGRAPHY: RETENTION_TIME",recVec)," ",fixed=TRUE)[[1]][1]
  ret_time_unit <- strsplit(getMBRecordEntry("AC$CHROMATOGRAPHY: RETENTION_TIME",recVec)," ",fixed=TRUE)[[1]][2]
  # export the data
  MetFragMBInfo <- list()
  MetFragMBInfo[['source_file']] <- MBrecord
  MetFragMBInfo[['peaklist']] <- peaklist_filepath
  MetFragMBInfo[['cmpd_Name']] <- cmpd_Name
  MetFragMBInfo[['cmpd_SMILES']] <- cmpd_SMILES
  MetFragMBInfo[['cmpd_InChIKey']] <- cmpd_InChIKey
  MetFragMBInfo[['mol_form']] <- mol_form
  MetFragMBInfo[['exact_mass']] <- exact_mass
  MetFragMBInfo[['prec_mass']] <- prec_mass
  MetFragMBInfo[['prec_type']] <- prec_type
  MetFragMBInfo[['ion_mode']] <- ion_mode
  MetFragMBInfo[['ret_time']] <- ret_time
  MetFragMBInfo[['ret_time_unit']] <- ret_time_unit

  return(MetFragMBInfo)

}


##### Fill in missing InChIKeys in CSV files ####


#' Fill in Missing InChIKeys in CSV files
#'
#' @description This is a small wrapper function to patch CSV files that are missing some
#' InChIKey entries, using Cactus then Open Babel to fill gaps from SMILES.
#'
#' @usage addMissingInChIKeys.CSV(csv_file, inchikey_col_number, smiles_col_number, babel_dir, write=TRUE, csv_file_out="")
#'
#' @param csv_file Full path and file name to CSV file to fill in InChIKey gaps
#' @param inchikey_col_number Column number of the column containing InChIKeys to check/fill
#' @param smiles_col_number Column number of the SMILES code needed to calculate missing InChIKeys
#' @param babel_dir Path to the directory containing \code{obabel.exe}, e.g. \code{"C:/Program Files (x86)/OpenBabel-2.3.2"}
#' @param write If \code{TRUE}, writes either to \code{csv_file_out} (of defined) or overwrites \code{csv_file} and
#' returns the file name. If \code{FALSE}, returns a list of InChIKeys and does not write to file.
#' @param csv_file_out Full path and file name of output CSV, define to avoid overwriting the original file.
#'
#' @return Returns path and file name to CSV file with fixed entries (if necessary) or list of InChIKeys.
#'
#' @author Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @seealso \code{\link{getInChIKey.obabel}}, \code{\link{getCactus}}, \code{\link{InChIKey_test}}
#'
#' @export
#'
#' @examples
#' testCSV <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2_wSelectMetaData.csv",package="ReSOLUTION")
#' babel_dir <- "C:/Program Files (x86)/OpenBabel-2.3.2"
#' InChIKeys <- addMissingInChIKeys.CSV(testCSV,13,5,babel_dir,write=FALSE)
#'
addMissingInChIKeys.CSV <- function(csv_file, inchikey_col_number, smiles_col_number, babel_dir, write=TRUE, csv_file_out="") {
  csv_content <- read.csv(csv_file,check.names=FALSE,colClasses = "character")
  cols <- colnames(csv_content)
  InChIKeys <- csv_content[,inchikey_col_number]
  # find out whether any InChIKeys are missing
  InChIKeyCheck <- sapply(InChIKeys, InChIKey_test)
  MissingInChIKey_index <- as.vector(which(InChIKeyCheck==FALSE))
  csv_file_name <- csv_file
  if (length(MissingInChIKey_index<1)) {
    print("No InChIKeys missing; no changes to file")
  } else {
    for (i in MissingInChIKey_index) {
      SMILES <- csv_content[i,smiles_col_number]
      InChIKey <- trimKey(getCactus(SMILES,stdinchikey))
      if (InChIKey_test(InChIKey)) {
        csv_content[i,inchikey_col_number] <- InChIKey
      } else {
        InChIKey <- getInChIKey.obabel(SMILES,babel_dir)
      }
      if (InChIKey_test(InChIKey)) {
        csv_content[i,inchikey_col_number] <- InChIKey
      } else {
        print("No InChIKey found using Cactus or OpenBabel")
      }
      if (write) {
        if (nchar(csv_file_out)>3) {
          write.csv(csv_content,csv_file_out,row.names=FALSE)
          print("New CSV file written with InChIKeys filled in ")
          csv_file_name <- csv_file_out
        } else {
          write.csv(csv_content,csv_file,row.names=FALSE)
          print("Input CSV file overwritten with InChIKeys filled in")
        }
      } else {
        print("Missing InChIKeys calculated but no output requested; returning InChIKey list")
      }
    }
  }
  if (write) {
    return(csv_file_name)
  } else {
    InChIKeys <- csv_content[,inchikey_col_number]
    return(InChIKeys)
  }
}
