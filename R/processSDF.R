# from Christoph Ruttkies, 05 April 2018.
# modified Emma Schymanski, 06 April 2018.
# documented Emma Schymanski, 18 April 2018


#' Process SDF File to retrieve tags and properties
#'
#' @description This parses an SDF file and retrieves information used to select relevant to MetFrag
#' scoring terms from CompTox Dashboard Export files. Used within \code{\link{CompToxSDFtoLocalSDFterms}}.
#'
#' @usage process.sdf.file(SDF_file)
#'
#' @param SDF_file Full path and file name to SDF file
#'
#' @author Christoph Ruttkies, Emma Schymanski <emma.schymanski@@uni.lu>
#'
#' @return This returns a set of SDF tags and information whether values are numeric or not
#'
#' @seealso \code{\link{CompToxSDFtoLocalSDFterms}}
#'
#' @export
#'
#' @examples
#' CompToxSDF <- system.file("extdata","CompToxBatchSearch_MetFrag_MSready_C10H14N2.sdf",package="ReSOLUTION")
#' process.sdf.file(CompToxSDF)
#'
process.sdf.file = function(SDF_file) {
  # a list to store possible scoring property names
  properties_to_number <- list()
  # test if file exists ...
  if(!file.exists(SDF_file))return(properties_to_number)
  # if it exists... open connection to SDF file
  con = file(SDF_file, "r")
  # store last found property name
  current.property.name = NULL
  while ( TRUE ) {
    line = readLines(con, n = 1)  # read next line
    if (length(line) == 0) {      # stop if end of file
      break
    }
    line = gsub("\\s*$", "", gsub("^\\s*", "", line))                      # remove starting/trailing spaces
    if(length(grep("^>", line)) > 0) {                                     # check for new property name
      current.property.name = gsub(">\\s*<(.*)>\\s*$", "\\1", line)        # store it if found
    } else if(!is.null(current.property.name)) {                           # if property name was found in previous round
      # check if current property was not yet found or is set to FALSE
      if(is.null(properties_to_number[[current.property.name]]) || !properties_to_number[[current.property.name]]) {
        properties_to_number[[current.property.name]] = !suppressWarnings(is.na(as.numeric(line))) # check for numeric

      }
      current.property.name = NULL
    }
  }
  close(con)

  # unset prohibited property values (not used for scoring)
  prohibited.properties = c("CompoundName", "PREFERRED_NAME_DTXSID", "Identifier", "DTXSID", "MAPPED_DTXSID", "InChIKey", "INCHIKEY_DTXCID", "InChIKey1", "InChIKey2", "InChIKey3", "InChI", "INCHI_STRING_DTXCID", "IUPACName", "MolecularFormula", "FORMULA_INDIVIDUAL_COMPONENT", "MonoisotopicMass", "MONOISOTOPIC_MASS_DTXCID", "XlogP3", "SMILES", "SMILES_INDIVIDUAL_COMPONENT", "MaximumTreeDepth", "CHEMSPIDER_ALOGP", "CHEMSPIDER_XLOGP", "Score", "NoExplPeaks", "ExplPeaks", "NumberPeaksUsed", "FormulasOfExplPeaks", "INCHI STRING", "MONOSIOTOPIC MASS", "MOL FORMULA", "MOLECULAR FORMULA")
  for (i in 1:length(properties_to_number)) {
    if(any(prohibited.properties == names(properties_to_number)[i])) properties_to_number[[names(properties_to_number)[i]]] = F
  }
  return(properties_to_number)
}
