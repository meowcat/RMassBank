# Script for writing MassBank files

#testtest change
#' Load MassBank compound information lists
#' 
#' Loads MassBank compound information lists (i.e. the lists which were created
#' in the first two steps of the MassBank \code{\link{mbWorkflow}} and
#' subsequently edited by hand.).
#' 
#' \code{resetInfolists} clears the information lists, i.e. it creates a new
#' empty list in \code{mbdata_archive}. \code{loadInfolist} loads a single CSV
#' file, whereas \code{loadInfolists} loads a whole directory.
#' 
#' @aliases loadInfolists loadInfolist resetInfolists
#' @usage loadInfolists(mb, path)
#' 
#'  loadInfolist(mb, fileName)
#' 
#'  resetInfolists(mb)
#' @param path Directory in which the namelists reside. All CSV files in this
#' directory will be loaded.
#' @param fileName A single namelist to be loaded.
#' @param mb The \code{mbWorkspace} to load/reset the lists in.
#' @return The new workspace with loaded/reset lists.
#' @author Michael Stravs
#' @examples
#' 
#' #
#' \dontrun{mb <- resetInfolists(mb)
#' 	mb <- loadInfolist(mb, "my_csv_infolist.csv")}
#' 
#' @export
loadInfolists <- function(mb, path)
{
  archivefiles <- list.files(path, ".csv", full.names=TRUE)
  for(afile in archivefiles)
    mb <- loadInfolist(mb, afile)
  return(mb)
}

# Load an "infolist". This loads a CSV file which should contain the entries
# edited and controlled by hand. All compound infos from fileName are added into the
# global mbdata_archive. Entries with a cpdID which was already present, are substituted
# by new entries from the fileName file.
#' @export
loadInfolist <- function(mb, fileName)
{
  # Prime a new infolist if it doesn't exist
  if(ncol(mb@mbdata_archive) == 0)
    mb <- resetInfolists(mb)
  mbdata_new <- read.csv(fileName, sep=",", stringsAsFactors=FALSE)
  # Legacy check for loading the Uchem format files.
  # Even if dbname_* are not used downstream of here, it's still good to keep them
  # for debugging reasons.
  n <- colnames(mbdata_new)
  cols <- c("id","dbcas","dataused")
  
  # Check if comma-separated or semicolon-separated
  d <- setdiff(cols, n)
  if(length(d)>0){
		mbdata_new <- read.csv2(fileName, stringsAsFactors=FALSE)
		n <- colnames(mbdata_new)
		d2 <- setdiff(cols, n)
		if(length(d2) > 0){
			stop("Some columns are missing in the infolist.")
		}
	}
  if("dbname_d" %in% colnames(mbdata_new))
  {
    colnames(mbdata_new)[[which(colnames(mbdata_new)=="dbname_d")]] <- "dbname"
    # dbname_e will be dropped because of the select= in the subset below.
  }
  if("COMMENT.EAWAG_UCHEM_ID" %in% colnames(mbdata_new))
    colnames(mbdata_new)[[which(colnames(mbdata_new)== "COMMENT.EAWAG_UCHEM_ID")]] <-
      "COMMENT.ID"
  
  # Clear from padding spaces and NAs
  mbdata_new <- as.data.frame(t(apply(mbdata_new, 1, function(r) 
    {
    # Substitute empty spaces by real NA values
    r[which(r == "")] <- NA
    # Trim spaces (in all non-NA fields)
    r[which(!is.na(r))] <- sub("^ *([^ ]+) *$", "\\1", r[which(!is.na(r))])
    return(r)
  })))
  # use only the columns present in mbdata_archive, no other columns added in excel
  mbdata_new <- mbdata_new[, colnames(mb@mbdata_archive)]
  # substitute the old entires with the ones from our files
  # then find the new (previously inexistent) entries, and rbind them to the table
  new_entries <- setdiff(mbdata_new$id, mb@mbdata_archive$id)
  old_entries <- intersect(mbdata_new$id, mb@mbdata_archive$id)
  for(entry in old_entries)
    mb@mbdata_archive[mb@mbdata_archive$id == entry,] <- mbdata_new[mbdata_new$id == entry,]
  mb@mbdata_archive <- rbind(mb@mbdata_archive, 
		  mbdata_new[mbdata_new$id==new_entries,])
  return(mb)

}


# Resets the mbdata_archive to an empty version.
#' @export
resetInfolists <- function(mb) 
{    
	mb@mbdata_archive <-
			structure(list(X = integer(0), id = integer(0), dbcas = character(0), 
							dbname = character(0), dataused = character(0), COMMENT.CONFIDENCE = character(0), 
							COMMENT.ID = integer(0), CH.NAME1 = character(0), 
							CH.NAME2 = character(0), CH.NAME3 = character(0), CH.COMPOUND_CLASS = character(0), 
							CH.FORMULA = character(0), CH.EXACT_MASS = numeric(0), CH.SMILES = character(0), 
							CH.IUPAC = character(0), CH.LINK.CAS = character(0), CH.LINK.CHEBI = integer(0), 
							CH.LINK.HMDB = character(0), CH.LINK.KEGG = character(0), CH.LINK.LIPIDMAPS = character(0), 
							CH.LINK.PUBCHEM = character(0), CH.LINK.INCHIKEY = character(0), 
							CH.LINK.CHEMSPIDER = integer(0)), .Names = c("X", "id", "dbcas", 
							"dbname", "dataused", "COMMENT.CONFIDENCE", "COMMENT.ID", 
							"CH.NAME1", "CH.NAME2", "CH.NAME3", "CH.COMPOUND_CLASS", "CH.FORMULA", 
							"CH.EXACT_MASS", "CH.SMILES", "CH.IUPAC", "CH.LINK.CAS", "CH.LINK.CHEBI", 
							"CH.LINK.HMDB", "CH.LINK.KEGG", "CH.LINK.LIPIDMAPS", "CH.LINK.PUBCHEM", 
							"CH.LINK.INCHIKEY", "CH.LINK.CHEMSPIDER"), row.names = integer(0), class = "data.frame")
	return(mb)
	
}

# The workflow function, i.e. (almost) the only thing you actually need to call. 
# See below for explanation of steps.
#' MassBank record creation workflow
#' 
#' Uses data generated by \code{\link{msmsWorkflow}} to create MassBank records.
#' 
#' See the vignette \code{vignette("RMassBank")} for detailed informations about the usage.
#' 
#' Steps:
#' 
#' Step 1: Find which compounds don't have annotation information yet. For these
#' 		 compounds, pull information from several databases (using gatherData).
#' 
#' Step 2: If new compounds were found, then export the infolist.csv and stop the workflow.
#' 		Otherwise, continue.
#' 
#' Step 3: Take the archive data (in table format) and reformat it to MassBank tree format.
#' 
#' Step 4: Compile the spectra. Using the skeletons from the archive data, create
#'   MassBank records per compound and fill them with peak data for each spectrum.
#'   Also, assign accession numbers based on scan mode and relative scan no.
#' 
#' Step 5: Convert the internal tree-like representation of the MassBank data into
#'  flat-text string arrays (basically, into text-file style, but still in memory)
#' 
#' Step 6: For all OK records, generate a corresponding molfile with the structure
#'   of the compound, based on the SMILES entry from the MassBank record. (This molfile
#'   is still in memory only, not yet a physical file)
#' 
#' Step 7: If necessary, generate the appropriate subdirectories, and actually write
#'   the files to disk.
#' 
#' Step 8: Create the list.tsv in the molfiles folder, which is required by MassBank
#'   to attribute substances to their corresponding structure molfiles. 
#' 
#' @param steps Which steps in the workflow to perform.
#' @param infolist_path A path where to store newly downloaded compound informations,
#' 			which should then be manually inspected.
#' @param mb The \code{mbWorkspace} to work in.
#' @param gatherData A variable denoting whether to retrieve information using several online databases \code{gatherData= "online"}
#' or to use the local babel installation \code{gatherData= "babel"}. Note that babel is used either way, if a directory is given 
#' in the settings
#' @return The processed \code{mbWorkspace}.
#' @seealso \code{\link{mbWorkspace-class}}
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @examples \dontrun{
#' 		mb <- newMbWorkspace(w) # w being a msmsWorkspace
#' 		mb <- loadInfolists(mb, "D:/myInfolistPath")
#' 		mb <- mbWorkflow(mb, steps=c(1:3), "newinfos.csv")
#' 		
#' }
#' @export
mbWorkflow <- function(mb, steps=c(1,2,3,4,5,6,7,8), infolist_path="./infolist.csv", gatherData = "online")
{
  # Step 1: Find which compounds don't have annotation information yet. For these
  # compounds, pull information from CTS (using gatherData).
  if(1 %in% steps)
  {
      mbdata_ids <- lapply(selectSpectra(mb@spectra, "found", "object"), function(spec) spec@id)
	  if(gatherData == "online"){
		message("mbWorkflow: Step 1. Gather info from several databases")
	  } 
	  if(gatherData == "babel"){
		message("mbWorkflow: Step 1. Gather info using babel")
	  }
      # Which IDs are not in mbdata_archive yet?
      new_ids <- setdiff(as.numeric(unlist(mbdata_ids)), mb@mbdata_archive$id)
      mb@mbdata <- lapply(new_ids, function(id) 
      {
        #print(id)
		if(gatherData == "online"){
			d <- gatherData(id)
		} 
		if(gatherData == "babel"){
			d <- gatherDataBabel(id)
		}
        #print(d$dataused)
		message(paste(id, ": ", d$dataused, sep=''))
        return(d)
      })
  }
  # Step 2: If new compounds were found, then export the infolist.csv and stop the workflow.
  # Otherwise, continue!
  if(2 %in% steps)
  {
	message("mbWorkflow: Step 2. Export infolist (if required)")
    if(length(mb@mbdata)>0)
    {
      mbdata_mat <- flatten(mb@mbdata)
      write.csv(as.data.frame(mbdata_mat),infolist_path, na="")
      message(paste("The file", infolist_path, "was generated with new compound information. Please check and edit the table, and add it to your infolist folder.")
            )
      return(mb)
    }
    else
      message("No new data added.")
  }
  # Step 3: Take the archive data (in table format) and reformat it to MassBank tree format.
  if(3 %in% steps)
  {
	message("mbWorkflow: Step 3. Data reformatting")
    mb@mbdata_relisted <- apply(mb@mbdata_archive, 1, readMbdata)
  }
  # Step 4: Compile the spectra! Using the skeletons from the archive data, create
  # MassBank records per compound and fill them with peak data for each spectrum.
  # Also, assign accession numbers based on scan mode and relative scan no.
  if(4 %in% steps)
  {
	  message("mbWorkflow: Step 4. Spectra compilation")
	  mb@compiled <- lapply(
			  selectSpectra(mb@spectra, "found", "object"),
			  function(r) {
				  message(paste("Compiling: ", r@name, sep=""))
				  mbdata <- mb@mbdata_relisted[[which(mb@mbdata_archive$id == as.numeric(r@id))]]
				  if(nrow(mb@additionalPeaks) > 0)
					  res <-compileRecord(r, mbdata, mb@aggregated, mb@additionalPeaks)
				  else
					  res <-compileRecord(r, mbdata, mb@aggregated, NULL)
				  return(res)
			  })
	  # check which compounds have useful spectra
	  mb@ok <- which(!is.na(mb@compiled) & !(lapply(mb@compiled, length)==0))
	  mb@problems <- which(is.na(mb@compiled))
	  mb@compiled_ok <- mb@compiled[mb@ok]
  }
  # Step 5: Convert the internal tree-like representation of the MassBank data into
  # flat-text string arrays (basically, into text-file style, but still in memory)
  if(5 %in% steps)
  {
	message("mbWorkflow: Step 5. Flattening records")
    mb@mbfiles <- lapply(mb@compiled_ok, function(c) lapply(c, toMassbank))
  }
  # Step 6: For all OK records, generate a corresponding molfile with the structure
  # of the compound, based on the SMILES entry from the MassBank record. (This molfile
  # is still in memory only, not yet a physical file)
  if(6 %in% steps)
  {
	message("mbWorkflow: Step 6. Generate molfiles")
    mb@molfile <- lapply(mb@compiled_ok, function(c) createMolfile(c[[1]][["CH$SMILES"]]))
  }
  # Step 7: If necessary, generate the appropriate subdirectories, and actually write
  # the files to disk.
  if(7 %in% steps)
  {
	message("mbWorkflow: Step 7. Generate subdirs and export")
    dir.create(paste(getOption("RMassBank")$annotations$entry_prefix, "moldata", sep='/'),recursive=TRUE)
    dir.create(paste(getOption("RMassBank")$annotations$entry_prefix, "recdata", sep='/'),recursive=TRUE)
    for(cnt in 1:length(mb@compiled_ok))
      exportMassbank(mb@compiled_ok[[cnt]], mb@mbfiles[[cnt]], mb@molfile[[cnt]])    
  }
  # Step 8: Create the list.tsv in the molfiles folder, which is required by MassBank
  # to attribute substances to their corresponding structure molfiles.
  if(8 %in% steps)
  {
	message("mbWorkflow: Step 8. Create list.tsv")
    makeMollist(mb@compiled_ok)
  }
  return(mb)
}


# Calls openbabel and converts the SMILES code string (or retrieves the SMILES code from
# the ID, and then calls openbabel) to create a molfile in text format.
# If fileName is given, the file is directly stored. Otherwise, it is returned as a 
# character array.
#' Create MOL file for a chemical structure
#' 
#' Creates a MOL file (in memory or on disk) for a compound specified by the
#' compound ID or by a SMILES code.
#' 
#' The function invokes OpenBabel (and therefore needs a correctly set
#' OpenBabel path in the RMassBank settings), using the SMILES code retrieved
#' with \code{findSmiles} or using the SMILES code directly. The current
#' implementation of the workflow uses the latter version, reading the SMILES
#' code directly from the MassBank record itself.
#' 
#' @usage createMolfile(id_or_smiles, fileName = FALSE)
#' @param id_or_smiles The compound ID or a SMILES code.
#' @param fileName If the filename is set, the file is written directly to disk
#' using the specified filename. Otherwise, it is returned as a text array.
#' @return A character array containing the MOL/SDF format file, ready to be
#' written to disk.
#' @author Michael Stravs
#' @seealso \code{\link{findSmiles}}
#' @references OpenBabel: \url{http://openbabel.org}
#' @examples
#' 
#' # Benzene:
#' \dontrun{
#' createMolfile("C1=CC=CC=C1")
#' }
#' 
#' @export
createMolfile <- function(id_or_smiles, fileName = FALSE)
{
	.checkMbSettings()
	babeldir <- getOption("RMassBank")$babeldir

	if(!is.numeric(id_or_smiles))
		smiles <- id_or_smiles
	else
		smiles <- findSmiles(id_or_smiles)
    # if no babeldir was set, get the result from cactus.
	if(is.na(babeldir))
	{
		res <- getCactus(smiles, "sdf")
		
		if(any(is.na(res))){
			res <- getPcSDF(smiles)
		}
		if(any(is.na(res))){
			stop("Pubchem and Cactus both seem to be down.")
		}
		if(is.character(fileName))
			writeLines(res, fileName)
	}
	# otherwise use the better-tested OpenBabel toolkit.
	else
	{
		if(!is.character(fileName))
			cmd <- paste(babeldir, "babel -ismi -osdf -d -b --gen2D", sep='')
		else
			cmd <- paste(babeldir, "babel -ismi -osdf ", fileName , " -d -b --gen2D", sep='')
		res <- system(cmd, intern=TRUE, input=smiles, ignore.stderr=TRUE)
		# If we wrote to a file, read it back as return value.
		if(is.character(fileName))
			res <- readLines(fileName)
	}
	return(res)
}



# Retrieve annotation data for a compound, from the internet service Pubchem
#' Retrieve supplemental annotation data from Pubchem
#' 
#' Retrieves annotation data for a compound from the internet service Pubchem 
#' based on the inchikey generated by babel or Cactus
#' 
#' The data retrieved is the Pubchem CID, a synonym from the Pubchem database,
#' the IUPAC name (using the preferred if available) and a Chebi link
#' 
#' @usage gatherPubChem(key)
#' @param key An Inchi-Key
#' @return Returns a list with 4 slots:
#' \code{PcID} The Pubchem CID
#' \code{Synonym} An arbitrary synonym for the compound name
#' \code{IUPAC} A IUPAC-name (preferred if available)
#' \code{Chebi} The identification number of the chebi database
#' @author Erik Mueller
#' @seealso \code{\link{mbWorkflow}}
#' @references Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' Chebi:
#' \url{http://www.ebi.ac.uk/chebi}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherPubChem("QEIXBXXKTUNWDK-UHFFFAOYSA-N")}
#' 
#' @export
gatherPubChem <- function(key){
	
	PubChemData <- list()
	
	##Trycatches are there because pubchem has connection issues 1 in 50 times.
	##Write NA into the respective fields if something goes wrong with the conenction or the data.
	
	##Retrieve Pubchem CID
	tryCatch(
		PubChemData$PcID <- getPcId(key),
		error=function(e){
		PubChemData$PcID <<- NA
	})
	
	##Retrieve a synonym to the name
	tryCatch(
		PubChemData$Synonym <- getPcSynonym(key),
		error=function(e){
		PubChemData$Synonym <<- NA
	})
	
	##Retrieve the IUPAC-name
	tryCatch(
		PubChemData$IUPAC <- getPcIUPAC(key),
		error=function(e){
		PubChemData$IUPAC <<- NA
	})
	
	##Retrieve the Chebi-ID
	tryCatch(
		PubChemData$Chebi <- getPcCHEBI(key),
		error=function(e){
		PubChemData$Chebi <<- NA
	})
	
	return(PubChemData)
}

# Retrieve annotation data for a compound, from the internet services Cactvs, Pubchem, Chemspider and CTS.
#' Retrieve annotation data
#' 
#' Retrieves annotation data for a compound from the internet services CTS, Pubchem, Chemspider and
#' Cactvs, based on the SMILES code and name of the compounds stored in the
#' compound list.
#' 
#' Composes the "upper part" of a MassBank record filled with chemical data
#' about the compound: name, exact mass, structure, CAS no., links to PubChem,
#' KEGG, ChemSpider.  The instrument type is also written into this block (even
#' if not strictly part of the chemical information). Additionally, index
#' fields are added at the start of the record, which will be removed later:
#' \code{id, dbcas, dbname} from the compound list, \code{dataused} to indicate
#' the used identifier for CTS search (\code{smiles} or \code{dbname}).
#' 
#' Additionally, the fields \code{ACCESSION} and \code{RECORD_TITLE} are
#' inserted empty and will be filled later on.
#' 
#' @usage gatherData(id)
#' @aliases gatherData
#' @param id The compound ID.
#' @return Returns a list of type \code{list(id= \var{compoundID}, ...,
#' 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc. %% ...
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}
#' @references Chemical Translation Service:
#' \url{http://uranus.fiehnlab.ucdavis.edu:8080/cts/homePage} 
#' cactus Chemical Identifier Resolver: 
#' \url{http://cactus.nci.nih.gov/chemical/structure}
#' MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' Chemspider InChI conversion:
#' \url{https://www.chemspider.com/InChI.asmx}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherData(131)}
#' 
#' @export
gatherData <- function(id)
{ 
	##Preamble: Is a babeldir supplied?
	##If yes, use it
	
	.checkMbSettings()
	usebabel=TRUE
	babeldir <- getOption("RMassBank")$babeldir
	
	if(is.na(babeldir)){
		usebabel=FALSE
	}
	
	
	##Get all useful information from the local "database" (from the CSV sheet)
	
	smiles <- findSmiles(id)
	mass <- findMass(smiles)
	dbcas <- findCAS(id)
	dbname <- findName(id)
	if(is.na(dbname)) dbname <- ""
	if(is.na(dbcas)) dbcas <- ""
	iupacName <- dbname
	synonym <- dbname
	formula <- findFormula(id)
	
	##Convert SMILES to InChI key via Cactvs or babel. CTS doesn't "interpret" the SMILES per se,
	##it just matches identical known SMILES, so we need to convert to a "searchable" and
	##standardized format beforehand. Other databases are able to interpret the smiles.
	
	if(usebabel){
		cmdinchikey <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchikey')
		inchikey_split <- system(cmdinchikey, intern=TRUE, input=smiles, ignore.stderr=TRUE)
	} else{
		inchikey <- getCactus(smiles, 'stdinchikey')
		if(!is.na(inchikey)){
			##Split the "InChiKey=" part off the key
			inchikey_split <- strsplit(inchikey, "=", fixed=TRUE)[[1]][[2]]
		} else{
		    inchikey_split <- getPcInchiKey(smiles)
		}
	}
	
	##Use Pubchem to retrieve information
	PcInfo <- gatherPubChem(inchikey_split)
	
	if(!is.null(PcInfo$Synonym) & !is.na(PcInfo$Synonym)){
		synonym <- PcInfo$Synonym
	}
	
	if(!is.null(PcInfo$IUPAC) & !is.na(PcInfo$IUPAC)){
		iupacName <- PcInfo$IUPAC
	}
	
	##Get Chemspider-ID
	csid <- getCSID(inchikey_split)
	
	if(is.na(csid)){
		##Get ChemSpider ID from Cactus if the Chemspider page is down
		csid <- getCactus(inchikey_split, 'chemspider_id')
	}
	
	##Use CTS to retrieve information
	CTSinfo <- getCtsRecord(inchikey_split)
		
	if((CTSinfo[1] == "Sorry, we couldn't find any matching results") || is.null(CTSinfo[1]))
	{
		CTSinfo <- NA
	}
	
	##List the names
	if(iupacName == ""){
		warning(paste0("Compound ID ",id,": no IUPAC name could be identified."))
	}

	if(toupper(dbname) == toupper(synonym)){
		synonym <- dbname
	}
	
	if(toupper(dbname) == toupper(iupacName)){
		iupacName <- dbname
	}
	
	if(toupper(synonym) == toupper(iupacName)){
		synonym <- iupacName
	}
	
	names <- as.list(unique(c(dbname, synonym, iupacName)))
	
	##If no name is found, it must be supplied in one way or another
	if(all(sapply(names, function(x) x == ""))){
		stop("RMassBank wasn't able to extract a usable name for this compound from any database. Please supply a name manually.")
	}
	
	# Start to fill the MassBank record.
	# The top 4 entries will not go into the final record; they are used to identify
	# the record and also to facilitate manual editing of the exported record table.
	mbdata <- list()
	mbdata[['id']] <- id
	mbdata[['dbcas']] <- dbcas
	mbdata[['dbname']] <- dbname
	mbdata[['dataused']] <- "smiles"
	mbdata[['ACCESSION']] <- ""
	mbdata[['RECORD_TITLE']] <- ""
	mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
	mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
	mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
	mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
	# Confidence annotation and internal ID annotation.
	# The ID of the compound will be written like:
	# COMMENT: EAWAG_UCHEM_ID 1234
	# if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
	mbdata[["COMMENT"]] <- list()
	mbdata[["COMMENT"]][["CONFIDENCE"]] <- getOption("RMassBank")$annotations$confidence_comment
	mbdata[["COMMENT"]][["ID"]] = id
    
	# here compound info starts
	mbdata[['CH$NAME']] <- names
	# Currently we use a fixed value for Compound Class, since there is no useful
	# convention of what should go there and what shouldn't, and the field is not used
	# in search queries.
	mbdata[['CH$COMPOUND_CLASS']] <- getOption("RMassBank")$annotations$compound_class
	mbdata[['CH$FORMULA']] <- formula
	mbdata[['CH$EXACT_MASS']] <- mass
	mbdata[['CH$SMILES']] <- smiles
	
	if(usebabel){
		cmdinchi <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchi')
		mbdata[['CH$IUPAC']] <- system(cmdinchi, intern=TRUE, input=smiles, ignore.stderr=TRUE)
	} else{
		mbdata[['CH$IUPAC']] <- getCactus(smiles, "stdinchi")
	}
	

	
	# Add all CH$LINK fields present in the compound datasets
	link <- list()
	# CAS
	if(!is.na(CTSinfo[1])){
		if("CAS" %in% CTS.externalIdTypes(CTSinfo))
		{
			# Prefer database CAS if it is also listed in the CTS results.
			# otherwise take the shortest one.
			cas <- CTS.externalIdSubset(CTSinfo,"CAS")
			if(dbcas %in% cas)
				link[["CAS"]] <- dbcas
			else
				link[["CAS"]] <- cas[[which.min(nchar(cas))]]
		} else{
			if(dbcas != ""){
				link[["CAS"]] <- dbcas
			}
		}
	} else{
		if(dbcas != ""){
			link[["CAS"]] <- dbcas
		}
	}
	
	
	# CHEBI
	if(is.na(PcInfo$Chebi[1])){
		if(!is.na(CTSinfo[1])){
			if("ChEBI" %in% CTS.externalIdTypes(CTSinfo))
			{
				# Cut off front "CHEBI:" if present
				chebi <- CTS.externalIdSubset(CTSinfo,"ChEBI")
				chebi <- chebi[[which.min(nchar(chebi))]]
				chebi <- strsplit(chebi,":")[[1]]
				link[["CHEBI"]] <- chebi[[length(chebi)]]
			}
		}
	} else{
		chebi <- PcInfo$Chebi
		chebi <- chebi[[which.min(nchar(chebi))]]
		chebi <- strsplit(chebi,":")[[1]]
		link[["CHEBI"]] <- chebi[[length(chebi)]]
	}
	# HMDB
	if(!is.na(CTSinfo[1])){
		if("Human Metabolome Database" %in% CTS.externalIdTypes(CTSinfo))
			link[["HMDB"]] <- CTS.externalIdSubset(CTSinfo,"HMDB")[[1]]
		# KEGG
		if("KEGG" %in% CTS.externalIdTypes(CTSinfo))
			link[["KEGG"]] <- CTS.externalIdSubset(CTSinfo,"KEGG")[[1]]
		# LipidMAPS
		if("LipidMAPS" %in% CTS.externalIdTypes(CTSinfo))
			link[["LIPIDMAPS"]] <- CTS.externalIdSubset(CTSinfo,"LipidMAPS")[[1]]
	}
	# PubChem CID
	if(is.na(PcInfo$PcID[1])){
		if(!is.na(CTSinfo[1])){
			if("PubChem CID" %in% CTS.externalIdTypes(CTSinfo))
			{
				pc <- CTS.externalIdSubset(CTSinfo,"PubChem CID")
				link[["PUBCHEM"]] <- paste0(min(pc))
			}
		}
	} else{
		link[["PUBCHEM"]] <- PcInfo$PcID[1]
	}
	
	
	if(!is.null(link[["PUBCHEM"]])){
		if(substr(link[["PUBCHEM"]],1,4) != "CID:"){
			link[["PUBCHEM"]] <- paste0("CID:", link[["PUBCHEM"]])
		}
	}
	
	link[["INCHIKEY"]] <- inchikey_split
	if(length(csid)>0) if(any(!is.na(csid))) link[["CHEMSPIDER"]] <- min(as.numeric(as.character(csid)))
	mbdata[['CH$LINK']] <- link
	
	mbdata[['AC$INSTRUMENT']] <- getOption("RMassBank")$annotations$instrument
	mbdata[['AC$INSTRUMENT_TYPE']] <- getOption("RMassBank")$annotations$instrument_type
	
	return(mbdata)  
}

# Retrieve annotation data for a compound, using only babel
#' Retrieve annotation data
#' 
#' Retrieves annotation data for a compound by using babel,
#' based on the SMILES code and name of the compounds stored in the
#' compound list.
#' 
#' Composes the "upper part" of a MassBank record filled with chemical data
#' about the compound: name, exact mass, structure, CAS no..  
#' The instrument type is also written into this block (even
#' if not strictly part of the chemical information). Additionally, index
#' fields are added at the start of the record, which will be removed later:
#' \code{id, dbcas, dbname} from the compound list.
#' 
#' Additionally, the fields \code{ACCESSION} and \code{RECORD_TITLE} are
#' inserted empty and will be filled later on.
#' 
#' This function is an alternative to gatherData, in case CTS is down or if information
#' on one or more of the compounds in the compound list are sparse
#'
#' @usage gatherDataBabel(id)
#' @param id The compound ID.
#' @return Returns a list of type \code{list(id= \var{compoundID}, ...,
#' 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc. %% ...
#' @author Michael Stravs, Erik Mueller
#' @seealso \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherDataBabel(131)}
#' 
#' @export
gatherDataBabel <- function(id){
		.checkMbSettings()
		babeldir <- getOption("RMassBank")$babeldir
		smiles <- findSmiles(id)
			
		
		# if no babeldir was set, throw an error that says that either CTS or babel have to be used
		if(is.na(babeldir))
		{
			stop("No babeldir supplied; It is currently not possible to convert the information without either babel or CTS")
		} else {
			###Babel conversion
			cmdinchikey <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchikey')
			inchikey <- system(cmdinchikey, intern=TRUE, input=smiles, ignore.stderr=TRUE)
			cmdinchi <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchi')
			inchi <- system(cmdinchi, intern=TRUE, input=smiles, ignore.stderr=TRUE)
			
			##Read from Compoundlist
			smiles <- findSmiles(id)
			mass <- findMass(smiles)
			dbcas <- findCAS(id)
			dbname <- findName(id)
			if(is.na(dbname)) dbname <- ""
			if(is.na(dbcas)) dbcas <- ""
			formula <- findFormula(id)
			
			##Create 
			mbdata <- list()
			mbdata[['id']] <- id
			mbdata[['dbcas']] <- dbcas
			mbdata[['dbname']] <- dbname
			mbdata[['dataused']] <- "smiles"
			mbdata[['ACCESSION']] <- ""
			mbdata[['RECORD_TITLE']] <- ""
			mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
			mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
			mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
			mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
			# Confidence annotation and internal ID annotation.
			# The ID of the compound will be written like:
			# COMMENT: EAWAG_UCHEM_ID 1234
			# if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
			mbdata[["COMMENT"]] <- list()
			mbdata[["COMMENT"]][["CONFIDENCE"]] <- getOption("RMassBank")$annotations$confidence_comment
			mbdata[["COMMENT"]][["ID"]] <- id

			# here compound info starts
			mbdata[['CH$NAME']] <- as.list(dbname)
			
			# Currently we use a fixed value for Compound Class, since there is no useful
			# convention of what should go there and what shouldn't, and the field is not used
			# in search queries.
			mbdata[['CH$COMPOUND_CLASS']] <- getOption("RMassBank")$annotations$compound_class
			mbdata[['CH$FORMULA']] <- formula
			mbdata[['CH$EXACT_MASS']] <- mass
			mbdata[['CH$SMILES']] <- smiles
			mbdata[['CH$IUPAC']] <- inchi
			
			link <- list()
			if(dbcas != "")
			link[["CAS"]] <- dbcas
			link[["INCHIKEY"]] <- inchikey
			mbdata[['CH$LINK']] <- link
			mbdata[['AC$INSTRUMENT']] <- getOption("RMassBank")$annotations$instrument
			mbdata[['AC$INSTRUMENT_TYPE']] <- getOption("RMassBank")$annotations$instrument_type
		}
		return(mbdata)
}


# Flatten the internal tree-like representation of MassBank data to a flat table.
# Note that this limits us, in that the fields should be constant over all records!
# Therefore, e.g. the fixed number of 3 names which may be filled.
# If anybody has a cooler solution, I'll be happy to hear from you :)
#
# Note: the records from gatherData have additional information which is discarded, like
# author, copyright etc. They will be re-filled automatically when reading the file.
#' Flatten, or re-read, MassBank header blocks
#' 
#' \code{flatten} converts a list of MassBank compound information sets (as
#' retrieved by \code{\link{gatherData}}) to a flat table, to be exported into
#' an \link[=loadInfolist]{infolist}. \code{readMbdata} reads a single record
#' from an infolist flat table back into a MassBank (half-)entry.
#' 
#' Neither the flattening system itself nor the implementation are particularly
#' fantastic, but since hand-checking of records is a necessary evil, there is
#' currently no alternative (short of coding a complete GUI for this and
#' working directly on the records.)
#' 
#' @aliases flatten readMbdata
#' @usage flatten(mbdata) 
#' 
#' readMbdata(row)
#' @param mbdata A list of MassBank compound information sets as returned from
#' \code{\link{gatherData}}.
#' @param row One row of MassBank compound information retrieved from an
#' infolist.
#' @return \code{flatten} returns a matrix (not a data frame) to be written to
#' CSV.
#' 
#' \code{readMbdata} returns a list of type \code{list(id= \var{compoundID},
#' ..., 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc.
#' @author Michael Stravs
#' @seealso \code{\link{gatherData}},\code{\link{loadInfolist}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples \dontrun{
#' 	# Collect some data to flatten
#' 	ids <- c(40,50,60,70)
#'  data <- lapply(ids, gatherData)
#'  # Flatten the data trees to a table
#'  flat.table <- flatten(data)
#'  # reimport the table into a tree
#'  data.reimported <- apply(flat.table, 1, readMbdata)
#' }
#' 
#' @export
#' 
flatten <- function(mbdata)
{
  .checkMbSettings()
  
  colList <- c(
              "id",
              "dbcas",
              "dbname",
              "dataused",
              "COMMENT.CONFIDENCE",
              # Note: The field name of the internal id field is replaced with the real name
              # at "compilation" time. Therefore, functions DOWNSTREAM from compileRecord() 
              # must use the full name including the info from options("RMassBank").
              "COMMENT.ID",
              "CH$NAME1",
              "CH$NAME2",
              "CH$NAME3",
              "CH$COMPOUND_CLASS",
              "CH$FORMULA",
              "CH$EXACT_MASS",
              "CH$SMILES",
              "CH$IUPAC",
              "CH$LINK.CAS",
              "CH$LINK.CHEBI",
              "CH$LINK.HMDB",
              "CH$LINK.KEGG",
              "CH$LINK.LIPIDMAPS",
              "CH$LINK.PUBCHEM",
              "CH$LINK.INCHIKEY",
              "CH$LINK.CHEMSPIDER")
  # make an empty data frame with the right length
  rows <- length(mbdata)
  cols <- length(colList)
  mbframe <- matrix(data=NA, nrow=rows, ncol=cols)
  colnames(mbframe) <- colList
  #browser()
  for(row in 1:rows)
  {
    # fill in all the data into the dataframe: all columns which 
    # a) exist in the target dataframe and b) exist in the (unlisted) MB record
    # are written into the dataframe.
    data <- unlist(mbdata[[row]])
	# bugfix for the case of only one name
	if(!("CH$NAME1" %in% names(data)))
		data[["CH$NAME1"]] <- data[["CH$NAME"]]
    datacols <- intersect(colList, names(data))
    mbframe[row,datacols] <- data[datacols]
  }
  return(mbframe)
  
}

# Read data from a flat-table MassBank record row and feed it into a
# MassBank tree-like record. Also, prime the ACCESSION and RECORD_TITLE fields in the
# correct position in the record.
#' @export
readMbdata <- function(row)
{
  .checkMbSettings()
  
  # Listify the table row. Lists are just cooler to work with :)
  row <- as.list(row)
  
  mbdata <- list()
  # Accession and title are added empty for now, to have them in the right place.
  # Constants are read from the options or generated.
  mbdata[['ACCESSION']] <- ""
  mbdata[['RECORD_TITLE']] <- ""
  mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
  mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
  mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
  mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
  mbdata[['PUBLICATION']] <- getOption("RMassBank")$annotations$publication
  
  # Read all determined fields from the file
  # This is not very flexible, as you can see...
    colList <- c(
              "COMMENT.CONFIDENCE",
              "COMMENT.ID",
              "CH$NAME1",
              "CH$NAME2",
              "CH$NAME3",
              "CH$COMPOUND_CLASS",
              "CH$FORMULA",
              "CH$EXACT_MASS",
              "CH$SMILES",
              "CH$IUPAC",
              "CH$LINK.CAS",
              "CH$LINK.CHEBI",
              "CH$LINK.HMDB",
              "CH$LINK.KEGG",
              "CH$LINK.LIPIDMAPS",
              "CH$LINK.PUBCHEM",
              "CH$LINK.INCHIKEY",
              "CH$LINK.CHEMSPIDER")
  mbdata[["COMMENT"]] = list()
  mbdata[["COMMENT"]][["CONFIDENCE"]] <- row[["COMMENT.CONFIDENCE"]]
  # Again, our ID field. 
  
  mbdata[["COMMENT"]][["ID"]]<-
            row[["COMMENT.ID"]]
  names = c(row[["CH.NAME1"]], row[["CH.NAME2"]], row[["CH.NAME3"]])
  names = names[which(!is.na(names))]
  
  names <- gsub("'", "`", names) 
  mbdata[["CH$NAME"]] = names
  mbdata[["CH$COMPOUND_CLASS"]] = row[["CH.COMPOUND_CLASS"]]
  mbdata[["CH$FORMULA"]] = row[["CH.FORMULA"]]
  mbdata[["CH$EXACT_MASS"]] = row[["CH.EXACT_MASS"]]
  mbdata[["CH$SMILES"]] = row[["CH.SMILES"]]
  mbdata[["CH$IUPAC"]] = row[["CH.IUPAC"]]
  # Add all links and then eliminate the NA values from the tree.
  link = list()
  link[["CAS"]] = row[["CH.LINK.CAS"]]
  link[["CHEBI"]] = row[["CH.LINK.CHEBI"]]
  link[["HMDB"]] = row[["CH.LINK.HMDB"]]
  link[["KEGG"]] = row[["CH.LINK.KEGG"]]
  link[["LIPIDMAPS"]] = row[["CH.LINK.LIPIDMAPS"]]
  link[["PUBCHEM"]] = row[["CH.LINK.PUBCHEM"]]
  link[["INCHIKEY"]] = row[["CH.LINK.INCHIKEY"]]
  link[["CHEMSPIDER"]] = row[["CH.LINK.CHEMSPIDER"]]
  link[which(is.na(link))] <- NULL
  mbdata[["CH$LINK"]] <- link
  # again, these constants are read from the options:
  mbdata[['AC$INSTRUMENT']] <- getOption("RMassBank")$annotations$instrument
  mbdata[['AC$INSTRUMENT_TYPE']] <- getOption("RMassBank")$annotations$instrument_type
  
  return(mbdata)
  
}

# For each compound, this function creates the "lower part" of the MassBank record, i.e.
# everything that comes after AC$INSTRUMENT_TYPE.
#' Compose data block of MassBank record
#' 
#' \code{gatherCompound} composes the data blocks (the "lower half") of all
#' MassBank records for a compound, using the annotation data in the RMassBank
#' options, spectrum info data from the \code{analyzedSpec}-type record and the
#' peaks from the reanalyzed, multiplicity-filtered peak table. It calls
#' \code{gatherSpectrum} for each child spectrum.
#' 
#' The returned data blocks are in format \code{list( "AC\$MASS_SPECTROMETRY" =
#' list('FRAGMENTATION_MODE' = 'CID', ...), ...)} etc.
#' 
#' @aliases gatherCompound gatherSpectrum
#' @usage gatherCompound(spec, aggregated, additionalPeaks = NULL)
#' 
#' 		gatherSpectrum(spec, msmsdata, ac_ms, ac_lc, aggregated, 
#'	 		additionalPeaks = NULL)
#' @param spec A \code{RmbSpectraSet} object, representing a compound with multiple spectra.
#' @param aggregated An aggregate peak table where the peaks are extracted from.
#' @param msmsdata A \code{RmbSpectrum2} object from the \code{spec} spectra set, representing a single spectrum to give a record.
#' @param ac_ms,ac_lc Information for the AC\$MASS_SPECTROMETRY and
#' AC\$CHROMATOGRAPHY fields in the MassBank record, created by
#' \code{gatherCompound} and then fed into \code{gatherSpectrum}.
#' @param additionalPeaks If present, a table with additional peaks to add into the spectra.
#' 		As loaded with \code{\link{addPeaks}}.
#' @return \code{gatherCompound} returns a list of tree-like MassBank data
#' blocks. \code{gatherSpectrum} returns one single MassBank data block or
#' \code{NA} if no useful peak is in the spectrum. 
#' @note Note that the global table \code{additionalPeaks} is also used as an
#' additional source of peaks.
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}, \code{\link{compileRecord}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples \dontrun{
#'      myspectrum <- w@@spectra[[1]]
#' 		massbankdata <- gatherCompound(myspectrum, w@@aggregated)
#' 		# Note: ac_lc and ac_ms are data blocks usually generated in gatherCompound and
#' 		# passed on from there. The call below gives a relatively useless result :)
#' 		ac_lc_dummy <- list()
#' 		ac_ms_dummy <- list() 
#' 		justOneSpectrum <- gatherSpectrum(myspectrum, myspectrum@@child[[2]],
#' 			ac_ms_dummy, ac_lc_dummy, w@@aggregated)
#' }
#' 
#' 
#' @export
gatherCompound <- function(spec, aggregated, additionalPeaks = NULL)
{
  # compound ID
  id <- spec@id
  # processing mode
  imode <- spec@mode
  
  # define positive or negative, based on processing mode.
  ion_modes <- list(
    "pH" = "POSITIVE",
    "pNa" = "POSITIVE",
    "mH" = "NEGATIVE",
    "mFA" = "NEGATIVE",
    "pM" = "POSITIVE",
    "mM" = "NEGATIVE",
	"pNH4" = "POSITIVE")
  mode <- ion_modes[[imode]]
  
  # for format 2.01
  ac_ms <- list();
  ac_ms[['MS_TYPE']] <- getOption("RMassBank")$annotations$ms_type
  ac_ms[['IONIZATION']] <- getOption("RMassBank")$annotations$ionization
  ac_ms[['ION_MODE']] <- mode
  
  # This list could be made customizable.
  ac_lc <- list();
  rt  <- spec@parent@rt / 60
  ac_lc[['COLUMN_NAME']] <- getOption("RMassBank")$annotations$lc_column
  ac_lc[['FLOW_GRADIENT']] <- getOption("RMassBank")$annotations$lc_gradient
  ac_lc[['FLOW_RATE']] <- getOption("RMassBank")$annotations$lc_flow
  ac_lc[['RETENTION_TIME']] <- sprintf("%.1f min", rt)  
  ac_lc[['SOLVENT A']] <- getOption("RMassBank")$annotations$lc_solvent_a
  ac_lc[['SOLVENT B']] <- getOption("RMassBank")$annotations$lc_solvent_b
  
  # Go through all child spectra, and fill our skeleton with scan data!
  # Pass them the AC_LC and AC_MS data, which are added at the right place
  # directly in there.
  allSpectra <- lapply(spec@children, function(m)
    gatherSpectrum(spec, m, ac_ms, ac_lc, aggregated, additionalPeaks))
  allSpectra <- allSpectra[which(!is.na(allSpectra))]
  return(allSpectra)
}



# Process one single MSMS child scan.
# spec: an object of "analyzedSpectrum" type (i.e. contains 
#       14x (or other number) msmsdata, info, mzrange,
#       compound ID, parent MS1, cpd id...)
# msmsdata: the msmsdata sub-object from the spec which is the child scan we want to process.
#       Contains childFilt, childBad, scan #, etc. Note that the peaks are actually not
#       taken from here! They were taken from msmsdata initially, but after introduction
#       of the refiltration and multiplicity filtering, this was changed. Now only the
#       scan information is actually taken from msmsdata.
# ac_ms, ac_lc: pre-filled info for the MassBank dataset (see above)
# refiltered: the refilteredRcSpecs dataset which contains our good peaks :)
#       Contains peaksOK, peaksReanOK, peaksFiltered, peaksFilteredReanalysis, 
#       peaksProblematic. Currently we use peaksOK and peaksReanOK to create the files.
#       (Also, the global additionalPeaks table is used.)
#' @export
gatherSpectrum <- function(spec, msmsdata, ac_ms, ac_lc, aggregated, additionalPeaks = NULL)
{
  # If the spectrum is not filled, return right now. All "NA" spectra will
  # not be treated further.
  if(msmsdata@ok == FALSE)
      return(NA)
  # get data
  scan <- msmsdata@acquisitionNum
  id <- spec@id
  # Further fill the ac_ms datasets, and add the ms$focused_ion with spectrum-specific data:
  precursor_types <- list(
    "pH" = "[M+H]+",
    "pNa" = "[M+Na]+",
    "mH" = "[M-H]-",
    "mFA" = "[M+HCOO-]-",
    "pM" = "[M]+",
    "mM" = "[M]-",
	"pNH4" = "[M+NH4]+")
  ac_ms[['FRAGMENTATION_MODE']] <- msmsdata@info$mode
  #ac_ms['PRECURSOR_TYPE'] <- precursor_types[spec$mode]
  ac_ms[['COLLISION_ENERGY']] <- msmsdata@info$ce
  ac_ms[['RESOLUTION']] <- msmsdata@info$res
  
  # Calculate exact precursor mass with Rcdk, and find the base peak from the parent
  # spectrum. (Yes, that's what belongs here, I think.)
  precursorMz <- findMz(spec@id, spec@mode)
  ms_fi <- list()
  ms_fi[['BASE_PEAK']] <- round(mz(spec@parent)[which.max(intensity(spec@parent))],4)
  ms_fi[['PRECURSOR_M/Z']] <- round(precursorMz$mzCenter,4)
  ms_fi[['PRECURSOR_TYPE']] <- precursor_types[spec@mode]
  
  # Select all peaks which belong to this spectrum (correct cpdID and scan no.)
  # from peaksOK
  # Note: Here and below it would be easy to customize the source of the peaks.
  # Originally the peaks came from msmsdata$childFilt, and the subset
  # was used where dppm == dppmBest (because childFilt still contains multiple formulas)
  # per peak.
	peaks <- aggregated[aggregated$filterOK,,drop=FALSE]
	peaks <- peaks[(peaks$cpdID == id) & (peaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
  
  # No peaks? Aha, bye
    if(nrow(peaks) == 0)
      return(NA)
  
  # If we don't include the reanalyzed peaks:
  if(!getOption("RMassBank")$use_rean_peaks)
	peaks <- peaks[is.na(peaks$matchedReanalysis),,drop=FALSE]
  # but if we include them:
  else
  {
	  # for info, the following data will be used in the default annotator:
	  # annotation <- annotation[,c("mzSpec","formula", "formulaCount", "mzCalc", "dppm")]
	  # and in the peaklist itself:
	  # c("mzSpec", "int", "intrel")
	  peaks[!is.na(peaks$matchedReanalysis),"formula"]  <- peaks[!is.na(peaks$matchedReanalysis),"reanalyzed.formula"]
	  peaks[!is.na(peaks$matchedReanalysis),"mzCalc"]  <- peaks[!is.na(peaks$matchedReanalysis),"reanalyzed.mzCalc"]
	  peaks[!is.na(peaks$matchedReanalysis),"dppm"]  <- peaks[!is.na(peaks$matchedReanalysis),"reanalyzed.dppm"]
	  peaks[!is.na(peaks$matchedReanalysis),"dbe"]  <- peaks[!is.na(peaks$matchedReanalysis),"reanalyzed.dbe"]
	  peaks[!is.na(peaks$matchedReanalysis),"formulaCount"]  <- peaks[!is.na(peaks$matchedReanalysis),"reanalyzed.formulaCount"]
  }

  # Calculate relative intensity and make a formatted m/z to use in the output
  # (mzSpec, for "spectrum")
  peaks$intrel <- floor(peaks$intensity / max(peaks$intensity) * 999)
  peaks$mzSpec <- round(peaks$mzFound, 4)
  # reorder peaks after addition of the reanalyzed ones
  peaks <- peaks[order(peaks$mzSpec),]
  
  # Also format the other values, which are used in the annotation
  peaks$dppm <- round(peaks$dppm, 2)
  peaks$mzCalc <- round(peaks$mzCalc, 4)
  peaks$intensity <- round(peaks$intensity, 1)
  # copy the peak table to the annotation table. (The peak table will then be extended
  # with peaks from the global "additional_peaks" table, which can be used to add peaks
  # to the spectra by hand.
  annotation <- peaks
  # Keep only peaks with relative intensity >= 1 o/oo, since the MassBank record
  # makes no sense otherwise. Also, keep only the columns needed in the output.
  peaks <- peaks[ peaks$intrel >= 1, c("mzSpec", "intensity", "intrel")]

  # Here add the additional peaks if there are any for this compound!
  # They are added without any annotation.
  if(!is.null(additionalPeaks))
  {
    # select the peaks from the corresponding spectrum which were marked with "OK=1" in the table.
    spec_add_peaks <- additionalPeaks[ 
			(additionalPeaks$OK == 1) & 
			(additionalPeaks$cpdID == spec@id) &
			(additionalPeaks$scan == msmsdata@acquisitionNum),
            c("mzFound", "intensity")]
    # If there are peaks to add:
    if(nrow(spec_add_peaks)>0)
    {
      # add the column for rel. int.
      spec_add_peaks$intrel <- 0
      # format m/z value
      spec_add_peaks$mzSpec <- round(spec_add_peaks$mzFound, 4)
      # bind tables together
      peaks <- rbind(peaks, spec_add_peaks[,c("mzSpec", "intensity", "intrel")])
      # recalculate rel.int.  and reorder list
      peaks$intrel <- floor(peaks$intensity / max(peaks$intensity) * 999)
      # Again, select the correct columns, and drop values with rel.int. <1 o/oo
      # NOTE: If the highest additional peak is > than the previous highest peak,
      # this can lead to the situation that a peak is in "annotation" but not in "peaks"!
      # See below.
      peaks <- peaks[ peaks$intrel >= 1, c("mzSpec", "intensity", "intrel")]
      # Reorder again.
      peaks <- peaks[order(peaks$mzSpec),]
    }
  }
  

  
  # add + or - to fragment formulas
  formula_tag <- list(
    "pH" = "+",
    "pNa" = "+",
    "mH" = "-",
    "mFA" = "-",
    "pM" = "+",
    "mM" = "-",
	"pNH4" = "+")
  type <- formula_tag[[spec@mode]]
  
  annotator <- getOption("RMassBank")$annotator
  if(is.null(annotator))
    annotator <- "annotator.default"
  
  
  
  # Here, the relative intensity is recalculated using the newly added additional
  # peaks from the peak list. Therefore, we throw superfluous peaks out again.
  # NOTE: It is a valid question whether or not we should kick peaks out at this stage.
  # The alternative would be to leave the survivors at 1 o/oo, but keep them in the spectrum.
  annotation$intrel <- floor(annotation$intensity / max(peaks$intensity) * 999)
  annotation <- annotation[annotation$intrel >= 1,]
  
  annotation <- do.call(annotator, list(annotation= annotation, type=type))
  
  
  # Name the columns correctly.
  colnames(peaks) <- c("m/z", "int.", "rel.int.")
  peaknum <- nrow(peaks)
  
  # Create the "lower part" of the record.  
  mbdata <- list()
  # Add the AC$MS, AC$LC info.
  if(getOption("RMassBank")$use_version == 2)
  {
    mbdata[["AC$MASS_SPECTROMETRY"]] <- ac_ms
    mbdata[["AC$CHROMATOGRAPHY"]] <- ac_lc
  }
  else
  {
      # Fix for MassBank data format 1, where ION_MODE must be renamed to MODE
      mbdata[["AC$ANALYTICAL_CONDITION"]] <- c(ac_ms, ac_lc)
      names(mbdata[["AC$ANALYTICAL_CONDITION"]])[[3]] <- "MODE"
  }
  # Add the MS$FOCUSED_ION info.
  mbdata[["MS$FOCUSED_ION"]] <- ms_fi

  ## The SPLASH is a hash value calculated across all peaks
  ## http://splash.fiehnlab.ucdavis.edu/
  ## Has to be temporarily added as "PK$SPLASH" in the "lower" part
  ## of the record, but will later be moved "up" when merging parts in compileRecord()  
  mbdata[["MS$RELATED_MS"]] <- list(SPLASH = getSplash(peaks[,c("m/z", "int.")]))
  
  # the data processing tag :)
  # Change by Tobias:
  # I suggest to add here the current version number of the clone due to better distinction between different makes of MB records
  # Could be automatised from DESCRIPTION file?
  if(getOption("RMassBank")$use_rean_peaks)
      processingComment <- list("REANALYZE" = "Peaks with additional N2/O included")
  else
      processingComment <- list()
  mbdata[["MS$DATA_PROCESSING"]] <- c(
    getOption("RMassBank")$annotations$ms_dataprocessing,
    processingComment,
    list("WHOLE" = paste("RMassBank", packageVersion("RMassBank")))
    )
  
  # Annotation:
  if(getOption("RMassBank")$add_annotation==TRUE)
    mbdata[["PK$ANNOTATION"]] <- annotation

  # Peak table
  mbdata[["PK$NUM_PEAK"]] <- peaknum
  mbdata[["PK$PEAK"]] <- peaks
  # These two entries will be thrown out later, but they are necessary to build the
  # record title and the accession number.
  mbdata[["RECORD_TITLE_CE"]] <- msmsdata@info$ces #formatted collision energy
  # Mode of relative scan calculation: by default it is calculated relative to the
  # parent scan. If a corresponding option is set, it will be calculated from the first
  # present child scan in the list.
  relativeScan <- "fromParent"
  if(!is.null(getOption("RMassBank")$recomputeRelativeScan))
	  if(getOption("RMassBank")$recomputeRelativeScan == "fromFirstChild")
		  relativeScan <- "fromFirstChild"
  if(relativeScan == "fromParent")
	  mbdata[["SUBSCAN"]] <- msmsdata@acquisitionNum - spec@parent@acquisitionNum #relative scan
  else if(relativeScan == "fromFirstChild")
  {
	  firstChild <- min(unlist(lapply(spec@children,function(d) d@acquisitionNum)))
	  mbdata[["SUBSCAN"]] <- msmsdata@acquisitionNum - firstChild + 1
  }
  return(mbdata)
}


# This compiles a MassBank record from the analyzedRcSpecs format (using the peaks from
# refilteredRcSpecs) together with the compound annotation data.
# Correspondingly:
# spec:       contains the analyzedRcSpec-format spectrum collection to be compiled
#             (i.e. a block of length(spectraList) child spectra)
# mbdata:     contains the corresponding MassBank "header" (the upper part of the record)
#             until INSTRUMENT TYPE.
# refiltered: the refilteredRcSpecs which contain our nice peaks.
#' Compile MassBank records
#' 
#' Takes a spectra block for a compound, as returned from
#' \code{\link{analyzeMsMs}}, and an aggregated cleaned peak table, together
#' with a MassBank information block, as stored in the infolists and loaded via
#' \code{\link{loadInfolist}}/\code{\link{readMbdata}} and processes them to a
#' MassBank record
#' 
#' \code{compileRecord} calls \code{\link{gatherCompound}} to create blocks of
#' spectrum data, and finally fills in the record title and accession number,
#' renames the "internal ID" comment field and removes dummy fields.
#' 
#' @usage compileRecord(spec, mbdata, aggregated, additionalPeaks = NULL)
#' @param spec A \code{RmbSpectraSet} for a compound, after analysis (\code{\link{analyzeMsMs}}).
#' Note that \bold{peaks are not read from this
#' object anymore}: Peaks come from the \code{aggregated} dataframe (and from
#' the global \code{additionalPeaks} dataframe; cf. \code{\link{addPeaks}} for
#' usage information.)
#' @param mbdata The information data block for the record header, as stored in
#' \code{mbdata_relisted} after loading an infolist.
#' @param aggregated An aggregated peak data table containing information about refiltered spectra etc.
#' @param additionalPeaks If present, a table with additional peaks to add into the spectra.
#' 		As loaded with \code{\link{addPeaks}}.
#' @return Returns a MassBank record in list format: e.g.
#' \code{list("ACCESSION" = "XX123456", "RECORD_TITLE" = "Cubane", ...,
#' "CH\$LINK" = list( "CAS" = "12-345-6", "CHEMSPIDER" = 1111, ...))}
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}, \code{\link{addPeaks}},
#' \code{\link{gatherCompound}}, \code{\link{toMassbank}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' 
#' #
#' \dontrun{myspec <- w@@spectra[[2]]}
#' # after having loaded an infolist:
#' \dontrun{mbdata <- mbdata_relisted[[which(mbdata_archive\$id == as.numeric(myspec\$id))]]}
#' \dontrun{compiled <- compileRecord(myspec, mbdata, w@@aggregated)}
#' 
#' @export
compileRecord <- function(spec, mbdata, aggregated, additionalPeaks = NULL)
{
  # gather the individual spectra data
  mblist <- gatherCompound(spec, aggregated, additionalPeaks)
  # this returns a n-member list of "lower parts" of spectra (one for each subscan).
  # (n being the number of child scans per parent scan.)
  # Now we put the two parts together.
  # (lapply on all n subscans, returns a list.)
  mblist_c <- lapply(mblist, function(l)
    {
      # This is the step which sticks together the upper and the lower part of the
      # record (the upper being compound-specific and the lower being scan-specific.)
      # Note that the accession number and record title (in the upper part) must of course
      # be filled in with scan-specific info.
      mbrecord <- c(mbdata, l)

      # Here is the right place to fix the name of the INTERNAL ID field.
      names(mbrecord[["COMMENT"]])[[which(names(mbrecord[["COMMENT"]]) == "ID")]] <-
        getOption("RMassBank")$annotations$internal_id_fieldname
	  # get mode parameter (for accession number generation) depending on version 
	  # of record definition
	  # Change by Tobias:
	  # I suggest to include fragmentation mode here for information
      if(getOption("RMassBank")$use_version == 2)
        mode <- mbrecord[["AC$MASS_SPECTROMETRY"]][["ION_MODE"]]
      else
        mode <- mbrecord[["AC$ANALYTICAL_CONDITION"]][["MODE"]]
	  # Generate the title and then delete the temprary RECORD_TITLE_CE field used before
	  mbrecord[["RECORD_TITLE"]] <- .parseTitleString(mbrecord)
      mbrecord[["RECORD_TITLE_CE"]] <- NULL
      # Calculate the accession number from the options.
      shift <- getOption("RMassBank")$accessionNumberShifts[[spec@mode]]
      mbrecord[["ACCESSION"]] <- sprintf("%s%04d%02d", getOption("RMassBank")$annotations$entry_prefix, as.numeric(spec@id), as.numeric(mbrecord[["SUBSCAN"]])+shift)
      # Clear the "SUBSCAN" field.
      mbrecord[["SUBSCAN"]] <- NULL
      # return the record.
      return(mbrecord)
  })
}



#' Generate peak annotation from peaklist
#' 
#' Generates the PK$ANNOTATION entry from the peaklist obtained. This function is
#' overridable by using the "annotator" option in the settings file.
#' 
#' @param annotation A peak list to be annotated. Contains columns:
#' \code{"cpdID","formula","mzFound" ,"scan","mzCalc","dppm",
#'      "dbe","mz","int","formulaCount","parentScan","fM_factor","dppmBest",
#'     "formulaMultiplicity","intrel","mzSpec"}
#' 
#' @param type The ion type to be added to annotated formulas ("+" or "-" usually)
#' 
#' @return The annotated peak table. Table \code{colnames()} will be used for the
#' 		titles (preferrably don't use spaces in the column titles; however no format is
#' 		strictly enforced by the MassBank data format.
#' 
#' @examples 
#' \dontrun{
#' annotation <- annotator.default(annotation)
#' }
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
annotator.default <- function(annotation, type)
{
  
  annotation$formula <- paste(annotation$formula, type, sep='')
  # Select the right columns and name them correctly for output.
  annotation <- annotation[,c("mzSpec","formula", "formulaCount", "mzCalc", "dppm")]
  colnames(annotation) <- c("m/z", "tentative_formula", "formula_count", "mass", "error(ppm)")
  return(annotation)
}

#' Parse record title
#' 
#' Parses a title for a single MassBank record using the title format
#' specified in the option titleFormat. Internally used, not exported.
#' 
#' If the option is not set, a standard title format is used (for record definition
#' version 1 or 2).
#' 
#' @usage .parseTitleString(mbrecord)
#' @param mbrecord A MassBank record in list format, as returned from
#' 	\code{\link{gatherSpectrum}}.
#' @return A string with the title.
#' @author Michael Stravs, Eawag
#' @seealso \code{\link{compileRecord}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' \dontrun{
#' 		# used in compileRecord()
#' 		title <- .parseTitleString(mbrecord)
#' }
#' 
#' 
#' 
.parseTitleString <- function(mbrecord)
{
	
	varlist <- getOption("RMassBank")$titleFormat
	
	# Set the standard title format.
	if(is.null(varlist))
	{
		if(getOption("RMassBank")$use_version == 2)
		{
			varlist <- c(
					"{CH$NAME}",
					"{AC$INSTRUMENT_TYPE}",
					"{AC$MASS_SPECTROMETRY: MS_TYPE}",
					"CE: {RECORD_TITLE_CE}",
					"R={AC$MASS_SPECTROMETRY: RESOLUTION}",
					"{MS$FOCUSED_ION: PRECURSOR_TYPE}"
			)
		}
		else
		{
			varlist <- c(
					"{CH$NAME}",
					"{AC$INSTRUMENT_TYPE}",
					"{AC$ANALYTICAL_CONDITION: MS_TYPE}",
					"CE: {RECORD_TITLE_CE}",
					"R={AC$ANALYTICAL_CONDITION: RESOLUTION}",
					"{MS$FOCUSED_ION: PRECURSOR_TYPE}"
			)
		}
	}
  
	
	# Extract a {XXX} argument from each title section.
	# check that every title has one and only one match
	args <- regexec("\\{(.*)\\}", varlist)
	arglist <- regmatches(varlist, args)
	if(any(unlist(lapply(arglist, length)) != 2))
		stop("Title format is incorrectly specified: a section with not exactly 1 parameters")
	
	parsedVars <- lapply(varlist, function(var)
			{
				# Extract the specified parameter inside the {}.
				# I.e. from a string like "R={BLA: BLUB}" return "BLA: BLUB"
				args <- regexec("\\{(.*)\\}", var)
				arg <- regmatches(var, args)[[1]][[2]]
				
				# Split the parameter by colon if necessary
				splitVar <- strsplit(arg, ": ")[[1]]
				# Read the parameter value from the record
				if(length(splitVar) == 2)
					replaceVar <- mbrecord[[splitVar[[1]]]][[splitVar[[2]]]]
				else if(length(splitVar) ==  1)
					replaceVar <- mbrecord[[splitVar]]
				else
					stop(paste(
									"Title format is incorrectly specified:", var)
					)
				# Fix problems: NULL returns
				if(is.null(replaceVar))
					replaceVar <- ""
				# Fix problems: Names will have >= 1 match. Take the first
				if(length(replaceVar) > 1)
					replaceVar <- replaceVar[[1]]
				# Substitute the parameter value into the string
				parsedVar <- sub("\\{(.*)\\}", replaceVar, var)	
				return(parsedVar)
			})
	title <- paste(parsedVars, collapse="; ")
	return(title)
}


# This converts the tree-like list (as obtained e.g. from compileRecord())
# into a plain text array, which can then be dumped to a file suitable for 
# MassBank upload.
#' Write MassBank record into character array
#' 
#' Writes a MassBank record in list format to a text array.
#' 
#' The function is a general conversion tool for the MassBank format; i.e. the
#' field names are not fixed. \code{mbdata} must be a named list, and the
#' entries can be as follows: \itemize{
#'  \item A single text line:
#' 
#' \code{'CH\$EXACT_MASS' = '329.1023'}
#' 
#'  is written as
#' 
#'  \code{CH\$EXACT_MASS: 329.1023} 
#' \item A character array:
#' 
#'  \code{'CH\$NAME' = c('2-Aminobenzimidazole', '1H-Benzimidazol-2-amine')} 
#' 
#' is written as
#' 
#' \code{CH\$NAME: 2-Aminobenzimidazole}
#' 
#' \code{CH\$NAME: 1H-Benzimidazol-2-amine}
#' 
#' \item A named list of strings: 
#' 
#' 	\code{'CH\$LINK' = list('CHEBI' = "27822", "KEGG" = "C10901")} 
#' 
#' is written as 
#' 
#' \code{CH\$LINK: CHEBI 27822}
#' 
#' \code{CH\$LINK: KEGG C10901} 
#' 
#' \item A data frame (e.g. the peak table) is written as specified in
#' the MassBank record format (Section 2.6.3): the column names are used as
#' headers for the first line, all data rows are printed space-separated. 
#' }
#' 
#' @usage toMassbank(mbdata)
#' @param mbdata A MassBank record in list format.
#' @return The result is a text array, which is ready to be written to the disk
#' as a file.
#' @note The function iterates over the list item names. \bold{This means that
#' duplicate entries in \code{mbdata} are (partially) discarded!} The correct
#' way to add them is by making a character array (as specified above): Instead
#' of \code{'CH\$NAME' = 'bla', 'CH\$NAME' = 'blub'} specify \code{'CH\$NAME' =
#' c('bla','blub')}.
#' @author Michael Stravs
#' @seealso \code{\link{compileRecord}}, \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' \dontrun{
#' # Read just the compound info skeleton from the Internet for some compound ID
#' id <- 35
#' mbdata <- gatherData(id)
#' #' # Export the mbdata blocks to line arrays 
#' # (there is no spectrum information, just the compound info...)
#' mbtext <- toMassbank(mbdata)
#' }
#' 
#' @export
toMassbank <- function (mbdata)
{
  # mbf is an array of lines and count is the line counter.
  # Very old-school, but it works. :)
  mbf <- character(0)
  count <- 1
  lapply(names(mbdata), function(entry)
    {
      # If entry is a char line, add it to the file.
      # If it is a named sublist, add each subentry with name
      # If it is an unnamed sublist, add each subentry without name
      # if it is a dataframe, write in PEAKS mode
    
      # Note: this is were I liked "lapply" a little too much. "for" would
      # be more idiomatic, and wouldn't need the <<- assignments.
      
      # Data frame: table mode. A header line and one space-separated line for
      # each data frame row.
      if(is.data.frame(mbdata[[entry]]))
      {
        mbf[[count]] <<- paste(entry,": " ,
                               paste(colnames(mbdata[[entry]]), collapse=" "),
                               sep='')
        count <<- count+1
        for(row in 1:nrow(mbdata[[entry]]))
        {
          mbf[[count]] <<- paste("  ", 
                                 paste(mbdata[[entry]][row,],collapse=" "), 
                                 sep="")
          count <<- count+1
        }
        #browser()
      }
      # List with named items: Write every entry like CH$LINK: CAS 12-345-678
      else if(is.list(mbdata[[entry]]) & !is.null(names(mbdata[[entry]])))
      {
        
        lapply(names(mbdata[[entry]]), function(subentry)
        {
          # todo
          mbf[[count]] <<- paste(entry,": ",subentry, " ", mbdata[[entry]][[subentry]], sep='')
          #print(mbf)
          count <<- count + 1
        })
      }
      # Array (or list) of unnamed items: Write every entry like CH$NAME: Paracetamol
      # (iterative entry without subindices)
      else if (length(mbdata[[entry]]) > 1 & is.null(names(mbdata[[entry]])))
      {
        lapply(mbdata[[entry]], function(subentry)
        {
          mbf[[count]] <<- paste(entry,": ",subentry, sep='')
          #print(mbf)
          count <<- count + 1
        })   
      }
      # Length is 1: just write the entry like PK$NUM_PEAKS: 131
      else
      {
        mbf[[count]] <<- paste(entry,": ",mbdata[[entry]], sep='')
        count <<- count + 1
      }
    }
    ) # End of lapply block (per child spectrum)
  # Add mandatory EOF marker
  mbf[[count]] <- "//"
  return(mbf)
}

# Exports compiled and massbanked spectra, with their associated molfiles, to physical files.
# "compiled" is still used here, because we need an accessible accession number.
# In the plain text arrays, the accession number is already "hidden".
# compiled: is ONE "compiled" entry, i.e. ONE compound with e.g. 14 spectra.
# files: is a return value from lapply(toMassbank), i.e. contains 14 plain-text arrays
#  (for a 14-spectra method)
# molfile: a molfile from createMolfile
#' Export internally stored MassBank data to files
#' 
#' Exports MassBank recfile data arrays and corresponding molfiles to physical
#' files on hard disk, for one compound.
#' 
#' The data from \code{compiled} is still used here, because it contains the
#' "visible" accession number. In the plain-text format contained in
#' \code{files}, the accession number is not "accessible" anymore since it's in
#' the file.
#' 
#' @usage exportMassbank(compiled, files, molfile)
#' @param compiled Is ONE "compiled" entry, i.e. ONE compound with e.g. 14
#' spectra, as returned from \code{\link{compileRecord}}.
#' @param files A n-membered array (usually a return value from
#' \code{lapply(\link{toMassbank})}), i.e. contains n plain-text arrays with
#' MassBank records.
#' @param molfile A molfile from \code{\link{createMolfile}}
#' @return No return value.
#' @note An improvement would be to write the accession numbers into
#' \code{names(compiled)} and later into \code{names(files)} so \code{compiled}
#' wouldn't be needed here anymore. (The compound ID would have to go into
#' \code{names(molfile)}, since it is also retrieved from \code{compiled}.)
#' @author Michael Stravs
#' @seealso \code{\link{createMolfile}}, \code{\link{compileRecord}},
#' \code{\link{toMassbank}}, \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' \dontrun{
#' 		compiled <- compileRecord(record, mbdata, refilteredRcSpecs)
#' 		mbfiles <- toMassbank(compiled)
#' 		molfile <- createMolfile(compiled[[1]][["CH$SMILES"]])
#' 		exportMassbank(compiled, mbfiles, molfile)
#' }
#' 
#' @export
exportMassbank <- function(compiled, files, molfile)
{
  molnames <- c()
  for(file in 1:length(compiled))
  {
    # Read the accession no. from the corresponding "compiled" entry
    filename <- compiled[[file]]["ACCESSION"]
    # use this accession no. as filename
    filename <- paste(filename, ".txt", sep="")
    write(files[[file]], 
          file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata",filename)
    )
  }
  # Use internal ID for naming the molfiles
  molname <- sprintf("%04d", as.numeric(
    compiled[[1]][["COMMENT"]][[getOption("RMassBank")$annotations$internal_id_fieldname]][[1]]))
  molname <- paste(molname, ".mol", sep="")
  write(molfile,
        file.path(getOption("RMassBank")$annotations$entry_prefix, "moldata",molname)
  )
}

# Makes a list.tsv with molfile -> massbank ch$name attribution.

#' Write list.tsv file
#' 
#' Makes a list.tsv file in the "moldata" folder.
#' 
#' Generates the list.tsv file which is needed by MassBank to connect records with
#' their respective molfiles. The first compound name is linked to a mol-file with
#' the compound ID (e.g. 2334.mol for ID 2334).
#' 
#' @param compiled A list of compiled spectra (in tree-format, as returned by \code{compileRecord}).
#' @return No return value.
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @examples \dontrun{
#' 		compiled <- compileRecord(record, mbdata, refilteredRcSpecs)
#' 		# a list.tsv for only one record:
#' 		clist <- list(compiled)
#' 		makeMollist(clist)
#' }
#' @export
makeMollist <- function(compiled)
{
  # For every "compiled" entry (here, compiled is not one "compiled" entry but the total
  # list of all compiled spectra), extract the uppermost CH$NAME and the ID (from the
  # first spectrum.) Make the ID into 0000 format.
  tsvlist <- t(sapply(compiled, function(entry)
    {
    name <- entry[[1]][["CH$NAME"]][[1]]
    id <- sprintf("%04d", as.numeric(entry[[1]][["COMMENT"]][[getOption("RMassBank")$annotations$internal_id_fieldname]][[1]]))
    molfilename <- paste(id,".mol",sep='')
    return(c(name,molfilename))
  }))
  # Write the file with the 
  write.table(tsvlist, 
              paste(getOption("RMassBank")$annotations$entry_prefix,"/moldata/list.tsv", sep=''),
              quote = FALSE,
              sep="\t",
              row.names=FALSE,
              col.names=FALSE
              )
}


# Load a dataframe or file into additional_peaks (or add additional points in there.)
# The columns cpdID, scan, mzFound, int, OK are mandatory. OK=1 means that the peaks
# will be added into the spectrum. mzFound and int will be taken for the table.
# No annotation will be written.
# Add peaks to the spectra by hand

#' Add additional peaks to spectra
#' 
#' Loads a table with additional peaks to add to the MassBank spectra. Required
#' columns are \code{cpdID, scan, int, mzFound, OK}.
#' 
#' All peaks with OK=1 will be included in the spectra.
#' 
#' @usage addPeaks(mb, filename_or_dataframe)
#' @param mb The \code{mbWorkspace} to load the peaks into.
#' @param filename_or_dataframe Filename of the csv file, or name of the R
#' dataframe containing the peaklist.
#' @return The \code{mbWorkspace} with loaded additional peaks.
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}
#' @examples
#' 
#' 	\dontrun{addPeaks("myrun_additionalPeaks.csv")}
#' 
#' @export 
addPeaks <- function(mb, filename_or_dataframe)
{
	
	errorvar <- 0
	currEnvir <- environment()
	d <- 1
	
	if(is.data.frame(filename_or_dataframe))
		df <- filename_or_dataframe
	else
	tryCatch(
		df <- read.csv(filename_or_dataframe),
		error=function(e){
		currEnvir$errorvar <- 1
	})
	# I change your heuristic fix to another heuristic fix, because I will have to test for a column name change...
	
	if(!errorvar){
	
		if(ncol(df) < 2)
			df <- read.csv(filename_or_dataframe, sep=";")
		# here: the column int was renamed to intensity, and we need to be able to read old files. sorry.
		if(!("intensity" %in% colnames(df)) & ("int" %in% colnames(df)))
			df$intensity <- df$int
		
		cols <- c("cpdID", "scan", "mzFound", "intensity", "OK")
		n <- colnames(df)
		# Check if comma-separated or semicolon-separated
		d <- setdiff(cols, n)
		if(length(d)>0){
			stop("Some columns are missing in the additional peak list. Needs at least cpdID, scan, mzFound, intensity, OK.")
		}
	}
	
	culled_df <- df[,c("cpdID", "scan", "mzFound", "intensity", "OK")]
	
	
	if(nrow(mb@additionalPeaks) == 0)
		mb@additionalPeaks <- culled_df
	else
		mb@additionalPeaks <- rbind(mb@additionalPeaks, culled_df)
	return(mb)
}
