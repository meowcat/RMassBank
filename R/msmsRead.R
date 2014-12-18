#' 
#' Extracts and processes spectra from a specified file list, according to 
#' loaded options and given parameters.
#' 
#' The filenames of the raw LC-MS runs are read from the array \code{files} 
#' in the global enviroment.
#' See the vignette \code{vignette("RMassBank")} for further details about the
#' workflow.
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param filetable The path to a .csv-file that contains the columns "Files" and "ID" supplying
#' 			the relationships between files and compound IDs. Either this or the parameter "files" need
#'			to be specified. 
#' @param files A vector or list containing the filenames of the files that are to be read as spectra. 
#'				For the IDs to be inferred from the filenames alone, there need to be exactly 2 underscores.
#' @param cpdids A vector or list containing the compound IDs of the files that are to be read as spectra.
#'				The ordering of this and \code{files} implicitly assigns each ID to the corresponding file.
#'				If this is supplied, then the IDs implicitly named in the filenames are ignored.
#' @param readMethod Several methods are available to get peak lists from the files.
#'        Currently supported are "mzR", "xcms", "MassBank" and "peaklist".
#'        The first two read MS/MS raw data, and differ in the strategy 
#'        used to extract peaks. MassBank will read existing records, 
#'        so that e.g. a recalibration can be performed, and "peaklist" 
#'        just requires a CSV with two columns and the column header "mz", "int".
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param confirmMode Defaults to false (use most intense precursor). Value 1 uses
#' 			the 2nd-most intense precursor for a chosen ion (and its data-dependent scans)
#' 			, etc.
#' @param useRtLimit Whether to enforce the given retention time window.
#' @param Args A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param MSe A boolean value that determines whether the spectra were recorded using MSe or not
#' @param plots A boolean value that determines whether the pseudospectra in XCMS should be plotted
#' @return The \code{msmsWorkspace} with msms-spectra read.
#' @seealso \code{\link{msmsWorkspace-class}}, \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
msmsRead <- function(w, filetable = NULL, files = NULL, cpdids = NULL, 
					readMethod, mode, confirmMode = FALSE, useRtLimit = TRUE, 
					Args = NULL, settings = getOption("RMassBank"), progressbar = "progressBarHook", MSe = FALSE, plots = FALSE){
	
	##Read the files and cpdids according to the definition
	##All cases are silently accepted, as long as they can be handled according to one definition
	if(!any(mode %in% c("pH","pNa","pM","mH","mFA","mM",""))) stop(paste("The ionization mode", mode, "is unknown."))
	
	if(is.null(filetable)){
		##If no filetable is supplied, filenames must be named explicitly
		if(is.null(files))
			stop("Please supply the files")
		
		##Assign the filenames to the workspace
		w@files <- unlist(files)
		
		##If no filetable is supplied, cpdids must be delivered explicitly or implicitly within the filenames
		if(is.null(cpdids)){
			splitfn <- strsplit(files,"_")
			splitsfn <- sapply(splitfn, function(x) x[length(x)-1])
			if(suppressWarnings(any(is.na(as.numeric(splitsfn)[1]))))
				stop("Please supply the cpdids corresponding to the files in the filetable or the filenames")
			cpdids <- splitsfn
		}
	} else{
		##If a filetable is supplied read it
		tab <- read.csv(filetable, stringsAsFactors = FALSE)
		w@files <- tab[,"Files"]
		cpdids <- tab[,"ID"]
	}
	
	##If there's more cpdids than filenames or the other way around, then abort
	if(length(w@files) != length(cpdids)){
		stop("There are a different number of cpdids than files")
	}
	if(!(readMethod %in% c("mzR","peaklist","xcms","simple"))){
		stop("The supplied method does not exist")
	}
	if(!all(file.exists(w@files))){
		stop("The supplied files ", paste(w@files[!file.exists(w@files)]), " don't exist")
	}

	##This should work
	if(readMethod == "simple"){
		##Edit options
		opt <- getOption("RMassBank")
		opt$recalibrator$MS1 <- "recalibrate.identity"
		opt$recalibrator$MS2 <- "recalibrate.identity"
		options(RMassBank=opt)
		##Edit analyzemethod
		analyzeMethod <- "intensity"
	}
	
	if(readMethod == "mzR"){
		##Progressbar
		nLen <- length(w@files)
		nProg <- 0
		pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
		
		count <- 1
		w@specs <-  lapply(w@files, function(fileName){
			spec <- findMsMsHR(fileName, cpdids[count], mode, confirmMode, useRtLimit,
		 		ppmFine = settings$findMsMsRawSettings$ppmFine,
		 		mzCoarse = settings$findMsMsRawSettings$mzCoarse,
		 		fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
		 		rtMargin = settings$rtMargin,
		 		deprofile = settings$deprofile)
			
			## Progress:
			nProg <<- nProg + 1
			pb <- do.call(progressbar, list(object=pb, value= nProg))
			
			##Counting the index of cpdids
			count <<- count + 1
			return(spec)
		})
		names(w@specs) <- basename(as.character(w@files))
		return(w)
	}
	
	##Peaklist-readmethod 
	if((readMethod == "peaklist") || (readMethod=="simple")){
		w <- createSpecsFromPeaklists(w, cpdids, filenames=w@files, mode=mode)
		uIDs <- unique(cpdids)
		files <- list()
		
		for(i in 1:length(uIDs)){
			indices <- sapply(cpdids,function(a){return(uIDs[i] %in% a)})
			files[[i]] <- w@files[indices]
		}
		
		w@files <- sapply(files,function(file){return(file[1])})
		specnames <- basename(as.character(w@files))
		if(length(unique(specnames)) == length(specnames)){
			names(w@specs) <- basename(as.character(w@files))
		} else {
			for(i in 1:length(specnames)){
				specnames[i] <- paste(i,"_",specnames[i],sep="")
			}
		}
		names(w@specs) <- specnames
		message("Peaks read")
		return(w)
	}
	
	##xcms-readmethod 
	if(readMethod == "xcms"){
		require(xcms)
		require(CAMERA)
		ufiles <- unique(w@files)
		uIDs <- unique(cpdids)
		##Routine for the case of multiple cpdIDs per file and multiple files per cpdID
		dummySpecs <- list()
		w@specs <- list()
		
		nLen <- length(ufiles)
		nProg <- 0
		pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
		i <- 1
		dummyapply <- lapply(ufiles, function(currfile){
			
			dummySpecs[[i]] <<- newMsmsWorkspace()
			dummySpecs[[i]]@specs <<- list()
			FileIDs <<- cpdids[which(w@files == currfile)]
			dummylinebreak <- capture.output(metaSpec <<- findMsMsHRperxcms.direct(currfile, FileIDs, mode=mode, findPeaksArgs=Args, plots, MSe = MSe))
			
			for(j in 1:length(FileIDs)){
				dummySpecs[[i]]@specs[[length(dummySpecs[[i]]@specs)+1]] <<- metaSpec[[j]]
			}
			
			##Progress
			nProg <<- nProg + 1
			pb <- do.call(progressbar, list(object=pb, value= nProg))
			i <<- i+1
			})
			
			if(length(dummySpecs) > 1){
				for(j in 2:length(dummySpecs)){
					dummySpecs[[1]] <- c.msmsWSspecs(dummySpecs[[1]],dummySpecs[[j]])
				}
			}
			
			##You need as many names as there were different IDs
			##And the Names and IDs have to go together in some way
			##Find out Names that make sense: (cpdID with Name of File that uses cpdID)
			FNames <- vector()
			for(i in uIDs){
				nindex <- min(which(i == cpdids))
				FNames <- c(FNames,paste(w@files[nindex],"_",cpdids[nindex],sep=""))
			}
			
			w@specs <- dummySpecs[[1]]@specs
			names(w@specs) <- basename(as.character(FNames))
			w@files <- basename(as.character(FNames))
			return(w)
	}
}

msmsRead.RAW <- function(w, xRAW = NULL, cpdids = NULL, mode, confirmMode = FALSE, useRtLimit = TRUE, 
						findPeaksArgs = NULL, settings = getOption("RMassBank"), progressbar = "progressBarHook", plots = FALSE){
	
	##xRAW will be coerced into a list of length 1 if it is an xcmsRaw-object
	if(class(xRAW) == "xcmsRaw"){
		xRAW <- list(xRAW)
	}
	
	##Error messages
	if((class(xRAW) != "list") || any(sapply(xRAW, function(x) class(x) != "xcmsRaw"))){
		stop("No list of xcmsRaw-objects supplied")
	}
	
	if(is.null(cpdids)){
		stop("No cpdids supplied")
	}
		
	#msnExist <- which(sapply(xRAW,function(x) length(x@msnPrecursorScan) != 0))
	#print(length(msnExist))
	#print(length(xRAW))
	
	#if(length(msnExist) != length(xRAW)){
	#	stop(paste("No msn data in list elements", setdiff(1:length(xRAW),msnExist)))
	#}
	
	require(CAMERA)
	
	parentMass <- findMz(cpdids[1], mode=mode)$mzCenter
	if(is.na(parentMass)){
		stop(paste("There was no matching entry to the supplied cpdID", cpdids[1] ,"\n Please check the cpdIDs and the compoundlist."))
	}
		
	RT <- findRt(cpdids[1])$RT * 60
	mzabs <- 0.1
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(peaks(xa@xcmsSet)[x, "rt"])})
	}
	
	suppressWarnings(setReplicate <- xcmsSet(files=xRAW[[1]]@filepath, method="MS1"))
	xsmsms <- as.list(replicate(length(xRAW),setReplicate))
	candidates <- list()
	anmsms <- list()
	psp <- list()
	spectra <- list()
	whichmissing <- vector()
	metaspec <- list()
	for(i in 1:length(xRAW)){
		devnull <- suppressWarnings(capture.output(peaks(xsmsms[[i]]) <- do.call(findPeaks,c(findPeaksArgs, object = xRAW[[i]]))))
		
		if (nrow(peaks(xsmsms[[i]])) == 0) { ##If there are no peaks
			spectra[[i]] <- matrix(0,2,7)
			next
		} else{	
			## Get pspec 
			pl <- peaks(xsmsms[[i]])[,c("mz", "rt"), drop=FALSE]

			## Best: find precursor peak
			candidates[[i]] <- which( pl[,"mz", drop=FALSE] < parentMass + mzabs & pl[,"mz", drop=FALSE] > parentMass - mzabs
							& pl[,"rt", drop=FALSE] < RT * 1.1 & pl[,"rt", drop=FALSE] > RT * 0.9 )
			devnull <- capture.output(anmsms[[i]] <- xsAnnotate(xsmsms[[i]]))
			devnull <- capture.output(anmsms[[i]] <- groupFWHM(anmsms[[i]]))

				if(length(candidates[[i]]) > 0){
				closestCandidate <- which.min (abs( RT - pl[candidates[[i]], "rt", drop=FALSE]))
				psp[[i]] <- which(sapply(anmsms[[i]]@pspectra, function(x) {candidates[[i]][closestCandidate] %in% x}))
				} else{psp[[i]] <- which.min( abs(getRT(anmsms[[i]]) - RT) )}
				## Now find the pspec for compound       

				## 2nd best: Spectrum closest to MS1
				##psp <- which.min( abs(getRT(anmsms) - actualRT))

				## 3rd Best: find pspec closest to RT from spreadsheet
				##psp <- which.min( abs(getRT(anmsms) - RT) )
				if((plots == TRUE) && (length(psp[[i]]) > 0)){
					plotPsSpectrum(anmsms[[i]], psp[[i]], log=TRUE,  mzrange=c(0, findMz(cpdids)[[3]]), maxlabel=10)
				}
				if(length(psp[[i]]) != 0){
					spectra[[i]] <- getpspectra(anmsms[[i]], psp[[i]])
				} else {whichmissing <- c(whichmissing,i)}
			}
		}
		if(length(spectra) != 0){
			for(i in whichmissing){
				spectra[[i]] <- matrix(0,2,7)
			}
		}
	if(length(w@specs) != 0){
		w@specs <- c(w@specs,list(toRMB(spectra,cpdids[1],mode)))
	} else {
		w@specs[[1]] <- toRMB(spectra,cpdids[1],mode)
	}
	w@files <- c(w@files,xRAW[[1]]@filepath)
	names(w@specs)[length(w@specs)] <- basename(xRAW[[1]]@filepath)
	return(w)
}