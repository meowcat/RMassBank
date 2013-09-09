#' @import mzR
#' @importClassesFrom mzR
#' @importMethodsFrom mzR
# library(mzR)
NULL

#' Extract MS/MS spectra for specified precursor
#' 
#' Extracts MS/MS spectra from LC-MS raw data for a specified precursor, specified
#' either via the RMassBank compound list (see \code{\link{loadList}}) or via a mass.
#' 
#' Different versions of the function get the data from different sources. Note that 
#' 		findMsMsHR and findMsMsHR.direct differ mainly in that findMsMsHR opens a file
#' 		whereas findMsMs.direct uses an open file handle - both are intended to be used
#' 		in a full process which involves compound lists etc. In contrast, findMsMsHR.mass
#' 		is a low-level function which uses the mass directly for lookup and is intended for
#' 		use as a standalone function in unrelated applications.
#' 
#' @usage findMsMsHR(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
#' 		ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
#' 		mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
#' 		fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
#' 		rtMargin = getOption("RMassBank")$rtMargin,
#' 		deprofile = getOption("RMassBank")$deprofile)
#' 		
#' 		findMsMsHR.mass(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
#' 		headerCache = NA, fillPrecursorScan = FALSE,
#' 		deprofile = getOption("RMassBank")$deprofile)
#'
#' findMsMsHR.direct(msRaw, cpdID, mode = "pH", confirmMode = 0, useRtLimit = TRUE, 
#'			ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
#'			mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
#'			fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
#'			rtMargin = getOption("RMassBank")$rtMargin,
#'			deprofile = getOption("RMassBank")$deprofile)
#' 
#' @aliases findMsMsHR.mass findMsMsHR.direct findMsMsHR
#' @param fileName The file to open and search the MS2 spectrum in.
#' @param msRaw The opened raw file (mzR file handle) to search the MS2 spectrum in.
#' @param cpdID The compound ID in the compound list (see \code{\link{loadList}})
#' 			to use for formula lookup.
#' @param mz The mass to use for spectrum search.
#' @param ppmFine The limit in ppm to use for fine limit (see below) calculation.
#' @param mzCoarse The coarse limit to use for locating potential MS2 scans:
#'			this tolerance is used when finding scans with a suitable precursor
#' 			ion value.  
#' @param limit.fine The fine limit to use for locating MS2 scans: this tolerance
#' 			is used when locating an appropriate analyte peak in the MS1 precursor
#' 			spectrum. 
#' @param limit.coarse Parameter in \code{findMsMsHR.mass} corresponding to \code{mzCoarse}.
#' 			(The parameters are distinct to clearly conceptually distinguish findMsMsHR.mass
#' 			(a standalone useful function) from the cpdID based functions (workflow functions).)
#' @param mode The processing mode (determines which ion/adduct is searched):
#' 			\code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-). 
#' @param confirmMode Whether to use the highest-intensity precursor (=0), second-
#' 			highest (=1), third-highest (=2)...
#' @param useRtLimit Whether to respect retention time limits from the compound list.
#' @param rtLimits \code{c(min, max)}: Minimum and maximum retention time to use
#' 			when locating the MS2 scans. 
#' @param headerCache If present, the complete \code{mzR::header(msRaw)}. Passing
#' 			this value is useful if spectra for multiple compounds should be 
#' 			extracted from the same mzML file, since it avoids getting the data
#' 			freshly from \code{msRaw} for every compound.
#' @param maxCount The maximal number of spectra groups to return. One spectra group
#' 			consists of all data-dependent scans from the same precursor whose precursor
#' 			mass matches the specified search mass.
#' @param fillPrecursorScan If \code{TRUE}, the precursor scan will be filled from MS1 data.
#' 			To be used for data where the precursor scan is not stored in the raw data.
#' @param rtMargin	The retention time tolerance to use.
#' @param deprofile	Whether deprofiling should take place, and what method should be
#' 			used (cf. \code{\link{deprofile}}) 
#' @return	For \code{findMsMsHR} and \code{findMsMsHR.direct}: A "spectrum set", a list with items:
#' 			\item{foundOK}{\code{TRUE} if a spectrum was found, \code{FALSE} otherwise.
#' 				Note: if \code{FALSE}, all other values can be missing!}
#' 			\item{parentScan}{The scan number of the precursor scan.}
#' 			\item{parentHeader}{The header row of the parent scan, as returned by 
#' 				\code{mzR::header}.}
#' 			\item{childScans}{The scan numbers of the data-dependent MS2 scans.}
#' 			\item{childHeaders}{The header rows of the MS2 scan, as returned by
#' 				\code{mzR::header}.}
#' 			\item{parentPeak}{The MS1 precursor spectrum as a 2-column matrix}
#' 			\item{peaks}{A list of  2-column \code{mz, int} matrices of the MS2 scans.}
#' 			For \code{findMsMsHR.mass}: a list of "spectrum sets" as defined above, sorted
#' 			by decreasing precursor intensity.
#' 
#' @examples \dontrun{
#' 			loadList("mycompoundlist.csv")
#' 			# if Atrazine has compound ID 1:
#' 			msms_atrazine <- findMsMsHR("Atrazine_0001_pos.mzML", 1, "pH")
#' 			# Or alternatively:
#' 			msRaw <- openMSfile("Atrazine_0001_pos.mzML")
#' 			msms_atrazine <- findMsMsHR.direct(msRaw, 1, "pH")
#' 			# Or directly by mass (this will return a list of spectra sets):
#' 			mz <- findMz(1)$mzCenter
#' 			msms_atrazine_all <- findMsMsHR.mass(msRaw, mz, 1, ppm(msRaw, 10, p=TRUE))
#' 			msms_atrazine <- msms_atrazine_all[[1]]
#' }
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findEIC
#' @export
findMsMsHR <- function(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
		ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
		mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
		fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
		rtMargin = getOption("RMassBank")$rtMargin,
		deprofile = getOption("RMassBank")$deprofile)
{
	
	# access data directly for finding the MS/MS data. This is done using
	# mzR.
	msRaw <- openMSfile(fileName)
	ret <- findMsMsHR.direct(msRaw, cpdID, mode, confirmMode, useRtLimit, ppmFine, mzCoarse, fillPrecursorScan,
				rtMargin, deprofile)
	mzR::close(msRaw)
	return(ret)
}

#' @export
findMsMsHR.mass <- function(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
		headerCache = NA, fillPrecursorScan = FALSE,
		deprofile = getOption("RMassBank")$deprofile)
{
	eic <- findEIC(msRaw, mz, limit.fine, rtLimits)
	#	if(!is.na(rtLimits))
	#	{  
	#		eic <- subset(eic, rt >= rtLimits[[1]] & rt <= rtLimits[[2]])
	#	}
	if(!all(is.na(headerCache)))
		headerData <- headerCache
	else
		headerData <- as.data.frame(header(msRaw))
	
	if(fillPrecursorScan == TRUE)
	{
		# reset the precursor scan number. first set to NA, then
		# carry forward the precursor scan number from the last parent scan
		headerData$precursorScanNum <- NA
		headerData[which(headerData$msLevel == 1),"precursorScanNum"] <-
				headerData[which(headerData$msLevel == 1),"acquisitionNum"]
		headerData[,"precursorScanNum"] <- .locf(headerData[,"precursorScanNum"])
		# Clear the actual MS1 precursor scan number again
		headerData[which(headerData$msLevel == 1),"precursorScanNum"] <- 0
	}
	
	# Find MS2 spectra with precursors which are in the allowed 
	# scan filter (coarse limit) range
	findValidPrecursors <- headerData[
			(headerData$precursorMZ > mz - limit.coarse) &
					(headerData$precursorMZ < mz + limit.coarse),]
	# Find the precursors for the found spectra
	validPrecursors <- unique(findValidPrecursors$precursorScanNum)
	# check whether the precursors are real: must be within fine limits!
	# previously even "bad" precursors were taken. e.g. 1-benzylpiperazine
	which_OK <- lapply(validPrecursors, function(pscan)
			{
				pplist <- as.data.frame(
						mzR::peaks(msRaw, which(headerData$acquisitionNum == pscan)))
				colnames(pplist) <- c("mz","int")
				pplist <- pplist[(pplist$mz >= mz -limit.fine)
								& (pplist$mz <= mz + limit.fine),]
				if(nrow(pplist) > 0)
					return(TRUE)
				return(FALSE)
			})
	validPrecursors <- validPrecursors[which(which_OK==TRUE)]
	# Crop the "EIC" to the valid precursor scans
	eic <- eic[eic$scan %in% validPrecursors,]
	# Order by intensity, descending
	eic <- eic[order(eic$intensity, decreasing=TRUE),]
	if(nrow(eic) == 0)
		return(list(list(foundOK = FALSE)))
	if(!is.na(maxCount))
	{
		spectraCount <- min(maxCount, nrow(eic))
		eic <- eic[1:spectraCount,]
	}
	# Construct all spectra groups in decreasing intensity order
	spectra <- lapply(eic$scan, function(masterScan)
			{
				masterHeader <- headerData[headerData$acquisitionNum == masterScan,]
				childHeaders <- headerData[(headerData$precursorScanNum == masterScan) 
								& (headerData$precursorMZ > mz - limit.coarse) 
								& (headerData$precursorMZ < mz + limit.coarse) ,]
				childScans <- childHeaders$acquisitionNum
				
				msPeaks <- mzR::peaks(msRaw, masterHeader$seqNum)
				# if deprofile option is set: run deprofiling
				deprofile.setting <- deprofile
				if(!is.na(deprofile.setting))
					msPeaks <- deprofile.scan(
							msPeaks, method = deprofile.setting, noise = NA, colnames = FALSE
					)
				colnames(msPeaks) <- c("mz","int")
				msmsPeaks <- lapply(childHeaders$seqNum, function(scan)
						{
							pks <- mzR::peaks(msRaw, scan)
							if(!is.na(deprofile.setting))
							{								
								pks <- deprofile.scan(
										pks, method = deprofile.setting, noise = NA, colnames = FALSE
								)
							}
							colnames(pks) <- c("mz","int")
							return(pks)
						}
				)
				return(list(
								foundOK = TRUE,
								parentScan = masterScan,
								parentHeader = masterHeader,
								childScans = childScans,
								childHeaders= childHeaders,
								parentPeak=msPeaks,
								peaks=msmsPeaks
						#xset=xset#,
						#msRaw=msRaw
						))
			})
	names(spectra) <- eic$acquisitionNum
	return(spectra)
}

#' @export
findMsMsHR.direct <- function(msRaw, cpdID, mode = "pH", confirmMode = 0, useRtLimit = TRUE, 
			ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
			mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
			fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
			rtMargin = getOption("RMassBank")$rtMargin,
			deprofile = getOption("RMassBank")$deprofile)
{
  # for finding the peak RT: use the gauss-fitted centwave peak
  # (centroid data converted with TOPP is necessary. save as
  # mzData, since this is correctly read :P)
  #xset <- xcmsSet(fileName, method="centWave",ppm=5, fitgauss=TRUE)

  # find cpd m/z
  mzLimits <- findMz(cpdID, mode)
  mz <- mzLimits$mzCenter
  limit.fine <- ppm(mz, ppmFine, p=TRUE)
  if(!useRtLimit)
	  rtLimits <- NA
  else
  {
	  dbRt <- findRt(cpdID)
	  rtLimits <- c(dbRt$RT - rtMargin, dbRt$RT + rtMargin) * 60
  }
  spectra <- findMsMsHR.mass(msRaw, mz, mzCoarse, limit.fine, rtLimits, confirmMode + 1,NA
  	,fillPrecursorScan, deprofile)
  # check whether a) spectrum was found and b) enough spectra were found
  if(length(spectra) < (confirmMode + 1))
    sp <- list(foundOK = FALSE)
  else
    sp <- spectra[[confirmMode + 1]]
  
  sp$mz <- mzLimits
  sp$id <- cpdID
  sp$formula <- findFormula(cpdID)
  return(sp)
}

#' Read in mz-files using XCMS
#' 
#' Picks peaks from mz-files and returns the pseudospectra that CAMERA creates with the help of XCMS
#'
#' @usage findMsMsHRperxcms.direct(fileName, cpdID, mode="pH", findPeaksArgs = NULL, plots = FALSE) 
#' @param fileName The path to the mz-file that should be read
#' @param cpdID The compoundID of the compound that has been used for the file
#' @param mode The ionization mode that has been used for the spectrum represented by the peaklist
#' @param findPeaksArgs A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param plots A parameter that determines whether the spectra should be plotted or not
#' @return The \code{msmsWorkspace} with the additional peaklist added to the right spectrum
#' @seealso \code{\link{msmsWorkflow}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		fileList <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin", full.names=TRUE)[3]
#'		loadList(system.file("XCMSinput/compoundList.csv",package="RMassBank"))
#'      psp <- findMsMsHRperxcms.direct(fileList,2184)
#' }
#' @export
findMsMsHRperxcms.direct <- function(fileName, cpdID, mode="pH", findPeaksArgs = NULL, plots = FALSE) {
	
	require(CAMERA)
	require(xcms)
	parentMass <- findMz(cpdID)$mzCenter
	RT <- findRt(cpdID)$RT * 60
	mzabs <- 0.1
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(peaks(xa@xcmsSet)[x, "rt"])})
	}
	##
	## MS
	##
	
	
	
	##
	## MSMS
	##
	xrmsms <- xcmsRaw(fileName, includeMSn=TRUE)
	print("File read")
	## Where is the wanted isolation ?
	precursorrange <- range(which(xrmsms@msnPrecursorMz == parentMass)) ## TODO: add ppm one day

	## Fake MS1 from MSn scans
	## xrmsmsAsMs <- msn2xcmsRaw(xrmsms)
	
	xrs <- split(msn2xcmsRaw(xrmsms), f=xrmsms@msnCollisionEnergy)
	## Fake s simplistic xcmsSet
	setReplicate <- xcmsSet(files=fileName, method="MS1")
	xsmsms <- as.list(replicate(length(xrs),setReplicate))
	candidates <- list()
	anmsms <- list()
	psp <- list()
	spectra <- list()
	whichmissing <- vector()
	
	for(i in 1:length(xrs)){
		peaks(xsmsms[[i]]) <- do.call(findPeaks,c(findPeaksArgs, object = xrs[[i]]))

                if (nrow(peaks(xsmsms[[i]])) == 0) {
                  spectra[[i]] <- matrix(0,2,7)
                  next 
                }
                
		## Get pspec 
		pl <- peaks(xsmsms[[i]])[,c("mz", "rt"), drop=FALSE]

		## Best: find precursor peak
		candidates[[i]] <- which( pl[,"mz", drop=FALSE] < parentMass + mzabs & pl[,"mz", drop=FALSE] > parentMass - mzabs
						& pl[,"rt", drop=FALSE] < RT * 1.1 & pl[,"rt", drop=FALSE] > RT * 0.9 )
                
		print(paste("Candidates:",candidates[[i]]))
		anmsms[[i]] <- xsAnnotate(xsmsms[[i]])
		anmsms[[i]] <- groupFWHM(anmsms[[i]])

        if(length(candidates[[i]]) > 0)
		closestCandidate <- which.min (abs( RT - pl[candidates[[i]], "rt", drop=FALSE]  ) )
        else(closestCandidate <- 1)
		## Now find the pspec for compound
		
		psp[[i]] <- which(sapply(anmsms[[i]]@pspectra, function(x) {candidates[[i]][closestCandidate] %in% x}))       
		print(paste("Pseudospectra:",psp[[i]]))

		## 2nd best: Spectrum closest to MS1
		##psp <- which.min( abs(getRT(anmsms) - actualRT))

		## 3rd Best: find pspec closest to RT from spreadsheet
		##psp <- which.min( abs(abs(getRT(anmsms) - RT) )
		if((plots == TRUE) && (length(psp[[i]]) > 0)){
			plotPsSpectrum(anmsms[[i]], psp[[i]], log=TRUE,  mzrange=c(0, findMz(cpdID)[[3]]), maxlabel=10)
		}
		if(length(psp[[i]]) != 0){
		spectra[[i]] <- getpspectra(anmsms[[i]], psp[[i]])
		} else {whichmissing <- c(whichmissing,i)}
	}
	if(length(spectra) != 0){
		for(i in whichmissing){
			spectra[[i]] <- matrix(0,2,7)
		}
	}
	return(spectra)
}

# Finds the EIC for a mass trace with a window of x ppm.
# (For ppm = 10, this is +5 / -5 ppm from the non-recalibrated mz.)
#' Extract EICs 
#' 
#' Extract EICs from raw data for a determined mass window.
#' 
#' @param msRaw The mzR file handle 
#' @param mz The mass or mass range to extract the EIC for: either a single mass
#' 			(with the range specified by \code{limit} below) or a mass range
#' 			in the form of \code{c(min, max)}. 
#' @param limit If a single mass was given for \code{mz}: the mass window to extract.
#' 			A limit of 0.001 means that the EIC will be returned for \code{[mz - 0.001, mz + 0.001]}.
#' @param rtLimit If given, the retention time limits in form \code{c(rtmin, rtmax)} in seconds.
#' @return A \code{[rt, intensity, scan]} matrix (\code{scan} being the scan number.) 
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findMsMsHR
#' @export
findEIC <- function(msRaw, mz, limit = NULL, rtLimit = NA)
{
	# calculate mz upper and lower limits for "integration"
	if(all(c("mzMin", "mzMax") %in% names(mz)))
		mzlimits <- c(mz$mzMin, mz$mzMax)
	else
		mzlimits <- c(mz - limit, mz + limit)
	# Find peaklists for all MS1 scans
	headerData <- as.data.frame(header(msRaw))
	# If RT limit is already given, retrieve only candidates in the first place,
	# since this makes everything much faster.
	if(all(!is.na(rtLimit)))
		headerMS1 <- headerData[
				(headerData$msLevel == 1) & (headerData$retentionTime >= rtLimit[[1]])
						& (headerData$retentionTime <= rtLimit[[2]])
				,]
	else
		headerMS1 <- headerData[headerData$msLevel == 1,]
	pks <- mzR::peaks(msRaw, headerMS1$seqNum)
	# Sum intensities in the given mass window for each scan
	pks_t <- unlist(lapply(pks, function(peaktable)
						sum(peaktable[which((peaktable[,1] >= mzlimits[[1]]) & (peaktable[,1] <= mzlimits[[2]])) ,2])))
	rt <- headerMS1$retentionTime
	scan <- headerMS1$acquisitionNum
	return(data.frame(rt = rt, intensity=pks_t, scan=scan))
}

#' Conversion of XCMS-pseudospectra into RMassBank-spectra
#' 
#' Converts a pseudospectrum extracted from XCMS using CAMERA into the msmsWorkspace(at)specs-format that RMassBank uses
#'
#' @usage toRMB(msmsXCMSspecs, cpdID, mode, MS1spec)
#' @param msmsXCMSspecs The compoundID of the compound that has been used for the peaklist
#' @param cpdID The compound ID of the substance of the given spectrum
#' @param mode The ionization mode that has been used for the spectrum
#' @param MS1spec The MS1-spectrum from XCMS, which can be optionally supplied
#' @return One list element of the (at)specs-entry from an msmsWorkspace
#' @seealso \code{\link{msmsWorkspace-class}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		XCMSpspectra <- findmsmsHRperxcms.direct("Glucolesquerellin_2184_1.mzdata", 2184)
#'      wspecs <- toRMB(XCMSpspectra)
#' }
#' @export
toRMB <- function(msmsXCMSspecs = NA, cpdID = NA, mode="pH", MS1spec = NA){
	
	ret <- list()
	ret$mz <- findMz(cpdID,mode=mode)
	ret$id <- cpdID
	ret$formula <- findFormula(cpdID)
	print(paste("Length of msmsXCMSspecs:",length(msmsXCMSspecs)))
	if(length(msmsXCMSspecs) == 0){
		ret$foundOK <- FALSE
		print("blabla")
		return(ret)
	}
	if(is.na(msmsXCMSspecs)){
			stop("You need a readable spectrum!")
	}
	if(is.na(cpdID)){
			stop("Please supply the compoundID!")
	}
	numScan <- length(msmsXCMSspecs)

	
	ret$foundOK <- TRUE
	ret$parentscan <- 1
	ret$parentHeader <- matrix(0, ncol = 20, nrow = 1)
	rownames(ret$parentHeader) <- 1
	colnames(ret$parentHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
									"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
									"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
									"mergedResultStartScanNum", "mergedResultEndScanNum")
	ret$parentHeader[1,1:3] <- 1
	##Write nothing in the parents if there is no MS1-spec
	if(is.na(MS1spec)){
		ret$parentHeader[1,4:20] <- 0
                ret$parentHeader[1,6] <- NA
	} else { ##Else use the MS1spec spec to write everything into the parents
		ret$parentHeader[1,4] <- length(MS1spec[,1])
		ret$parentHeader[1,5] <- 0
		ret$parentHeader[1,6] <- findRt(cpdID)
		ret$parentHeader[1,7] <- MS1spec[which.max(MS1spec[,7]),1]
		ret$parentHeader[1,8] <- max(MS1spec[,7])
		ret$parentHeader[1,9] <- 0
		ret$parentHeader[1,10] <- 0
		ret$parentHeader[1,11] <- min(MS1spec[,1])
		ret$parentHeader[1,12] <- max(MS1spec[,1])
		ret$parentHeader[1,13:20] <- 0 ##Has no precursor and merge is not yet implemented
	}
	ret$parentHeader <- as.data.frame(ret$parentHeader)
	
	##Write the peaks into the childscans
	ret$childScans <- 2:(numScan+1)

	childHeader <- t(sapply(msmsXCMSspecs, function(spec){
		header <- vector()
		header[3] <- 2
		header[4] <- length(spec[,1])
		header[5] <- 0 ##Does this matter?
		header[6] <- median(spec[,4])
		header[7] <- spec[which.max(spec[,7]),1]
		header[8] <- max(spec[,7])
		header[9] <- 0 ##Does this matter?
		header[10] <- 0 ##Does this matter?
		header[11] <- min(spec[,1])
		header[12] <- max(spec[,1]) 
		header[13] <- 1
		header[14] <- findMz(cpdID)[[3]]
		header[15] <- 1 ##Will be changed for different charges
		header[16] <- 0 ##There sadly isnt any precursor intensity to find in the msms-scans. Workaround? msmsXCMS@files[1]
		header[17:20] <- 0 ##Will be changed if merge is wanted
		return(header)
		}))
		childHeader[,1:2] <- 2:(length(msmsXCMSspecs)+1)
	
	
	ret$childHeader <- as.data.frame(childHeader)
	rownames(ret$childHeader) <- 2:(numScan+1)
	colnames(ret$childHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
									"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
									"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
									"mergedResultStartScanNum", "mergedResultEndScanNum")
        if (is.na(ret$parentHeader[1,"retentionTime"])) {
          ## Overwrite MS1 RT with average from MS2 
          ret$parentHeader[1,"retentionTime"] <- median(ret$childHeader[which(ret$childHeader[,"retentionTime"] != 0), "retentionTime"])
        }
	
	ret$parentPeak <- matrix(nrow = 1, ncol = 2)
	colnames(ret$parentPeak) <- c("mz","int")
	ret$parentPeak[1,] <- c(findMz(cpdID,mode=mode)$mzCenter,100)
	ret$peaks <- list()
	ret$peaks <- lapply (msmsXCMSspecs, function(specs){
									peaks <- matrix(nrow = length(specs[,1]), ncol = 2)
									colnames(peaks) <- c("mz","int")
									peaks[,1] <- specs[,1]
									peaks[,2] <- specs[,7]
									return(peaks)
								})
	return(ret)
}

#' Addition of manual peaklists
#' 
#' Adds a manual peaklist in matrix-format
#'
#' @usage addPeaksManually(w, cpdID, handSpec, mode)
#' @param w The msmsWorkspace that the peaklist should be added to.
#' @param cpdID The compoundID of the compound that has been used for the peaklist
#' @param handSpec A peaklist with 2 columns, 
#' @param mode The ionization mode that has been used for the spectrum represented by the peaklist
#' @return The \code{msmsWorkspace} with the additional peaklist added to the right spectrum
#' @seealso \code{\link{msmsWorkflow}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		handSpec <- matrix(0,4,2)
#' 		handSpec[,1] <- c(274.986685367956, 259.012401087427, 95.9493025990907, 96.9573002472772)
#' 		handSpec[,2] <- c(357,761, 2821, 3446)
#' 		addPeaksManually(w, cpdID, handSpec)
#' }
#' @export
addPeaksManually <- function(w, cpdID, handSpec, mode = "pH"){
	childHeaderAddition <- t(sapply(handSpec, function(spec){
			header <- vector()
			header[1:3] <- 2
			header[4] <- length(spec[,1])
			header[5] <- 0 ##Does this matter?
			header[6] <- findRt(cpdID)$RT * 60
			header[7] <- spec[which.max(spec[,2]),1]
			header[8] <- max(spec[,2])
			header[9] <- 0 ##Does this matter?
			header[10] <- 0 ##Does this matter?
			header[11] <- min(spec[,1])
			header[12] <- max(spec[,1])
			header[13] <- 1
			header[14] <- findMz(cpdID)[[3]]
			header[15] <- -1 ##Will be changed for different charges
			header[16] <- 0 ##There sadly isnt any precursor intensity to find in the msms-scans. Workaround? msmsXCMS@files[1]
			header[17:20] <- 0 ##Will be changed if merge is wanted
			return(header)
		}))
	##Set colnames and rownames
	colnames(childHeaderAddition) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
										"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
										"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
										"mergedResultStartScanNum", "mergedResultEndScanNum")
	##Convert the manual peaklists
	peaksHand <- lapply (handSpec, function(specs){
							peaks <- matrix(nrow = length(specs[,1]), ncol = 2)
							colnames(peaks) <- c("mz","int")
							peaks <- specs
							return(peaks)
						})
	
	##Where do the peaks and the header need to be added?
	pos <- sapply(w@specs,function(spec){cpdID %in% spec$id})
	##If the compound for the cpdID isn't in specs yet, add a new spectrum
	if(length(pos) == 0){
		pos <- length(w@specs) + 1
		childHeaderAddition[,1:2] <- 1
		w@specs[[pos]] <- list()
		w@specs[[pos]]$foundOK <- 1
		w@specs[[pos]]$parentscan <- 1
		w@specs[[pos]]$parentHeader <- matrix(0, ncol = 20, nrow = 1)
		rownames(w@specs[[pos]]$parentHeader) <- 1
		colnames(w@specs[[pos]]$parentHeader) <- c("seqNum", "acquisitionNum", "msLevel", "peaksCount", "totIonCurrent", "retentionTime", "basepeakMZ", 
									"basePeakIntensity", "collisionEnergy", "ionisationEnergy", "lowMZ", "highMZ", "precursorScanNum",
									"precursorMZ", "precursorCharge", "precursorIntensity", "mergedScan", "mergedResultScanNum", 
									"mergedResultStartScanNum", "mergedResultEndScanNum")
		w@specs[[pos]]$parentHeader[1,1:3] <- 1
		w@specs[[pos]]$parentHeader[1,4:20] <- 0
		w@specs[[pos]]$parentHeader <- as.data.frame(w@specs[[pos]]$parentHeader)
		w@specs[[pos]]$childScans <- 1
		w@specs[[pos]]$childHeader <- as.data.frame(childHeaderAddition)
		w@specs[[pos]]$parentPeak <- matrix(nrow = 1, ncol = 2)
		colnames(w@specs[[pos]]$parentPeak) <- c("mz","int")
		w@specs[[pos]]$parentPeak[1,] <- c(findMz(cpdID,mode=mode)$mzCenter,100)
		w@specs[[pos]]$peaks <- peaksHand
		w@specs[[pos]]$mz <- findMz(cpdID,mode=mode)
		w@specs[[pos]]$id <- cpdID
		w@specs[[pos]]$formula <- findFormula(cpdID)
	} else { pos <- which(pos)
			w@specs[[pos]]$childHeader <- rbind(w@specs[[pos]]$childHeader,childHeaderAddition)
			w@specs[[pos]]$peaks <- c(w@specs[[pos]]$peaks, peaksHand) }
		
		return(w)
}

#' MassBank-record Addition
#' 
#' Adds the peaklist of a MassBank-Record to the specs of an msmsWorkspace
#'
#' @aliases addMB
#' @usage addMB(w, cpdID, fileName, mode)
#' @param w The msmsWorkspace that the peaklist should be added to.
#' @param cpdID The compoundID of the compound that has been used for the record
#' @param fileName The path to the record
#' @param mode The ionization mode that has been used to create the record
#' @return The \code{msmsWorkspace} with the additional peaklist from the record
#' @seealso \code{\link{addPeaksManually}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		addMB("filepath_to_records/RC00001.txt")
#' }
#' @export
addMB <- function(w, cpdID, fileName, mode){
	mb <- parseMassBank(fileName)
	peaklist <- list()
	peaklist[[1]] <- mb@compiled_ok[[1]][["PK$PEAK"]][,1:2]
	w <- addPeaksManually(w, cpdID, peaklist[[1]], mode)
	return(w)
}
