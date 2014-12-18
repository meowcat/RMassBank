
#library(xcms)

#' Backup \code{msmsWorkflow} results
#' 
#' Writes the results from different \code{msmsWorkflow} steps to a file.
#' 
#' @aliases archiveResults
#' @usage archiveResults(w, fileName, settings = getOption("RMassBank"))
#' @param w The \code{msmsWorkspace} to be saved.
#' @param fileName The filename to store the results under.
#' @param settings The settings to be stored into the msmsWorkspace image.
#' @examples 
#' 
#' 		# This doesn't really make a lot of sense,
#' 		# it stores an empty workspace.
#' 		RmbDefaultSettings()
#' 		w <- newMsmsWorkspace()
#' 		archiveResults(w, "narcotics.RData")
#' 
#' @export
archiveResults <- function(w, fileName, settings = getOption("RMassBank"))
{
  # save the settings into the settings slot
  w@settings <- settings
    # save
  save(w, file=fileName)

}


#' RMassBank mass spectrometry pipeline
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
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param steps Which steps of the workflow to process. See the vignette 
#' 			\code{vignette("RMassBank")} for details.
#' @param confirmMode Defaults to false (use most intense precursor). Value 1 uses
#' 			the 2nd-most intense precursor for a chosen ion (and its data-dependent scans)
#' 			, etc.
#' @param newRecalibration Whether to generate a new recalibration curve (\code{TRUE}, default) or
#' 			to reuse the currently stored curve (\code{FALSE}, useful e.g. for adduct-processing runs.) 
#' @param useRtLimit Whether to enforce the given retention time window.
#' @param archivename The prefix under which to store the analyzed result files.
#' @param readMethod Several methods are available to get peak lists from the files.
#'        Currently supported are "mzR", "xcms", "MassBank" and "peaklist".
#'        The first two read MS/MS raw data, and differ in the strategy 
#'        used to extract peaks. MassBank will read existing records, 
#'        so that e.g. a recalibration can be performed, and "peaklist" 
#'        just requires a CSV with two columns and the column header "mz", "int".
#' @param findPeaksArgs A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param plots A parameter that determines whether the spectra should be plotted or not (This parameter is only used for the xcms-method)
#' @param precursorscan.cf Whether to fill precursor scans. To be used with files which for
#' 		some reasons do not contain precursor scan IDs in the mzML, e.g. AB Sciex converted
#' 		files.
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param analyzeMethod The "method" parameter to pass to \code{\link{analyzeMsMs}}.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param MSe A boolean value that determines whether the spectra were recorded using MSe or not
#' @return The processed \code{msmsWorkspace}.
#' @seealso \code{\link{msmsWorkspace-class}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
msmsWorkflow <- function(w, mode="pH", steps=c(1:8), confirmMode = FALSE, newRecalibration = TRUE, 
		useRtLimit = TRUE, archivename=NA, readMethod = "mzR", findPeaksArgs = NULL, plots = FALSE,
		precursorscan.cf = FALSE,
		settings = getOption("RMassBank"), analyzeMethod = "formula",
		progressbar = "progressBarHook", MSe = FALSE)
{
    .checkMbSettings()
  if(!any(mode %in% c("pH","pNa","pM","mH","mFA","mM",""))) stop(paste("The ionization mode", mode, "is unknown."))
  
  if(!is.na(archivename))
	  w@archivename <- archivename
  
  # Make a progress bar:
  nProg <- 0
  nLen <- length(w@files)
  
  if(readMethod == "simple"){
	##Edit options
	opt <- getOption("RMassBank")
	opt$recalibrator$MS1 <- "recalibrate.identity"
	opt$recalibrator$MS2 <- "recalibrate.identity"
	options(RMassBank=opt)
	##Edit analyzemethod
	analyzeMethod <- "intensity"
  }
  
  # Step 1: acquire all MSMS spectra from files
  if(1 %in% steps)
  {
	if(readMethod == "mzR"){
		nProg <- 0
		message("msmsWorkflow: Step 1. Acquire all MSMS spectra from files")
		pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
		w@specs <- lapply(w@files, function(fileName) {
			
			# Find compound ID
			splitfn <- strsplit(fileName,'_')
			splitsfn <- splitfn[[1]]
			cpdID <- splitsfn[[length(splitsfn)-1]]
			# Retrieve spectrum data
			spec <- findMsMsHR(fileName, cpdID, mode, confirmMode, useRtLimit,
		 		ppmFine = settings$findMsMsRawSettings$ppmFine,
		 		mzCoarse = settings$findMsMsRawSettings$mzCoarse,
		 		fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
		 		rtMargin = settings$rtMargin,
		 		deprofile = settings$deprofile)
			gc()
		
			# Progress:
			nProg <<- nProg + 1
			pb <- do.call(progressbar, list(object=pb, value= nProg))
		
			return(spec)
		} )
		names(w@specs) <- basename(as.character(w@files))
		# close progress bar
		do.call(progressbar, list(object=pb, close=TRUE))
	}
	
	if(readMethod == "xcms"){
		splitfn <- strsplit(w@files,'_')
		cpdIDs <- sapply(splitfn, function(splitted){as.numeric(return(splitted[length(splitted)-1]))})
		files <- list()
		wfiles <- vector()
			for(i in 1:length(unique(cpdIDs))) {
			indices <- sapply(splitfn,function(a){return(unique(cpdIDs)[i] %in% a)})
			files[[i]] <- w@files[indices]
		}
		
		w@files <- sapply(files,function(files){return(files[1])})
		
		for(i in 1:length(unique(cpdIDs))){
			specs <- list()
			for(j in 1:length(files[[i]])){
				specs[[j]] <- findMsMsHRperxcms.direct(files[[i]][j], unique(cpdIDs)[i], mode=mode, findPeaksArgs=findPeaksArgs, plots, MSe=MSe)
			}
			w@specs[[i]] <- unlist(specs,recursive=FALSE)
		}
		w@specs <- unlist(w@specs, recursive=FALSE)
		names(w@specs) <- basename(as.character(w@files))
	}
	
	##if(readMethod == "MassBank"){
	##	for(i in 1:length(w@files)){
	##		w <- addMB(w, w@files[i], mode)
	##	}
	##	names(w@specs) <- basename(as.character(w@files))
	##}
		
	if((readMethod == "peaklist") || (readMethod == "simple")){
		splitfn <- strsplit(basename(w@files), "_")
		cpdIDs <- sapply(splitfn, function(splitted) {
			as.numeric(return(splitted[2]))
		})
		files <- list()
		for(i in 1:length(unique(cpdIDs))) {
			indices <- sapply(splitfn,function(a){return(unique(cpdIDs)[i] %in% a)})
			files[[i]] <- w@files[indices]
		}
			
		peaklist <- list()
			
		for(i in 1:length(w@files)){
			peaklist <- read.csv(w@files[i], header = TRUE)
			w <- addPeaksManually(w, cpdIDs[i], peaklist, mode=mode)
		}
		w@files <- sapply(files,function(file){return(file[1])})
		names(w@specs) <- basename(as.character(w@files))
	}
  }
  # Step 2: first run analysis before recalibration
  if(2 %in% steps)
  {
	  nProg <- 0
	  message("msmsWorkflow: Step 2. First analysis pre recalibration")
	  pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
	  w@analyzedSpecs <- lapply(w@specs, function(spec) {
				  #print(spec$id)
				  s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="preliminary",
						  filterSettings = settings$filterSettings,
						  spectraList = settings$spectraList, method = analyzeMethod)
				  # Progress:
				  nProg <<- nProg + 1
				  pb <- do.call(progressbar, list(object=pb, value= nProg))
				  
				  return(s)
			  })
	  for(f in w@files)
		  w@analyzedSpecs[[basename(as.character(f))]]$name <- basename(as.character(f))
	  do.call(progressbar, list(object=pb, close=TRUE))
  }
  # Step 3: aggregate all spectra
  if(3 %in% steps)
  {
	message("msmsWorkflow: Step 3. Aggregate all spectra")
    w@aggregatedSpecs <- aggregateSpectra(w@analyzedSpecs, addIncomplete=TRUE,
			spectraList = settings$spectraList)
  }
  # Step 4: recalibrate all m/z values in raw spectra
  if(4 %in% steps)
  {
	message("msmsWorkflow: Step 4. Recalibrate m/z values in raw spectra")
	if(newRecalibration)
	{
		recal <- makeRecalibration(w@aggregatedSpecs, mode,
				recalibrateBy = settings$recalibrateBy,
				recalibrateMS1 = settings$recalibrateMS1,
				recalibrator = settings$recalibrator,
				recalibrateMS1Window = settings$recalibrateMS1Window)
		w@rc <- recal$rc
		w@rc.ms1 <- recal$rc.ms1
	}
    w@recalibratedSpecs <- recalibrateSpectra(mode, w@specs, w = w,
			recalibrateBy = settings$recalibrateBy,
			recalibrateMS1 = settings$recalibrateMS1)
  }
  # Step 5: re-analysis on recalibrated spectra
  if(5 %in% steps)
  {
	nProg <- 0
	message("msmsWorkflow: Step 5. Reanalyze recalibrated spectra")
	pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    w@analyzedRcSpecs <- lapply(w@recalibratedSpecs, function(spec) {
      s <- analyzeMsMs(spec, mode=mode, detail=TRUE, run="recalibrated",
			  filterSettings = settings$filterSettings,
			  spectraList = settings$spectraList, method = analyzeMethod)
	  # Progress
	  nProg <<- nProg + 1
	  pb <- do.call(progressbar, list(object=pb, value= nProg))
	  
      return(s)
    }
      )
    for(f in w@files)
      w@analyzedRcSpecs[[basename(as.character(f))]]$name <- basename(as.character(f))
  do.call(progressbar, list(object=pb, close=TRUE))
  }
  # Step 6: aggregate recalibrated results
  if(6 %in% steps)
  {
    message("msmsWorkflow: Step 6. Aggregate recalibrated results")
    w@aggregatedRcSpecs <- aggregateSpectra(w@analyzedRcSpecs, addIncomplete=TRUE,
			spectraList = settings$spectraList)
    if(!is.na(archivename))
      archiveResults(w, paste(archivename, ".RData", sep=''), settings)
    w@aggregatedRcSpecs$peaksUnmatchedC <- 
			cleanElnoise(w@aggregatedRcSpecs$peaksUnmatched,
					settings$electronicNoise, settings$electronicNoiseWidth)
  }
  # Step 7: reanalyze failpeaks for (mono)oxidation and N2 adduct peaks
  if(7 %in% steps)
  {
	message("msmsWorkflow: Step 7. Reanalyze fail peaks for N2 + O")
    w@reanalyzedRcSpecs <- reanalyzeFailpeaks(
			w@aggregatedRcSpecs, custom_additions="N2O", mode=mode,
				filterSettings=settings$filterSettings,
				progressbar=progressbar)
    if(!is.na(archivename))
      archiveResults(w, paste(archivename, "_RA.RData", sep=''), settings)
  }
  # Step 8: heuristic filtering based on peak multiplicity;
  #         creation of failpeak list
  if(8 %in% steps)
  {
	message("msmsWorkflow: Step 8. Peak multiplicity filtering")
    if (is.null(settings$multiplicityFilter)) {
      message("msmsWorkflow: Step 8. Peak multiplicity filtering skipped because multiplicityFilter parameter is not set.")
    } else {
      # apply heuristic filter      
      w@refilteredRcSpecs <- filterMultiplicity(
			w@reanalyzedRcSpecs, archivename, mode, settings$multiplicityFilter )

      if(!is.na(archivename))
        archiveResults(w, paste(archivename, "_RF.RData", sep=''), settings)   
    }
  }
  message("msmsWorkflow: Done.")
  return(w)
}

#' Analyze MSMS spectra
#' 
#' Analyzes MSMS spectra of a compound by fitting formulas to each subpeak.
#' 
#' The analysis function uses Rcdk. Note
#' that in this step, \emph{satellite peaks} are removed by a simple heuristic
#' rule (refer to the documentation of \code{\link{filterPeakSatellites}} for details.)
#' 
#' @usage analyzeMsMs(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
#' 			filterSettings = getOption("RMassBank")$filterSettings,
#' 			spectraList = getOption("RMassBank")$spectraList, method="formula")
#' 
#' 		analyzeMsMs.formula(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
#' 			filterSettings = getOption("RMassBank")$filterSettings,
#' 			spectraList = getOption("RMassBank")$spectraList)
#' 
#' 		analyzeMsMs.intensity(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
#' 			filterSettings = getOption("RMassBank")$filterSettings,
#' 			spectraList = getOption("RMassBank")$spectraList)
#' 
#' @param msmsPeaks A group of parent spectrum and data-dependent MSMS spectra
#' as returned from \code{\link{findMsMsHR}} (refer to the corresponding
#' documentation for the precise format specifications).
#' @param mode Specifies the processing mode, i.e. which molecule species the
#' spectra contain. \code{\var{pH}} (positive H) specifies [M+H]+,
#' \code{\var{pNa}} specifies [M+Na]+, \code{\var{pM}} specifies [M]+,
#' \code{\var{mH}} and \code{\var{mNa}} specify [M-H]- and [M-Na]-,
#' respectively. (I apologize for the naming of \code{\var{pH}} which has
#' absolutely nothing to do with chemical \emph{pH} values.)
#' @param detail Whether detailed return information should be provided
#' (defaults to \code{FALSE}). See below.
#' @param run \code{"preliminary"} or \code{"recalibrated"}. In the
#' \code{preliminary} run, mass tolerance is set to 10 ppm (above m/z 120) and
#' 15 ppm (below m/z 120), the default intensity cutoff is $10^4$ for positive
#' mode (no default cutoff in negative mode), and the column \code{"mz"} from
#' the spectra is used as data source.  In the \code{recalibrated} run, the
#' mass tolerance is set to 5 ppm over the whole mass range, the default cutoff
#' is 0 and the column \code{"mzRecal"} is used as source for the m/z values.
#' Defaults to \code{"preliminary"}.
#' @param filterSettings
#' 		Settings for the filter parameters, by default loaded from the RMassBank settings
#' 		set with e.g. \code{\link{loadRmbSettings}}. Must contain:
#' 		\itemize{
#' 			\item \code{ppmHighMass}, allowed ppm deviation before recalibration
#' 				for high mass range
#' 			\item \code{ppmLowMass}, allowed ppm deviation before recalibration
#' 				for low mass range
#' 			\item \code{massRangeDivision}, division point between high and low mass
#' 				range (before recalibration)
#' 			\item \code{ppmFine}, allowed ppm deviation overall after recalibration
#' 			\item \code{prelimCut}, intensity cutoff for peaks in preliminary run
#' 			\item \code{prelimCutRatio}, relative intensity cutoff for peaks in 
#' 				preliminary run, e.g. 0.01 = 1%
#' 			\item \code{fineCut}, intensity cutoff for peaks in second run
#' 			\item \code{fineCutRatio}, relative intensity cutoff for peaks in 
#' 				second run
#' 			\item \code{specOkLimit}, minimum intensity of base peak for spectrum
#' 				to be accepted for processing
#' 			\item \code{dbeMinLimit}, minimum double bond equivalent for accepted
#' 				molecular subformula.
#' 			\item \code{satelliteMzLimit}, for satellite peak filtering 
#' 				(\code{\link{filterPeakSatellites}}: mass window to use for satellite
#' 				removal
#' 			\item \code{satelliteIntLimit}, the relative intensity below which to 
#' 				discard "satellites". (refer to  \code{\link{filterPeakSatellites}}).
#' 	}
#' @param spectraList The list of MS/MS spectra present in each data block. As also
#' 		defined in the settings file.  
#' @param method Selects which function to actually use for data evaluation. The default
#' 		"formula" runs a full analysis via formula assignment to fragment peaks. The
#' 		alternative setting "intensity" calls a "mock" implementation which circumvents
#' 		formula assignment and filters peaks purely based on intensity cutoffs and the
#' 		satellite filtering. (In this case, the ppm and dbe related settings in filterSettings
#' 		are ignored.)
#' @return \item{list("foundOK")}{
#'  	Boolean. Whether or not child spectra are
#' 		present for this compound (inherited from \code{msmsdata}).}
#' \item{list("mzrange")}{ 
#' 		The maximum m/z range over all child spectra.}
#' \item{list("id")}{ 
#' 		The compound ID (inherited from \code{msmsdata})}
#' \item{list("mode")}{
#' 		processing mode} 
#' \item{list("parentHeader")}{ 
#' 		Parent spectrum header data (ex \code{msmsdata})} 
#' \item{list("parentMs")}{ 
#' 		Parent spectrum (ex \code{msmsdata}) in matrix format} 
#' \item{list("msmsdata")}{
#'		 Analysis results for all child spectra: \itemize{
#' 			 \item\code{specOK} Boolean. Whether or not the spectrum contains
#' 			any useful peaks. If \code{specOK = FALSE}, all other information
#' 			(except scan info and compound ID) may be missing!
#' 			\item\code{parent} Parent mass and formula in a one-row data frame
#' 			format. Currently rather obsolete, originally contained data from 
#' 			MolgenMsMs results.  
#' 			\item \code{childFilt} Annotated peaks of the MSMS spectrum (after 
#' 			filtering by accuracy) 
#'			\item \code{childRaw} Raw (\code{mz, int}) spectrum before any 
#' 			treatment. (With recalibrated data, this is (\code{mz, int, mzRecal}).  
#' 			}
#'  		For \code{detail = TRUE}, additionally:
#' 			\itemize{
#'			\item\code{childRawLow} Peaks cut away because of low (absolute or
#' 			relative) intensity
#' 			\item\code{childRawSatellite} Peaks cut away as"satellites" 
#' 			\item\code{childRawOK} Peaks after cutting away low/satellite
#' 			peaks. Used for further analysis steps 
#' 			\item\code{child} Annotated peaks of the MSMS spectrum before filtering 
#' 			by accuracy 
#' 			\item \code{childBad} Annotated peaks of the MSMS spectrum which didn't
#' 			pass the accuracy threshold
#' 			\item\code{childUnmatched} Peaks of the MSMS spectrum with no annotated
#' 			formula 
#' }}
#' @aliases analyzeMsMs analyzeMsMs.formula analyzeMsMs.intensity
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{filterLowaccResults}},
#' \code{\link{filterPeakSatellites}}, \code{\link{reanalyzeFailpeaks}}
#' @examples
#' 
#' 	\dontrun{analyzed <- analyzeMsMs(spec, "pH", TRUE)}
#' 
#' @export
analyzeMsMs <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
		filterSettings = getOption("RMassBank")$filterSettings,
		spectraList = getOption("RMassBank")$spectraList, method="formula")
{
	.checkMbSettings()
	
	if(msmsPeaks$foundOK == FALSE)
		return(list(foundOK = FALSE, id=msmsPeaks$id))
	
	if(method=="formula")
	{
		return(analyzeMsMs.formula(msmsPeaks, mode, detail, run, filterSettings,
						spectraList
						))
	}
	else if(method == "intensity")
	{
		return(analyzeMsMs.intensity(msmsPeaks, mode, detail, run, filterSettings,
						spectraList
				))
	}
}


#' @export
analyzeMsMs.formula <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
			filterSettings = getOption("RMassBank")$filterSettings,
			spectraList = getOption("RMassBank")$spectraList)
{
  cut <- 0
  cut_ratio <- 0
  if(run=="preliminary")
  {
    mzColname <- "mz"
    filterMode <- "coarse"
	cut <- filterSettings$prelimCut
    if(is.na(cut))
    {
      if(mode %in% c("pH", "pM", "pNa"))
        cut <- 1e4
      else if(mode %in% c("mH", "mFA","mM"))
        cut <- 0
	  else stop(paste("The ionization mode", mode, "is unknown."))
    }
	cutRatio <- filterSettings$prelimCutRatio
  }
  else
  {
    mzColname <- "mzRecal"
    filterMode <- "fine"
	cut <- filterSettings$fineCut
	cut_ratio <- filterSettings$fineCutRatio
    if(is.na(cut)) cut <- 0
  }

  # find whole spectrum of parent peak, so we have reasonable data to feed into
  # MolgenMsMs
  parentSpectrum <- msmsPeaks$parentPeak


  # Check whether the spectra can be fitted to the spectra list correctly!
  if(nrow(msmsPeaks$childHeaders) != length(spectraList))
  {
    warning(paste0(
            "The spectra count of the substance ", msmsPeaks$id, " (", nrow(msmsPeaks$childHeaders), " spectra) doesn't match the provided spectra list (", length(spectraList), " spectra)."
                ))
    return(list(specOK=FALSE))
    
  }
  
  # On each spectrum the following function analyzeTandemShot will be applied.
  # It takes the raw peaks matrix as argument (mz, int) and processes the spectrum by
  # filtering out low-intensity (<1e4) and shoulder peaks (deltam/z < 0.5, intensity
  # < 5%) and subsequently matching the peaks to formulas using Rcdk, discarding peaks
  # with insufficient match accuracy or no match.
  analyzeTandemShot <- function(shot_mat)
  {
    shot <- as.data.frame(shot_mat)
    shot_orig <- shot
    # Filter out low intensity peaks:
    shot_lo <- shot[(shot$int < cut) | (shot$int < max(shot$int)*cut_ratio),]
    shot <- shot[(shot$int >= cut) & (shot$int > max(shot$int) * cut_ratio),]
    shot_full <- shot
    
    # Is there still anything left?
    if(nrow(shot)==0)
      return(list(specOK=FALSE))
    
    # Filter out satellite peaks:
    shot <- filterPeakSatellites(shot, filterSettings)
    shot_satellite_n <- setdiff(row.names(shot_full), row.names(shot))
    shot_satellite <- shot_full[shot_satellite_n,]

    # Is there still anything left?
    if(nrow(shot)==0)
      return(list(specOK=FALSE))
    
    if(max(shot$int) < as.numeric(filterSettings$specOkLimit))
      return(list(specOK=FALSE))
    # Crop to 4 digits (necessary because of the recalibrated values)
    shot[,mzColname] <- round(shot[,mzColname], 5)
    
	
	# here follows the Rcdk analysis
	#------------------------------------
	parentPeaks <- data.frame(mzFound=msmsPeaks$mz$mzCenter, 
			formula=msmsPeaks$formula,
			dppm=0,
			x1=0,x2=0,x3=0)
	
	# define the adduct additions
	if(mode == "pH") {
		allowed_additions <- "H"
		mode.charge <- 1
	} else if(mode == "pNa") {
		allowed_additions <- "Na"
		mode.charge <- 1
	} else if(mode == "pM") {
		allowed_additions <- ""
		mode.charge <- 1
	} else if(mode == "mM") {
		allowed_additions <- ""
		mode.charge <- -1
	} else if(mode == "mH") {
		allowed_additions <- "H-1"
		mode.charge <- -1
	} else if(mode == "mFA") {
		allowed_additions <- "C2H3O2"
		mode.charge <- -1
	} else {
          stop("mode = \"", mode, "\" not defined")
        }
    
	
	# the ppm range is two-sided here.
	# The range is slightly expanded because dppm calculation of
	# generate.formula starts from empirical mass, but dppm cal-
	# culation of the evaluation starts from theoretical mass.
	# So we don't miss the points on 'the border'.
	
	if(run=="preliminary")
		ppmlimit <- 2 * max(filterSettings$ppmLowMass, filterSettings$ppmHighMass)
	else
		ppmlimit <- 2.25 * filterSettings$ppmFine
	
	parent_formula <- add.formula(msmsPeaks$formula, allowed_additions)
	dbe_parent <- dbe(parent_formula)
	# check whether the formula is valid, i.e. has no negative or zero element numbers.
	#print(parent_formula)
	if(!is.valid.formula(parent_formula))
		return(list(specOK=FALSE))
	limits <- to.limits.rcdk(parent_formula)
	
	
	peakmatrix <- lapply(shot[,mzColname], function(mass) {
				# Circumvent bug in rcdk: correct the mass for the charge first, then calculate uncharged formulae
				# finally back-correct calculated masses for the charge
				mass.calc <- mass + mode.charge * .emass
				peakformula <- tryCatch(generate.formula(mass.calc, ppm(mass.calc, ppmlimit, p=TRUE), 
								limits, charge=0), error=function(e) NA)
				#peakformula <- tryCatch( 
				#  generate.formula(mass, 
				#                   ppm(mass, ppmlimit, p=TRUE),
				#                   limits, charge=1),
				#error= function(e) list())
				if(!is.list(peakformula))
					return(t(c(mzFound=as.numeric(as.character(mass)),
											formula=NA, mzCalc=NA)))
				else
				{
					return(t(sapply(peakformula, function(f)
											{
												mzCalc <- f@mass - mode.charge * .emass 
												c(mzFound=mass,
														formula=f@string, 
														mzCalc=mzCalc)
											})))
				}
			})
	
	childPeaks <- as.data.frame(do.call(rbind, peakmatrix))
	childPeaks$mzFound <- as.numeric(as.character(childPeaks$mzFound))
	childPeaks$formula <- as.character(childPeaks$formula)
	childPeaks$mzCalc <- as.numeric(as.character(childPeaks$mzCalc))
	childPeaks$dppm <- (childPeaks$mzFound / childPeaks$mzCalc - 1) * 1e6
	# delete the NA data out again, because MolgenMsMs doesn't have them
	# here and they will be re-added later
	# (this is just left like this for "historical" reasons)
	childPeaks <- childPeaks[!is.na(childPeaks$formula),]
	# check if a peak was recognized (here for the first time,
	# otherwise the next command would fail)
	if(nrow(childPeaks)==0)
		return(list(specOK=FALSE))
	
	# now apply the rule-based filters to get rid of total junk:
	# dbe >= -0.5, dbe excess over mother cpd < 3
	childPeaks$dbe <- unlist(lapply(childPeaks$formula, dbe))
	#iff_rcdk_pM_eln$maxvalence <- unlist(lapply(diff_rcdk_pM_eln$formula.rcdk, maxvalence))
	childPeaks <- childPeaks[childPeaks$dbe >= filterSettings$dbeMinLimit,] 
		# & dbe < dbe_parent + 3)
	
	# check if a peak was recognized
	if(nrow(childPeaks)==0)
		return(list(specOK=FALSE))
	
	# trim mz to 5 digits
	shot[,mzColname] <- round(shot[,mzColname], 5)
    
    childPeaksInt <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, all.x = TRUE, all.y = FALSE )
    # find the best ppm value
    bestPpm <- aggregate(childPeaksInt$dppm, list(childPeaksInt$mzFound),
                         function(dppm) dppm[[which.min(abs(dppm))]])
    colnames(bestPpm) <- c("mzFound", "dppmBest")
    childPeaksInt <- merge(childPeaksInt, bestPpm, by="mzFound", all.x=TRUE)
    # count formulas found per mass
    countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksInt)
    countFormulas <- colSums(countFormulasTab)
    childPeaksInt$formulaCount <- countFormulas[as.character(childPeaksInt$mzFound)]
    # filter results
    childPeaksFilt <- filterLowaccResults(childPeaksInt, filterMode, filterSettings)
    childPeaksGood <- childPeaksFilt[["TRUE"]]
    childPeaksBad <- childPeaksFilt[["FALSE"]]
    # count formulas within new limits
    # (the results of the "old" count stay in childPeaksInt and are returned
    # in $childPeaks)
    if(!is.null(childPeaksGood))
    {
      countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksGood)
      countFormulas <- colSums(countFormulasTab)
      childPeaksGood$formulaCount <- countFormulas[as.character(childPeaksGood$mzFound)]
    }
    
    childPeaksUnmatched <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, 
                                 all.x = TRUE, all.y = TRUE )
    childPeaksUnmatched$dppmBest <- NA
    childPeaksUnmatched$formulaCount <- 0
    childPeaksUnmatched$good <- FALSE
    childPeaksUnmatched <- childPeaksUnmatched[is.na(childPeaksUnmatched$mzCalc),]
    
    # return list:
    rl <- list(
      specOK = !is.null(childPeaksGood),
      parent = parentPeaks,
      childFilt = childPeaksGood,
      childRaw=shot_orig
      )
    # if "detail" is set to TRUE, return more detailed results including
    # all the deleted peaks and the stages when they were culled
    if(detail)
    {
      rl$childRawLow = shot_lo
      rl$childRawSatellite = shot_satellite
      rl$childRawOK = shot
      rl$child =childPeaksInt
      rl$childBad = childPeaksBad
      rl$childUnmatched = childPeaksUnmatched
    }
    return(rl)
  }
  shots <- lapply(msmsPeaks$peaks, analyzeTandemShot)
  #browser()
  shots <- mapply(function(shot, scan, info)
    {
      shot$scan <- scan
      shot$info <- info
      shot$header <- msmsPeaks$childHeaders[as.character(scan),]
      return(shot)
  }, shots, msmsPeaks$childScans, spectraList, SIMPLIFY=FALSE)
  
  mzranges <- t(sapply(shots, function(p) {return(range(p$childRaw[,mzColname]))}))
  mzmin <- min(mzranges[,1], na.rm=TRUE)
  mzmax <- max(mzranges[,2], na.rm=TRUE)
  
  return(list(
          msmsdata=shots,
          mzrange=c(mzmin, mzmax),
          id=msmsPeaks$id,
          mode=mode,
          parentHeader = msmsPeaks$parentHeader,
          parentMs = msmsPeaks$parentPeak,
		  formula = msmsPeaks$formula,
          foundOK = TRUE))
}



#' @export
analyzeMsMs.intensity <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
		filterSettings = getOption("RMassBank")$filterSettings,
		spectraList = getOption("RMassBank")$spectraList)
{
	cut <- 0
	cut_ratio <- 0
	if(run=="preliminary")
	{
		mzColname <- "mz"
		filterMode <- "coarse"
		cut <- filterSettings$prelimCut
		if(is.na(cut))
		{
			if(mode %in% c("pH", "pM", "pNa"))
				cut <- 1e4
			else if(mode %in% c("mH", "mFA", "mM"))
				cut <- 0
			else stop(paste("The ionization mode", mode, "is unknown."))
		}
		cutRatio <- filterSettings$prelimCutRatio
	}
	else
	{
		mzColname <- "mzRecal"
		filterMode <- "fine"
		cut <- filterSettings$fineCut
		cut_ratio <- filterSettings$fineCutRatio
		if(is.na(cut)) cut <- 0
	}
	
	# find whole spectrum of parent peak, so we have reasonable data to feed into
	# MolgenMsMs
	parentSpectrum <- msmsPeaks$parentPeak
	
	
	
	# On each spectrum the following function analyzeTandemShot will be applied.
	# It takes the raw peaks matrix as argument (mz, int) and processes the spectrum by
	# filtering out low-intensity (<1e4) and shoulder peaks (deltam/z < 0.5, intensity
	# < 5%) and subsequently matching the peaks to formulas using Rcdk, discarding peaks
	# with insufficient match accuracy or no match.
	analyzeTandemShot <- function(shot_mat)
	{
		shot <- as.data.frame(shot_mat)
		shot_orig <- shot
		# Filter out low intensity peaks:
		shot_lo <- shot[(shot$int < cut) | (shot$int < max(shot$int)*cut_ratio),]
		shot <- shot[(shot$int >= cut) & (shot$int > max(shot$int) * cut_ratio),]
		shot_full <- shot
		
		# Is there still anything left?
		if(nrow(shot)==0)
			return(list(specOK=FALSE))
		
		# Filter out satellite peaks:
		shot <- filterPeakSatellites(shot, filterSettings)
		shot_satellite_n <- setdiff(row.names(shot_full), row.names(shot))
		shot_satellite <- shot_full[shot_satellite_n,]
		
		# Is there still anything left?
		if(nrow(shot)==0)
			return(list(specOK=FALSE))
		
		if(max(shot$int) < as.numeric(filterSettings$specOkLimit))
			return(list(specOK=FALSE))
		# Crop to 4 digits (necessary because of the recalibrated values)
		shot[,mzColname] <- round(shot[,mzColname], 5)
		
		parentPeaks <- data.frame(mzFound=msmsPeaks$mz$mzCenter, 
				formula=msmsPeaks$formula,
				dppm=0,
				x1=0,x2=0,x3=0)
		
		childPeaks <- data.frame(mzFound = as.numeric(shot[,mzColname]),
						formula = "",
						mzCalc = as.numeric(shot[,mzColname])
						)
		
		childPeaks$dppm <- 0
		childPeaks <- childPeaks[!is.na(childPeaks$formula),]
		# check if a peak was recognized (here for the first time,
		# otherwise the next command would fail)
		if(nrow(childPeaks)==0)
			return(list(specOK=FALSE))
		
		# write dummy values for everything
		childPeaks$dbe <- 0
		
		# trim mz to 5 digits
		shot[,mzColname] <- round(shot[,mzColname], 5)
		
		childPeaksInt <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, all.x = TRUE, all.y = FALSE )
		# find the best ppm value
		bestPpm <- aggregate(childPeaksInt$dppm, list(childPeaksInt$mzFound),
				function(dppm) dppm[[which.min(abs(dppm))]])
		colnames(bestPpm) <- c("mzFound", "dppmBest")
		childPeaksInt <- merge(childPeaksInt, bestPpm, by="mzFound", all.x=TRUE)
		# count formulas found per mass
		countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksInt)
		countFormulas <- colSums(countFormulasTab)
		childPeaksInt$formulaCount <- countFormulas[as.character(childPeaksInt$mzFound)]
		# filter results
		childPeaksFilt <- filterLowaccResults(childPeaksInt, filterMode, filterSettings)
		childPeaksGood <- childPeaksFilt[["TRUE"]]
		childPeaksBad <- childPeaksFilt[["FALSE"]]
		# count formulas within new limits
		# (the results of the "old" count stay in childPeaksInt and are returned
		# in $childPeaks)
		if(!is.null(childPeaksGood))
		{
			countFormulasTab <- xtabs( ~formula + mzFound, data=childPeaksGood)
			countFormulas <- colSums(countFormulasTab)
			childPeaksGood$formulaCount <- countFormulas[as.character(childPeaksGood$mzFound)]
		}
		
		childPeaksUnmatched <- merge(childPeaks, shot, by.x = "mzFound", by.y = mzColname, 
				all.x = TRUE, all.y = TRUE )
		childPeaksUnmatched$dppmBest <- NA
		childPeaksUnmatched$formulaCount <- 0
		childPeaksUnmatched$good <- FALSE
		childPeaksUnmatched <- childPeaksUnmatched[is.na(childPeaksUnmatched$mzCalc),]
		
		# return list:
		rl <- list(
				specOK = !is.null(childPeaksGood),
				parent = parentPeaks,
				childFilt = childPeaksGood,
				childRaw=shot_orig
		)
		# if "detail" is set to TRUE, return more detailed results including
		# all the deleted peaks and the stages when they were culled
		if(detail)
		{
			rl$childRawLow = shot_lo
			rl$childRawSatellite = shot_satellite
			rl$childRawOK = shot
			rl$child =childPeaksInt
			rl$childBad = childPeaksBad
			rl$childUnmatched = childPeaksUnmatched
		}
		return(rl)
	}
	shots <- lapply(msmsPeaks$peaks, analyzeTandemShot)
	#browser()
	shots <- mapply(function(shot, scan, info)
			{
				shot$scan <- scan
				shot$info <- info
				shot$header <- msmsPeaks$childHeaders[as.character(scan),]
				return(shot)
			}, shots, msmsPeaks$childScans, spectraList, SIMPLIFY=FALSE)
	
	mzranges <- t(sapply(shots, function(p) {return(range(p$childRaw[,mzColname]))}))
	mzmin <- min(mzranges[,1], na.rm=TRUE)
	mzmax <- max(mzranges[,2], na.rm=TRUE)
	
	return(list(
					msmsdata=shots,
					mzrange=c(mzmin, mzmax),
					id=msmsPeaks$id,
					mode=mode,
					parentHeader = msmsPeaks$parentHeader,
					parentMs = msmsPeaks$parentPeak,
					formula = msmsPeaks$formula,
					foundOK = TRUE))
}


#' Filter peaks with low accuracy
#' 
#' Filters a peak table (with annotated formulas) for accuracy. Low-accuracy
#' peaks are removed.
#' 
#' In the \code{coarse} mode, mass tolerance is set to 10 ppm (above m/z 120)
#' and 15 ppm (below m/z 120). This is useful for formula assignment before
#' recalibration, where a wide window is desirable to accomodate the high mass
#' deviations at low m/z values, so we get a nice recalibration curve.
#' 
#' In the \code{fine} run, the mass tolerance is set to 5 ppm over the whole
#' mass range. This should be applied after recalibration.
#' 
#' @usage filterLowaccResults(peaks, mode="fine", filterSettings  = getOption("RMassBank")$filterSettings)
#' @param peaks A data frame with at least the columns \code{mzFound} and
#' \code{dppm}.
#' @param mode \code{coarse} or \code{fine}, see below.
#' @param filterSettings Settings for filtering. For details, see documentation of
#' 		\code{\link{analyzeMsMs}}
#' @return A \code{list(TRUE = goodPeakDataframe, FALSE = badPeakDataframe)} is
#' returned: A data frame with all peaks which are "good" is in
#' \code{return[["TRUE"]]}.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{filterPeakSatellites}}
#' @examples
#' 
#' # from analyzeMsMs:
#' \dontrun{childPeaksFilt <- filterLowaccResults(childPeaksInt, filterMode)}
#' 
#'
filterLowaccResults <- function(peaks, mode="fine", filterSettings  = getOption("RMassBank")$filterSettings)
{
  # Check if filter settings are properly set, otherwise use defaults 
  if(is.null(filterSettings))
  {
	  filterSettings <- list(
			ppmHighMass = 10,
	  		ppmLowMass = 15,
	  		massRangeDivision = 120,
	  		ppmFine = 5)
  }
	  
  peaks$good = TRUE
  # coarse mode: to use for determinating the recalibration function
  if(mode=="coarse")
  {
    if(nrow(peaks[which(abs(peaks$dppm) > filterSettings$ppmHighMass),])>0)
    	peaks[which(abs(peaks$dppm) > filterSettings$ppmHighMass), "good"] <- FALSE
	if(nrow(peaks[which(peaks$mzFound > filterSettings$massRangeDivision & abs(peaks$dppm) > filterSettings$ppmLowMass),])>0)
    	peaks[which(peaks$mzFound > filterSettings$massRangeDivision & abs(peaks$dppm) > filterSettings$ppmLowMass), "good"] <- FALSE
  }
  # fine mode: for use after recalibration
  else
  {
	if(nrow(peaks[which(abs(peaks$dppm) > filterSettings$ppmFine),]) > 0)
    	peaks[which(abs(peaks$dppm) > filterSettings$ppmFine), "good"] <- FALSE
  }
  return(split(peaks, peaks$good))
}

#' Aggregate analyzed spectra
#' 
#' Groups an array of analyzed spectra and creates aggregated peak tables
#' 
#' \code{\var{addIncomplete}} is relevant for recalibration. For recalibration,
#' we want to use only high-confidence peaks, therefore we set
#' \code{\var{addIncomplete}} to \code{FALSE}. When we want to generate a peak
#' list for actually generating MassBank records, we want to include all peaks
#' into the peak tables.
#' 
#' @usage aggregateSpectra(spec,  addIncomplete=FALSE, spectraList = getOption("RMassBank")$spectraList)
#' @param spec The set of spectra to aggregate
#' @param addIncomplete Whether or not the peaks from incomplete files (files
#' for which less than the maximal number of spectra are present)
#' @param spectraList The list of MS/MS spectra present in each data block. As also
#' 		defined in the settings file.  
#' @return 
#' \item{foundOK }{ A numeric vector with the compound IDs of all files
#' for which spectra were found. \code{names(foundOK)} are the filenames.}
#' \item{foundFail }{ A numeric vector with the compound IDs of all files for
#' which no spectra were found. \code{names(foundOK)} are the filenames.}
#' \item{spectraFound }{ A numeric vector indicated the number of found spectra
#' per compound} 
#' \item{specFound }{ A list of processed spectral data for all
#' 		compounds with at least 1 found spectrum, as returned by
#' 		\code{\link{analyzeMsMs}}.} 
#' \item{specEmpty }{ A list of (not-really-)processed spectral data for
#'  	compounds without spectra.}
#' \item{specComplete }{ A list of processed spectral data for all compounds
#' 		with the full spectrum count (i.e. \code{length(getOption("RMassBank")$spectraList)}
#'  	spectra.)  As such, \code{specComplete} is a subset of \code{specFound}.}
#' \item{specIncomplete}{ A list of processed spectral data for all compounds 
#' 		with incomplete spectrum count. The complement to \code{specComplete}.}
#' \item{peaksMatched }{ A dataframe of all peaks with a matched formula, 
#' 		which survived the elimination criteria. } 
#' \item{peaksUnmatched }{ A dataframe of all peaks without a matched formula, 
#' 		or with a formula which failed the filter criteria.} 
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{analyzeMsMs}}
#' @examples
#' 
#' ## As used in the workflow:
#' \dontrun{%
#' 	analyzedRcSpecs <- lapply(recalibratedSpecs, function(spec)
#' 		analyzeMsMs(spec, mode="pH", detail=TRUE, run="recalibrated", cut=0, cut_ratio=0 ) )
#' 	aggregatedSpecs <- aggregateSpectra(analyzedSpecs)
#' }
#' 
#' @export
aggregateSpectra <- function(spec,  addIncomplete=FALSE, spectraList = getOption("RMassBank")$spectraList)
{
  # Filter all spectra datasets: split into not-identified spectra, 
  # incomplete spectra and complete spectra
  foundOK <- which(lapply(spec, function(f) f$foundOK) == TRUE)
  foundFail <- which(lapply(spec, function(f) f$foundOK) == FALSE)
  resFOK <- spec[foundOK]
  resFFail <- spec[foundFail]
  specOK <- lapply(resFOK, function(f) 
      sum(unlist(lapply(f$msmsdata, function(f) ifelse(f$specOK==TRUE, 1, 0)))))
  # complete spectra have length(spectraList) annotated peaklists
  nspectra <- length(spectraList)
  resSFail <- resFOK[which(specOK!=nspectra)]
  resSOK <- resFOK[which(specOK==nspectra)]
  
  # Aggregate the incomplete spectra into the set?
  if(addIncomplete)
      resOK <- resFOK
  else
      resOK <- resSOK
  
  # Aggregate all identified and unidentified peaks into tables
  # collect all unmatched peaks from all examined samples:
  failpeaks <- lapply(resOK, function(f) lapply(f$msmsdata, function(g)
    {
		
      if(!is.null(g$childBad))
      {
        g$childBad$cpdID <- f$id
        g$childBad$scan <- g$scan
        g$childBad$parentScan <- f$parentHeader$acquisitionNum
      }
      if(!is.null(g$childUnmatched))
      {
		  g$childUnmatched$cpdID <- rep(f$id, nrow(g$childUnmatched))
		  g$childUnmatched$scan <- rep(g$scan, nrow(g$childUnmatched))
		  g$childUnmatched$parentScan <- rep(f$parentHeader$acquisitionNum, nrow(g$childUnmatched))
      }
      return(rbind(g$childBad, g$childUnmatched))
      }))
  # returns a n x length(spectraList) list of lists with the unmatched/excluded peaks from each sample and spectrum
  # aggregate to one big list:
  failpeaks_s <- lapply(failpeaks, function(f) do.call(rbind, f))
  failpeaks_t <- do.call(rbind, failpeaks_s)
  
  # generate list of all winpeaks for recalibration
  winpeaks <- lapply(resOK, function(f) do.call(rbind, lapply(f$msmsdata, function(g)
    {
    if(!is.null(g$childFilt))
    {
      g$childFilt$cpdID <- f$id
      g$childFilt$scan <- g$scan
      g$childFilt$parentScan <- f$parentHeader$acquisitionNum
    }
    return(g$childFilt)
    })))
  winpeaks_t <- do.call(rbind, winpeaks)
  
  # calculate dppm values
  winpeaks_t$dppmRc <- (winpeaks_t$mzFound/winpeaks_t$mzCalc - 1)*1e6

  return(list(
      foundOK = foundOK,
      foundFail = foundFail,
      spectraFound = specOK,
      specFound = resFOK,
      specEmpty = resFFail,
      specComplete = resSOK,
      specIncomplete = resSFail,
      peaksMatched = winpeaks_t,
      peaksUnmatched = failpeaks_t
      ))
}

#' Remove electronic noise
#' 
#' Removes known electronic noise peaks from a peak table
#' 
#' @usage cleanElnoise(peaks, noise=getOption("RMassBank")$electronicNoise,
#' 		width = getOption("RMassBank")$electronicNoiseWidth)
#' @param peaks A data frame with peaks containing at least the columns
#' \code{mzFound}, \code{dppm} and \code{dppmBest}.
#' @param noise A numeric vector of known m/z of electronic noise peaks from the instrument
#' 		Defaults to the	entries in the RMassBank settings.
#' @param width The window for the noise peak in m/z units. Defaults to the entries in 
#' 		the RMassBank settings.
#' @return Returns a dataframe where the rows matching electronic noise
#'		criteria are removed.  %% ...
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}
#' @examples
#' # As used in the workflow:
#' \dontrun{
#' 	    aggregatedRcSpecs$peaksUnmatchedC <- 
#' 		cleanElnoise(aggregatedRcSpecs$peaksUnmatched)	
#' }
#' @export
cleanElnoise <- function(peaks, noise=getOption("RMassBank")$electronicNoise,
		width = getOption("RMassBank")$electronicNoiseWidth)
{
      # use only best peaks
      p_best <- peaks[is.na(peaks$dppmBest) | (peaks$dppm == peaks$dppmBest),]
      
      # remove known electronic noise
      p_eln <- p_best
      for(noisePeak in noise)
      {
        p_eln <- p_eln[
				(p_eln$mzFound > noisePeak + width)
                | (p_eln$mzFound < noisePeak - width),]
      }
      return(p_eln)
}

#' Identify intense peaks (in a list of unmatched peaks)
#' 
#' Finds a list of peaks in spectra with a high relative intensity (>10% and
#' 1e4, or >1% and 1e5) to write a list of peaks which must be manually
#' checked.  Peaks orbiting around the parent peak mass (calculated from the
#' compound ID), which are very likely co-isolated substances, are ignored.
#' 
#' 
#' @usage problematicPeaks(peaks_unmatched, peaks_matched, mode = "pH")
#' @param peaks_unmatched Table of unmatched peaks, with at least \code{cpdID,
#' scan, mzFound, int}.
#' @param peaks_matched Table of matched peaks (used for base peak reference),
#' with at least \code{cpdID, scan, int}.
#' @param mode Processing mode (\code{"pH", "pNa"} etc.)
#' @return A filtered table with the potentially problematic peaks, including
#' the precursor mass and MSMS base peak intensity (\code{aMax}) for reference.
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}
#' @examples \dontrun{
#' # As used in the workflow: 
#'     fp_rean <-  problematicPeaks(
#'                 peaksNoformula,
#'                 specs$peaksMatched,
#'                 mode)
#' }
#' @export
problematicPeaks <- function(peaks_unmatched, peaks_matched, mode="pH")
{
  # find spectrum maximum for each peak, and merge into table
  assIntMax <- as.data.frame(aggregate(peaks_matched$int, 
        by=list(peaks_matched$cpdID, peaks_matched$scan), max))
  colnames(assIntMax) <- c("cpdID", "scan", "aMax")
  peaks_unmatched <- merge(peaks_unmatched, assIntMax)
  # which of these peaks are intense?
  p_control <- peaks_unmatched[
  	( (peaks_unmatched$int > 1e5) &
	(peaks_unmatched$int > 0.01*peaks_unmatched$aMax)) 
			| ( (peaks_unmatched$int > 1e4) &
				(peaks_unmatched$int > 0.1* peaks_unmatched$aMax)) ,]
  # find parent m/z to exclude co-isolated peaks
  #p_control$mzCenter <- numeric(nrow(p_control))
  p_control$mzCenter <- as.numeric(
    unlist(lapply(p_control$cpdID, function(id) findMz(id, mode)$mzCenter)) )
  p_control_noMH <- p_control[
		  (p_control$mzFound < p_control$mzCenter - 1) |
		  (p_control$mzFound > p_control$mzCenter + 1),]
  return(p_control_noMH)
}


#' Recalibrate MS/MS spectra
#' 
#' Recalibrates MS/MS spectra by building a recalibration curve of the
#' assigned putative fragments of all spectra in \code{aggregatedSpecs}
#' (measured mass vs. mass of putative associated fragment) and additionally
#' the parent ion peaks. 
#' 
#' Note that the actually used recalibration functions are governed by the
#' general MassBank settings (see \code{\link{recalibrate}}).
#' 
#' If a set of acquired LC-MS runs contains spectra for two different ion types
#' (e.g. [M+H]+ and [M+Na]+) which should both be processed by RMassBank, it is
#' necessary to do this in two separate runs. Since it is likely that one ion type
#' will be the vast majority of spectra (e.g. most in [M+H]+ mode), and only few
#' spectra will be present for other specific adducts (e.g. only few [M+Na]+ spectra),
#' it is possible that too few spectra are present to build a good recalibration curve
#' using only e.g. the [M+Na]+ ions. Therefore we recommend, for one set of LC/MS runs,
#' to build the recalibration curve for one ion type 
#' (\code{msmsWorkflow(mode="pH", steps=c(1:8), newRecalibration=TRUE)})
#' and reuse the same curve for processing different ion types 
#' (\code{msmsWorkflow(mode="pNa", steps=c(1:8), newRecalibration=FALSE)}).
#' This also ensures a consistent recalibration across all spectra of the same batch. 
#' 
#' @usage makeRecalibration(spec, mode, 
#'  	recalibrateBy = getOption("RMassBank")$recalibrateBy,
#' 		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1,
#' 		recalibrator = getOption("RMassBank")$recalibrator,
#' 		recalibrateMS1Window = getOption("RMassBank")$recalibrateMS1Window 
#' 		)
#' 
#'  recalibrateSpectra(mode, rawspec = NULL, rc = NULL, rc.ms1=NULL, w = NULL,
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy,
#' 		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1)
#' 
#'  recalibrateSingleSpec(spectrum, rc, 
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' @aliases makeRecalibration recalibrateSpectra recalibrateSingleSpec
#' @param spec For \code{recalibrateSpectra}: a list of \code{aggregatedSpecs} type
#' 			(i.e. as returned by \code{aggregateSpectra}). 
#' @param spectrum For \code{recalibrateSingleSpec}:
#' 			a matrix with columns \code{mz, int} to be recalibrated.
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param rawspec For \code{recalibrateSpectra}:a \code{list} of \code{specs}-type
#' 			object, i.e. as returned by the \code{\link{findMsMsHR}} function family.
#' 			If empty, no spectra are recalibrated, but the recalibration curve is
#' 			returned.  
#' @param rc,rc.ms1 The recalibration curves to be used in the recalibration.
#' @param recalibrateBy Whether recalibration should be done by ppm ("ppm") or by m/z ("mz").
#' @param recalibrateMS1 Whether MS1 spectra should be recalibrated separately ("separate"),
#' 		together with MS2 ("common") or not at all ("none"). Usually taken from settings.
#' @param recalibrator The recalibrator functions to be used.
#' 		 Refer to \code{\link{recalibrate}} for details. Usually taken from settings.
#' @param recalibrateMS1Window Window width to look for MS1 peaks to recalibrate (in ppm).
#' @param w The \code{msmsWorkspace} to write the calibration to or to get the calibration from.
#' @return \code{makeRecalibration}: a \code{list(rc, rc.ms1)} with recalibration curves
#' 			for the MS2 and MS1 spectra.
#' 
#' 			\code{recalibrateSpectra}: if \code{rawspec} is not \code{NULL}, returns the recalibrated
#' 			spectra in the same structure as the input spectra. Each spectrum matrix has
#' 			an additional column \code{mzRecal} with the recalibrated mass.
#' 
#' 			\code{recalibrateSingleSpec}: a matrix with the single recalibrated spectrum. 
#' 			Column \code{mzRecal} contains the recalibrated value.
#' 
#' @examples \dontrun{ 
#' 			rcCurve <- recalibrateSpectra(aggregatedSpecs, "pH")
#' 			recalibratedSpecs <- recalibrateSpectra(aggregatedSpecs, "pH", specs, w=myWorkspace)
#' 			recalibratedSpecs <- recalibrateSpectra(aggregatedSpecs, "pH", specs,
#' 				rcCurve$rc, rcCurve$rc.ms1)
#' 			s <- matrix(c(100,150,200,88.8887,95.0005,222.2223), ncol=2)
#' 			colnames(s) <- c("mz", "int")
#' 			recalS <- recalibrateSingleSpec(s, rcCurve$rc)
#' 			}
#' 
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export 
makeRecalibration <- function(spec, mode, 
		recalibrateBy = getOption("RMassBank")$recalibrateBy,
		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1,
		recalibrator = getOption("RMassBank")$recalibrator,
		recalibrateMS1Window = getOption("RMassBank")$recalibrateMS1Window 
		)
{
	if(is.null(spec))
		stop("No spectra present to generate recalibration curve.")

	if(length(spec$peaksMatched$formulaCount)==0)
		stop("No peaks matched to generate recalibration curve.")
	
	rcdata <- spec$peaksMatched[spec$peaksMatched$formulaCount==1,,drop=FALSE]
	ms1data <- recalibrate.addMS1data(spec, mode, recalibrateMS1Window)

	if (recalibrateMS1 != "none") {
          ## Add m/z values from MS1 to calibration datapoints
          rcdata <- rbind(rcdata, ms1data)
        }
        
	rcdata$dmz <- rcdata$mzFound - rcdata$mzCalc
	ms1data$dmz <- ms1data$mzFound - ms1data$mzCalc
	
	if(recalibrateBy == "dppm")
	{
		rcdata$recalfield <- rcdata$dppm
		ms1data$recalfield <- ms1data$dppm
	}
	else
	{
		rcdata$recalfield <- rcdata$dmz
		ms1data$recalfield <- ms1data$dmz
	}
	
	# generate recalibration model
	rc <- do.call(recalibrator$MS2, list(rcdata)) 
	if(recalibrateMS1 == "separate")
		rc.ms1 <- do.call(recalibrator$MS1, list(ms1data)) 
	else
		rc.ms1 <- rc
	
	# plot the model
	par(mfrow=c(2,2))
	if(nrow(rcdata)>0)
		plotRecalibration.direct(rcdata, rc, rc.ms1, "MS2", 
				range(spec$peaksMatched$mzFound),
				recalibrateBy)
	if(nrow(ms1data)>0)
		plotRecalibration.direct(ms1data, rc, rc.ms1, "MS1",
				range(ms1data$mzFound),
				recalibrateBy)
	# Return the computed recalibration curves
	return(list(rc=rc, rc.ms1=rc.ms1))
}



#' Plot the recalibration graph.
#' 
#' @aliases plotRecalibration plotRecalibration.direct
#' @usage plotRecalibration(w, recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' 
#' 		plotRecalibration.direct(rcdata, rc, rc.ms1, title, mzrange,
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' 
#' @param w The workspace to plot the calibration graph from
#' @param rcdata A data frame with columns \code{recalfield} and \code{mzFound}.
#' @param rc Predictor for MS2 data
#' @param rc.ms1 Predictor for MS1 data
#' @param title Prefix for the graph titles
#' @param mzrange m/z value range for the graph
#' @param recalibrateBy Whether recalibration was done by ppm ("ppm") or by m/z ("mz").
#' 		Important only for graph labeling here.
#' 
#' @author Michele Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
plotRecalibration <- function(w, recalibrateBy = getOption("RMassBank")$recalibrateBy)
{
	spec <- w@aggregatedSpecs
	
	rcdata <- data.frame(mzFound = w@rc$x, recalfield = w@rc$y)
	ms1data <- data.frame(mzFound = w@rc.ms1$x, recalfield = w@rc.ms1$y)
	
	par(mfrow=c(2,2))
	if(nrow(rcdata)>0)
		plotRecalibration.direct(rcdata, w@rc, w@rc.ms1, "MS2", 
				range(spec$peaksMatched$mzFound),recalibrateBy)
	if(nrow(ms1data)>0)
		plotRecalibration.direct(ms1data, w@rc, w@rc.ms1, "MS1",
				range(ms1data$mzFound),recalibrateBy)
	
}

#' @export 
plotRecalibration.direct <- function(rcdata, rc, rc.ms1, title, mzrange,
		recalibrateBy = getOption("RMassBank")$recalibrateBy
		)
{
	if(recalibrateBy == "dppm")
		ylab.plot <- expression(paste(delta, "ppm"))
	else
		ylab.plot <- expression(paste(delta, "m/z"))	
	
	plot(recalfield ~ mzFound, data=rcdata,
			xlab = "m/z", ylab = ylab.plot, main=paste(title, "scatterplot"))
	RcModelMz <- seq(mzrange[[1]], mzrange[[2]], by=0.2)
	RcModelRecal <- predict(rc, newdata= data.frame(mzFound =RcModelMz))
	RcModelRecalMs1 <- predict(rc.ms1, newdata= data.frame(mzFound =RcModelMz))
	lines(RcModelMz, RcModelRecal, col="blue")
	lines(RcModelMz, RcModelRecalMs1, col="yellow")
	if((length(unique(rcdata$mzFound))>1) & 
			(length(unique(rcdata$recalfield))>1))
	{
		if(require(gplots))
		{
			
			hist2d(rcdata$mzFound, rcdata$recalfield, 
					col=c("white", heat.colors(12)), xlab="m/z", 
					ylab = ylab.plot, main=paste(title, "density"))
			lines(RcModelMz, RcModelRecal, col="blue")
			lines(RcModelMz, RcModelRecalMs1, col="yellow")
		}
		else
		{
			message("Package gplots not installed. The recalibration density plot will not be displayed.")
			message("To install gplots: install.packages('gplots')")
		}
	}
}


#' @export
recalibrateSpectra <- function(mode, rawspec = NULL, rc = NULL, rc.ms1=NULL, w = NULL,
		recalibrateBy = getOption("RMassBank")$recalibrateBy,
		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1)
{
	# Load the recal curves from the workspace if one is specified.
  if(!is.null(w))
  {
	  rc <- w@rc
	  rc.ms1 <- w@rc.ms1
  }
  if(is.null(rc) || is.null(rc.ms1))
	  stop("Please specify the recalibration curves either via workspace (w) or via parameters rc, rc.ms1.")

  # Do the recalibration
  if(!is.null(rawspec))
  {
	  # go through all raw spectra and recalculate m/z values
	  recalibratedSpecs <- lapply(rawspec, function(s)
			  {
				  # recalculate tandem spectrum peaks
				  s$peaks <- lapply(s$peaks, function(p)
						  {
							  recalibrateSingleSpec(p, rc, recalibrateBy)
						  })
				  # recalculate MS1 spectrum if required
				  if(recalibrateMS1 != "none")
				  {
					  p <- s$parentPeak;
					  p <- as.data.frame(p)
					  if(nrow(p) > 0)
					  {
						  colnames(p) <- c("mzFound", "int")
						  drecal <- predict(rc.ms1, newdata= p)
						  if(recalibrateBy == "dppm")
							  p$mzRecal <- p$mzFound / ( 1 + drecal/1e6 )
						  else
							  p$mzRecal <- p$mzFound - drecal
						  colnames(p) <- c("mz", "int", "mzRecal")
					  }
					  p <- as.matrix(p)
					  s$parentPeak <- p
				  }
				  return(s)
			  })
	  return(recalibratedSpecs)
  }
  else # no rawspec passed
	  return(list())
}

#' @export
recalibrateSingleSpec <- function(spectrum, rc, 
		recalibrateBy = getOption("RMassBank")$recalibrateBy)
{
	p <- as.data.frame(spectrum)
	if(nrow(p) > 0)
	{
		# Fix the column names so our
		# prediction functions choose the right
		# rows. 
		colnames(p) <- c("mzFound", "int")
		drecal <- predict(rc, newdata= p)
		if(recalibrateBy == "dppm")
			p$mzRecal <- p$mzFound / ( 1 + drecal/1e6 )
		else
			p$mzRecal <- p$mzFound - drecal
		# And rename them back so our "mz" column is
		# called "mz" again
		colnames(p) <- c("mz", "int", "mzRecal")
	}
	p <- as.matrix(p)
	return(p)
}





#' Filter satellite peaks
#' 
#' Filters satellite peaks in FT spectra which arise from FT artifacts and from
#' conversion to stick mode. A very simple rule is used which holds mostly true
#' for MSMS spectra (and shouldn't be applied to MS1 spectra which contain
#' isotope structures...)
#' 
#' The function cuts off all peaks within 0.5 m/z from every peak, in
#' decreasing intensity order, which are below 5% of the referring peak's
#' intensity.  E.g. for peaks m/z=100, int=100; m/z=100.2, int=2, m/z=100.3,
#' int=6, m/z 150, int=10: The most intense peak (m/z=100) is selected, all
#' neighborhood peaks below 5% are removed (in this case, only the m/z=100.2
#' peak) and the next less intense peak is selected. Here this is the m/z=150
#' peak. All low-intensity neighborhood peaks are removed (nothing). The next
#' less intense peak is selected (m/z=100.3) and again neighborhood peaks are
#' cut away (nothing to cut here. Note that the m/z = 100.2 peak was alredy
#' removed.)
#' 
#' @usage filterPeakSatellites(peaks, filterSettings = getOption("RMassBank")$filterSettings)
#' @param peaks A peak dataframe with at least the columns \code{mz, int}. Note
#' that \code{mz} is used even for the recalibrated spectra, i.e. the
#' desatellited spectrum is identical for both the unrecalibrated and the
#' recalibrated spectra.
#' @param filterSettings The settings used for filtering. Refer to \code{\link{analyzeMsMs}}
#' 		documentation for filter settings.
#' @return Returns the peak table with satellite peaks removed.
#' @note This is a very crude rule, but works remarkably well for our spectra.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{filterLowaccResults}}
#' @examples
#' 
#' # From the workflow:
#' \dontrun{
#'     # Filter out satellite peaks:
#'     shot <- filterPeakSatellites(shot)
#'     shot_satellite_n <- setdiff(row.names(shot_full), row.names(shot))
#'     shot_satellite <- shot_full[shot_satellite_n,]
#'     # shot_satellite contains the peaks which were eliminated as satellites.
#' }
#' 
#' @export
filterPeakSatellites <- function(peaks, filterSettings = getOption("RMassBank")$filterSettings)
{
 cutoff_int_limit <- filterSettings$satelliteIntLimit
 cutoff_mz_limit <- filterSettings$satelliteMzLimit
  # Order by intensity (descending)
  peaks_o <- peaks[order(peaks$int, decreasing=TRUE),]
  n <- 1
  # As long as there are peaks left AND the last peak is small enough (relative
  # to selected), move to the next peak
  while(n < nrow(peaks_o))
  {
    if(peaks_o[nrow(peaks_o),"int"] >= cutoff_int_limit *peaks_o[n,"int"])
      break
    # remove all peaks within cutoff_mz_limit (std. m/z = 0.5) which have intensity
    # of less than 5% relative to their "parent" peak
    #
	peaks_l <- peaks_o[ (peaks_o$mz > peaks_o[n,"mz"] - cutoff_mz_limit)
							& (peaks_o$mz < peaks_o[n,"mz"] + cutoff_mz_limit)
							& (peaks_o$int < cutoff_int_limit * peaks_o[n,"int"]),]		 
	peaks_o <- peaks_o[ !((peaks_o$mz > peaks_o[n,"mz"] - cutoff_mz_limit)
								& (peaks_o$mz < peaks_o[n,"mz"] + cutoff_mz_limit)
								& (peaks_o$int < cutoff_int_limit * peaks_o[n,"int"])
								),]		 
    n <- n+1
  }
  return(peaks_o[order(peaks_o$mz),])
}


#' Reanalyze unmatched peaks
#' 
#' Reanalysis of peaks with no matching molecular formula by allowing
#' additional elements (e.g. "N2O").
#' 
#' \code{reanalyzeFailpeaks} examines the \code{unmatchedPeaksC} table in
#' \code{specs} and sends every peak through \code{reanalyzeFailpeak}.
#' 
#' @aliases reanalyzeFailpeaks reanalyzeFailpeak
#' @usage reanalyzeFailpeaks(specs, custom_additions, mode, filterSettings =
#' 				getOption("RMassBank")$filterSettings, progressbar = "progressBarHook")
#' reanalyzeFailpeak(custom_additions, mass, cpdID, counter, pb = NULL, mode,
#' 				filterSettings = getOption("RMassBank")$filterSettings)
#' @param specs An \code{aggregatedRcSpecs} object (after the electronic noise
#' was cleared from the unmatched peaks).
#' @param custom_additions The allowed additions, e.g. "N2O".
#' @param mode Processing mode (\code{"pH", "pNa", "mH"} etc.)
#' @param mass (Usually recalibrated) m/z value of the peak.
#' @param cpdID Compound ID of this spectrum.
#' @param counter Current peak index (used exclusively for the progress
#' indicator)
#' @param pb A progressbar object to display progress on, as passed by
#'  \code{reanalyzeFailpeaks} to \code{reanalyzeFailpeak}. No progress 
#' is displayed if NULL.
#' @param progressbar The progress bar callback to use. Only needed for specialized
#'  applications.	Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param filterSettings Settings for filtering data. Refer to\code{\link{analyzeMsMs}} for settings.
#' @return The returning list contains two tables: 
#' \item{peaksReanalyzed}{All reanalyzed peaks with or without matching formula.}
#' \item{peaksMatchedReanalysis}{Only the peaks with a matched reanalysis
#' formula.}
#' 
#' It would be good to merge the analysis functions of \code{analyzeMsMs} with
#' the one used here, to simplify code changes.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{msmsWorkflow}}
#' @examples
#' 
#' ## As used in the workflow:
#' \dontrun{    
#' 	reanalyzedRcSpecs <- reanalyzeFailpeaks(aggregatedRcSpecs, custom_additions="N2O", mode="pH")
#' # A single peak:
#' reanalyzeFailpeak("N2O", 105.0447, 1234, 1, 1, "pH")
#' }
#' 
#' @export
reanalyzeFailpeaks <- function(specs, custom_additions, mode, filterSettings =
				getOption("RMassBank")$filterSettings, progressbar = "progressBarHook")
{
	
	
  fp <- specs$peaksUnmatchedC
  custom_additions_l <- as.list(rep(x=custom_additions, times=nrow(fp)))
  mode_l <- as.list(rep(x=mode, times=nrow(fp)))
  nLen <- nrow(fp)
  
  pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=max(nLen,1)))
  temp <- data.frame()
  if(nLen == 0)
  {
	  message("reanalyzeFailpeaks: No peaks to reanalyze.")
	  temp <- data.frame(
			  "reanalyzed.formula" = character(),
			  "reanalyzed.mzCalc" = numeric(),
			  "reanalyzed.dppm" = numeric(),
			  "reanalyzed.formulaCount" = numeric(),
			  "reanalyzed.dbe" = numeric())
  }
  else
  {
	  counter <- as.list(1:nrow(fp))
	  # this is the reanalysis step: run reanalyze.failpeak (with the relevant parameters)
	  # on each failpeak.
	  temp <- mapply(reanalyzeFailpeak, custom_additions_l, fp$mzFound, fp$cpdID, counter, 
			  MoreArgs=list(mode=mode, pb=list(hook=progressbar, bar=pb), filterSettings=filterSettings))
	  # reformat the result and attach it to specs
	  temp <- as.data.frame(t(temp))
	  temp <- temp[,c("reanalyzed.formula", "reanalyzed.mzCalc", "reanalyzed.dppm", 
					  "reanalyzed.formulaCount", "reanalyzed.dbe")]	  
  }

  specs$peaksReanalyzed <- cbind(fp, temp)
  
  # Since some columns are in "list" type, they disturb later on.
  # therefore, fix them and make them normal vectors.
  listcols <- unlist(lapply(colnames(specs$peaksReanalyzed), function(col) 
    is.list(specs$peaksReanalyzed[,col])))
  for(col in colnames(specs$peaksReanalyzed)[which(listcols==TRUE)])
    specs$peaksReanalyzed[,col] <- 
      unlist(specs$peaksReanalyzed[,col])
  
  # Now only the matching ones:
  specs$peaksMatchedReanalysis <- specs$peaksReanalyzed[
		  !is.na(specs$peaksReanalyzed$reanalyzed.dppm),]
  
  do.call(progressbar, list(object=pb, close=TRUE))
  return(specs)
}


#' @export
reanalyzeFailpeak <- function(custom_additions, mass, cpdID, counter, pb = NULL, mode,
		filterSettings = getOption("RMassBank")$filterSettings)
{
	# the counter to show the progress
	if(!is.null(pb))
	{
		do.call(pb$hook, list(object=pb$bar, value=counter))
	}
	# here follows the Rcdk analysis
	#------------------------------------
	
	# define the adduct additions
	if(mode == "pH") {
		allowed_additions <- "H"
		mode.charge <- 1
	} else if(mode == "pNa") {
		allowed_additions <- "Na"
		mode.charge <- 1
	} else if(mode == "pM") {
		allowed_additions <- ""
		mode.charge <- 1
	} else if(mode == "mM") {
		allowed_additions <- ""
		mode.charge <- -1
	} else if(mode == "mH") {
		allowed_additions <- "H-1"
		mode.charge <- -1
	} else if(mode == "mFA") {
		allowed_additions <- "C2H3O2"
		mode.charge <- -1
	} else {
          stop("mode = \"", mode, "\" not defined")
        }
	
	# the ppm range is two-sided here.
	# The range is slightly expanded because dppm calculation of
	# generate.formula starts from empirical mass, but dppm cal-
	# culation of the evaluation starts from theoretical mass.
	# So we don't miss the points on 'the border'.
	
	db_formula <- findFormula(cpdID)
	
	ppmlimit <- 2.25 * filterSettings$ppmFine
	parent_formula <- add.formula(db_formula, allowed_additions)
	parent_formula <- add.formula(parent_formula, custom_additions)
	dbe_parent <- dbe(parent_formula)
	# check whether the formula is valid, i.e. has no negative or zero element numbers.
	#print(parent_formula)
	limits <- to.limits.rcdk(parent_formula)        
	
	peakformula <- tryCatch(generate.formula(mass, ppm(mass, ppmlimit, p=TRUE), 
					limits, charge=mode.charge), error=function(e) NA)
	# was a formula found? If not, return empty result
	if(!is.list(peakformula))
		return(as.data.frame(
						t(c(mzFound=as.numeric(as.character(mass)),
										reanalyzed.formula=NA, reanalyzed.mzCalc=NA, reanalyzed.dppm=NA,
										reanalyzed.formulaCount=0,
										reanalyzed.dbe=NA))))
	else # if is.list(peakformula)
	# formula found? then return the one with lowest dppm
	{
		# calculate dppm for all formulas
		peakformula <- sapply(peakformula, function(f)
				{
					l <- list(mzFound=as.numeric(as.character(mass)),
							reanalyzed.formula=as.character(f@string),
							reanalyzed.mzCalc=as.numeric(as.character(f@mass))
					)
					
					return(unlist(l))
				})
		
		# filter out bad dbe stuff
		peakformula <- as.data.frame(t(peakformula))
		# for some reason completely oblivious to me, the columns in peakformula
		# are still factors, even though i de-factored them by hand.
		# Therefore, convert them again...
		peakformula$mzFound <- as.numeric(as.character(peakformula$mzFound))
		peakformula$reanalyzed.formula <- as.character(peakformula$reanalyzed.formula)
		peakformula$reanalyzed.mzCalc <- as.numeric(as.character(peakformula$reanalyzed.mzCalc))
		
		peakformula$reanalyzed.dppm <- (peakformula$mzFound / peakformula$reanalyzed.mzCalc - 1) * 1e6
		peakformula$reanalyzed.formulaCount=nrow(peakformula)
		
		# filter out bad dbe and high ppm stuff          
		peakformula$reanalyzed.dbe <- unlist(lapply(peakformula$reanalyzed.formula, dbe))
		peakformula <- peakformula[(peakformula$reanalyzed.dbe >= filterSettings$dbeMinLimit) 
						& (abs(peakformula$reanalyzed.dppm) < filterSettings$ppmFine),]
		# is there still something left?
		if(nrow(peakformula) == 0)
			return(as.data.frame(
							t(c(mzFound=as.numeric(as.character(mass)),
											reanalyzed.formula=NA, reanalyzed.mzCalc=NA, reanalyzed.dppm=NA,
											reanalyzed.formulaCount=0, reanalyzed.dbe = NA))))
		else
		{
			#update formula count to the remaining formulas
			peakformula$reanalyzed.formulaCount=nrow(peakformula)
			return(peakformula[which.min(abs(peakformula$reanalyzed.dppm)),])
		}
		
	} # endif is.list(peakformula)
      

    
    }

#' Multiplicity filtering: Removes peaks which occur only once in a n-spectra set.
#' 
#' For every compound, every peak (with annotated formula) is compared 
#' across all spectra. Peaks whose formula occurs only once for all collision energies
#' / spectra types, are discarded. This eliminates "stochastic formula hits" of pure
#' electronic noise peaks efficiently from the spectra. Note that in the author's 
#' experimental setup two spectra were recorded at every collision energy,
#' and therefore every peak-formula should appear
#' at least twice if it is real, even if it is by chance a fragment which appears
#' on only one collision energy setting. The function was not tested in a different
#' setup. Therefore, use with a bit of caution.
#' @usage filterPeaksMultiplicity(peaks, formulacol, recalcBest = TRUE)
#' @param peaks A data frame containing all peaks to be analyzed; with at least
#' 			the columns \code{cpdID, scan, mzFound} and one column for the formula
#' 			specified with the \code{formulacol} parameter. 
#' @param formulacol Which column the assigned formula is stored in.
#' @param recalcBest Whether the best formula for each peak should be re-determined.
#' 			This is necessary for results from the ordinary \code{\link{analyzeMsMs}}
#' 			analysis which allows multiple potential formulas per peak - the old best match
#' 			could potentially have been dropped because of multiplicity filtering. For results
#' 			from \code{\link{reanalyzeFailpeak}} this is not necessary, since only one potential
#' 			formula is assigned in this case.
#' @return The peak table is returned, enriched with columns:
#' 			\itemize{
#' 				\item{\code{formulaMultiplicity}}{The # of occurrences of this formula
#' 					in the spectra of its compounds.}
#' 				\item{\code{fM_factor}}{\code{formulaMultiplicity} converted to \code{factor} type
#' 					for use with \code{\link{split}}}
#' 			}
#' @examples \dontrun{
#' 		peaksFiltered <- filterPeaksMultiplicity(aggregatedRcSpecs$peaksMatched, 
#' 			"formula", TRUE)
#' 		peaksOK <- subset(peaksFiltered, formulaMultiplicity > 1)
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
filterPeaksMultiplicity <- function(peaks, formulacol, recalcBest = TRUE)
{
	# create dummy for the case that we have no rows
	multInfo <- data.frame(cpdID = character(), 
			formulacol = character(),
			formulaMultiplicity = numeric())
	# rename (because "formulacol" is not the actually correct name)
	colnames(multInfo) <- c("cpdID", formulacol, "formulaMultiplicity")
	
	if(!is.data.frame(peaks))
	{
		message("filterPeaksMultiplicity: no peaks to aggregate")
		peaks <- cbind(peaks, data.frame(formulaMultiplicity=numeric()))
	}
	else
	{
		# calculate duplicity info
		multInfo <- aggregate(peaks$scan, list(peaks$cpdID, peaks[,formulacol]), FUN=length)
		# just for comparison:
		# nform <- unique(paste(pks$cpdID,pks$formula))
		
		# merge the duplicity info into the peak table
		colnames(multInfo) <- c("cpdID", formulacol, "formulaMultiplicity")
		peaks <- merge(peaks, multInfo)
	}

  # separate log intensity data by duplicity (needs duplicity as a factor)
  # and boxplot
  peaks$fM_factor <- as.factor(peaks$formulaMultiplicity)
  
  # nostalgy: dppmBest first, to compare :)
  # now we prioritize the most frequent formula instead, and only then apply the
  # dppmBest rule
  #pks2 <- subset(pks, dppm==dppmBest)
  
  # split peak intensity by multiplicity
  peakMultiplicitySets <- split(log(peaks$int,10), peaks$fM_factor)
  #boxplot(peakMultiplicitySets)
  # nice plot :)
  #if(length(peakMultiplicitySets) > 0)
  #	q <- quantile(peakMultiplicitySets[[1]], c(0,.25,.5,.75,.95,1))
  pk_data <- lapply(peakMultiplicitySets, length)

  # now by formula, not by peak:
  multInfo$fM_factor <- as.factor(multInfo$formulaMultiplicity)
  # the formulas are split into bins with their multiplicity 
  # (14 bins for our 14-spectra method)
  formulaMultiplicitySets <- split(multInfo[,formulacol], multInfo$fM_factor)
  formulaMultiplicityHist <- lapply(formulaMultiplicitySets, length)

  # if we use recalcBest, then we recalculate which peak in the
  # list was best. We do this for the peaks matched in the first analysis.
  # The peaks from the reanalysis are single anyway and don't get this additional
  # treatment.
  
  if(recalcBest == FALSE)
      return(peaks)
  
  # prioritize duplicate peaks
  # get unique peaks with their maximum-multiplicity formula attached
  best_mult <- aggregate(peaks$formulaMultiplicity, 
                         list(peaks$cpdID, peaks$scan, peaks$mzFound), 
                         max)
  colnames(best_mult) <- c("cpdID", "scan", "mzFound", "bestMultiplicity")
  peaks <- merge(peaks, best_mult)
  peaks <- peaks[peaks$formulaMultiplicity==peaks$bestMultiplicity,]
  
  # now we also have to recalculate dppmBest since the "old best" may have been
  # dropped.
  peaks$dppmBest <- NULL
  bestPpm <- aggregate(peaks$dppm, 
                       list(peaks$cpdID, peaks$scan, peaks$mzFound),
                        function(dppm) dppm[[which.min(abs(dppm))]])
  colnames(bestPpm) <- c("cpdID", "scan", "mzFound", "dppmBest")
  peaks <- merge(peaks, bestPpm)
  pks_best <- peaks[peaks$dppm==peaks$dppmBest,]
  
  # And, iteratively, the multiplicity also must be recalculated, because we dropped
  # some peaks and the multiplicites of some of the formulas will have decreased.
    
  pks_best$formulaMultiplicity <- NULL
  pks_best$bestMultiplicity <- NULL
  multInfo_best <- aggregate(pks_best$scan, 
                             list(pks_best$cpdID, pks_best[,formulacol]),
                             FUN=length)
  colnames(multInfo_best) <- c("cpdID", formulacol, "formulaMultiplicity")
  pks_best <- merge(pks_best, multInfo_best)
  pks_best$fM_factor <- as.factor(pks_best$formulaMultiplicity)
  multInfo_best$fM_factor <- as.factor(multInfo_best$formulaMultiplicity)
  
  formulaMultplicitySets_best <- split(multInfo_best[,formulacol], multInfo_best$fM_factor)
  formulaMultplicityHist_best <- lapply(formulaMultplicitySets_best, length)
  
  peakMultiplicitySets_best <- split(log(pks_best$int,10), pks_best$fM_factor)
  #boxplot(peakMultiplicitySets_best)
  #q <- quantile(peakMultiplicitySets_best[[1]], c(0,.25,.5,.75,.95,1))
  #peakMultiplicityHist_best <- lapply(peakMultiplicitySets_best, length)
  #q
  
  # this returns the "best" peaks (first by formula multiplicity, then by dppm)
  # before actually cutting the bad ones off.
  return(pks_best)  
}


#' filterMultiplicity
#' 
#' Multiplicity filtering: Removes peaks which occur only once in a n-spectra
#' set.
#' 
#' This function executes multiplicity filtering for a set of spectra using the
#' workhorse function \code{\link{filterPeaksMultiplicity}} (see details there)
#' and retrieves problematic filtered peaks (peaks which are of high intensity
#' but were discarded, because either no formula was assigned or it was not
#' present at least 2x), using the workhorse function
#' \code{\link{problematicPeaks}}. The results are returned in a format ready
#' for further processing with \code{\link{mbWorkflow}}.
#' 
#' @usage filterMultiplicity(specs, archivename=NA, mode="pH", recalcBest = TRUE,
#' 		multiplicityFilter = getOption("RMassBank")$multiplicityFilter)
#' @param specs aggregatedSpecs object whose peaks should be filtered
#' @param archivename The archive name, used for generation of
#' archivename_failpeaks.csv
#' @param mode Mode of ion analysis
#' @param recalcBest Boolean, whether to recalculate the formula multiplicity 
#' 		after the first multiplicity filtering step. Sometimes, setting this
#' 		to FALSE can be a solution if you have many compounds with e.g. fluorine
#' 		atoms, which often have multiple assigned formulas per peak and might occasionally
#' 		lose peaks because of that. 
#' @param multiplicityFilter Threshold for the multiplicity filter. If set to 1,
#' 		no filtering will apply (minimum 1 occurrence of peak). 2 equals minimum
#' 		2 occurrences etc. 
#' @return A list object with values: 
#' \item{peaksOK}{ Peaks with >1-fold formula multiplicity from the
#' 		"normal" peak analysis.  } 
#' \item{peaksReanOK}{ Peaks with >1-fold formula multiplicity from
#' 		peak reanalysis.  }
#' \item{peaksFiltered}{ All peaks with annotated formula multiplicity from
#' 		first analysis.  } 
#' \item{peaksFilteredReanalysis}{ All peaks with annotated
#' 		formula multiplicity from peak reanalysis.  } 
#' \item{peaksProblematic}{ Peaks with high intensity which do not match 
#' 		inclusion criteria -> possible false negatives. The list will be
#' 		exported into archivename_failpeaks.csv.
#' }
#' @author Michael Stravs
#' @seealso
#' \code{\link{filterPeaksMultiplicity}},\code{\link{problematicPeaks}}
#' @examples
#' \dontrun{
#'     refilteredRcSpecs <- filterMultiplicity(
#' 			reanalyzedRcSpecs, "myarchive", "pH")
#' }
#' @export
filterMultiplicity <- function(specs, archivename=NA, mode="pH", recalcBest = TRUE,
		multiplicityFilter = getOption("RMassBank")$multiplicityFilter)
{
    # Read multiplicity filter setting
    # For backwards compatibility: If the option is not set, define as 2
    # (which was the behaviour before we introduced the option)
    if(is.null(multiplicityFilter))
      multiplicityFilter <- 2
    
    peaksFiltered <- filterPeaksMultiplicity(specs$peaksMatched,
                                                        "formula", recalcBest)
    peaksFilteredReanalysis <- 
      filterPeaksMultiplicity(specs$peaksMatchedReanalysis, "reanalyzed.formula", FALSE)
    
    peaksNoformula <- specs$peaksReanalyzed[is.na(specs$peaksReanalyzed$reanalyzed.formula),]
    peaksNoformula$formulaMultiplicity <- rep(0, nrow(peaksNoformula))
    peaksNoformula$fM_factor <- as.factor( rep(0, nrow(peaksNoformula)))
    
	# Reorder the columns of peaksNoformula such that they match the columns
	# of peaksFilteredReanalysis; such that rbind gives an identical result
	# when peaksFilteredReanalysis is empty. (Otherwise peaksFilterReanalysis
	# would be dropped as 0x0, and rbind's output column order would be the one originally
	# in peaksNoformula. See ?cbind)
	peaksNoformula <- peaksNoformula[,colnames(peaksFilteredReanalysis)]
	
    # export the peaks which drop through reanalysis or filter criteria
    fp_rean <-  problematicPeaks(
                rbind(
                  peaksFilteredReanalysis[ 
						  peaksFilteredReanalysis$formulaMultiplicity < multiplicityFilter,],
                  peaksNoformula),
                specs$peaksMatched,
                mode)
    fp_mult <-  problematicPeaks(
                peaksFiltered[
						peaksFiltered$formulaMultiplicity < multiplicityFilter,],
                specs$peaksMatched,
                mode)
    fp_mult$good <- NULL
    fp_mult$dppmRc <- NULL
    colnames(fp_rean) <- c("cpdID", "scan", "formula", "mzFound", "formula.old", "mzCalc.old",
                           "dppm.old", "dbe.old", "mz", "int", "dppmBest", "formulaCount.old", 
                           "good.old", "parentScan", 
                           "mzCalc", "dppm",
                           "formulaCount", "dbe", "formulaMultiplicity", "fM_factor", 
                           "aMax", "mzCenter")
    fp_tot <- rbind(fp_rean[,colnames(fp_mult)], fp_mult)
    
    # Reorder and reformat for output
    fp_tot <- fp_tot[with(fp_tot, 
                    order(cpdID, mzCalc, scan)),
               ]

   fp_tot$OK <- rep('', nrow(fp_tot))
   fp_tot$name <- rownames(fp_tot)

   # Select the columns for output into the failpeaks file
    fp_tot <- fp_tot[,c("OK", "name", "cpdID", "scan", "mzFound", "formula", "mzCalc", "dppm", "dbe", "mz", "int",
                 "formulaCount", "parentScan", "aMax", "mzCenter")]
 	# Select the correct precursor scans. This serves to filter the list
	# for the cases where multiple workspaces were combined after step 7
	# with combineMultiplicities.
	# Note that this has drawbacks. Leaving the "duplicates" in would make it more easy
	# to identify legitimate unformulaed peaks. We might experiment by marking them up
	# somehow. 
	precursors <- unlist(lapply(specs$specFound, function(s) s$parentHeader$acquisitionNum))
	fp_tot <- fp_tot[
			fp_tot$parentScan %in% precursors
			,]

    
    peaksOK <- peaksFiltered[
                      peaksFiltered$formulaMultiplicity > (multiplicityFilter - 1),]
    peaksReanOK <- peaksFilteredReanalysis[
						(peaksFilteredReanalysis$formulaMultiplicity > (multiplicityFilter - 1)) &
						!is.na(peaksFilteredReanalysis$reanalyzed.formula),]
    # Kick the M+H+ satellites out of peaksReanOK:
    peaksReanOK$mzCenter <- as.numeric(
      unlist(lapply(peaksReanOK$cpdID, function(id) findMz(id)$mzCenter)) )
    peaksReanOK <- peaksReanOK[
			(peaksReanOK$mzFound < peaksReanOK$mzCenter - 1) |
			(peaksReanOK$mzFound > peaksReanOK$mzCenter + 1) ,]
    
    if(!is.na(archivename))
      write.csv(fp_tot, file=
        paste(archivename,"_Failpeaks.csv", sep=''), row.names=FALSE)
    return(list(peaksOK = peaksOK, peaksReanOK = peaksReanOK, peaksFiltered = peaksFiltered,
                peaksFilteredReanalysis = peaksFilteredReanalysis,
                peaksProblematic = fp_tot))

}

#' Return MS1 peaks to be used for recalibration
#' 
#' Returns the precursor peaks for all MS1 spectra in the \code{spec} dataset
#' with annotated formula to be used in recalibration.
#'  
#' For all spectra in \code{spec$specFound}, the precursor ion is extracted from
#' the MS1 precursor spectrum. All found ions are returned in a data frame with a
#' format matching \code{spec$peaksMatched} and therefore suitable for \code{rbind}ing
#' to the \code{spec$peaksMatched} table. However, only minimal information needed for
#' recalibration is returned. 
#' 
#' @usage  recalibrate.addMS1data(spec,mode="pH", recalibrateMS1Window = 
#' 				getOption("RMassBank")$recalibrateMS1Window)
#' @param spec A \code{aggregatedSpecs}-like object.
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param recalibrateMS1Window Window width to look for MS1 peaks to recalibrate (in ppm).
#' @return A dataframe with columns \code{mzFound, formula, mzCalc, dppm, dbe, int,
#' 		dppmBest, formulaCount, good, cpdID, scan, parentScan, dppmRc}. However,
#' 		columns \code{dbe, int, formulaCount, good, scan, parentScan} do not contain
#' 		real information and are provided only as fillers.
#' @examples \dontrun{
#' # More or less as used in recalibrateSpectra:
#' 		rcdata <- subset(aggregatedSpecs$peaksMatched, formulaCount==1)
#' 		ms1data <- recalibrate.addMS1data(aggregatedSpecs, "pH", 15)
#' 		rcdata <- rbind(rcdata, ms1data)
#'  # ... continue constructing recalibration curve with rcdata
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
recalibrate.addMS1data <- function(spec,mode="pH", recalibrateMS1Window = 
				getOption("RMassBank")$recalibrateMS1Window)
{
	## which_OK <- lapply(validPrecursors, function(pscan)
	##         {
	##             pplist <- as.data.frame(
	##                     mzR::peaks(msRaw, which(headerData$acquisitionNum == pscan)))
	##             colnames(pplist) <- c("mz","int")
	##             pplist <- subset(pplist, mz >= mzLimits$mzMin & mz <= mzLimits$mzMax)
	##             if(nrow(pplist) > 0)
	##                 return(TRUE)
	##             return(FALSE)
	##         })
	
	ms1peaks <- lapply(spec$specFound, function(cpd)
			{
				mzL <- findMz.formula(cpd$formula,mode,recalibrateMS1Window,0)
				mzCalc <- mzL$mzCenter
				ms1 <- as.data.frame(cpd$parentMs)
				pplist <- ms1[(ms1$mz >= mzL$mzMin) & (ms1$mz <= mzL$mzMax),]
				mzFound <- pplist[which.max(pplist$int),"mz"]
				dppmRc <- (mzFound/mzCalc - 1)*1e6
				return(c(
								mzFound = mzFound,
								formula = "",
								mzCalc = mzCalc,
								dppm = dppmRc,
								dbe = 0,
								int = 100,
								dppmBest = dppmRc,
								formulaCount = 1,
								good = TRUE,
								cpdID = cpd$id,
								scan = 0,
								parentScan = 0,
								dppmRc = dppmRc
						))
			})
	ms1peaks <- as.data.frame(do.call(rbind, ms1peaks), stringsAsFactors=FALSE)
	# convert numbers to numeric
	tonum <- c("mzFound", "dppm", "dppmRc", "mzCalc", "dppmBest")
	ms1peaks[,tonum] <- as.numeric(unlist(ms1peaks[,tonum]))
	# throw out NA stuff
	ms1peaks <- ms1peaks[!is.na(ms1peaks$mzFound),]
	return(ms1peaks)
}


# Custom recalibration function: You can overwrite the recal function by
# making any function which takes rcdata$recalfield ~ rcdata$mzFound.
# The settings define which recal function is used
# getOption("RMassBank")$recalibrator = list(
#	MS1 = "recalibrate.loess",
#	MS2 = "recalibrate.loess")

#' Predefined recalibration functions.
#' 
#' Predefined fits to use for recalibration: Loess fit and GAM fit.
#' 
#' \code{recalibrate.loess()} provides a Loess fit (\code{recalibrate.loess}) 
#' to a given recalibration parameter.  
#' If MS and MS/MS data should be fit together, recalibrate.loess 
#' provides good default settings for Orbitrap instruments.
#' 
#' \code{recalibrate.identity()} returns a non-recalibration, i.e. a predictor
#' which predicts 0 for all input values. This can be used if the user wants to
#' skip recalibration in the RMassBank workflow.
#' 
#' #' \code{recalibrate.mean()} and \code{recalibrate.linear()} are simple recalibrations
#' which return a constant shift or a linear recalibration. They will be only useful
#' in particular cases.
#' 
#' \code{recalibrate()} itself is only a dummy function and does not do anything.
#' 
#' Alternatively other functions can be defined. Which functions are used for recalibration
#' is specified by the RMassBank options file. (Note: if \code{recalibrateMS1: common}, the
#' \code{recalibrator: MS1} value is irrelevant, since for a common curve generated with
#' the function specified in \code{recalibrator: MS2} will be used.)
#' 
#' @aliases recalibrate.loess recalibrate recalibrate.identity recalibrate.mean recalibrate.linear
#' @usage recalibrate.loess(rcdata)
#' 
#' 		recalibrate.identity(rcdata)
#' 
#' 		recalibrate.mean(rcdata)
#' 
#' 		recalibrate.linear(rcdata)
#' 
#' @param rcdata A data frame with at least the columns \code{recalfield} and
#' 			\code{mzFound}. \code{recalfield} will usually contain delta(ppm) or
#' 			delta(mz) values and is the target parameter for the recalibration.
#' @return Returns a model for recalibration to be used with \code{predict} and the like.
#' @examples \dontrun{
#' rcdata <- subset(spec$peaksMatched, formulaCount==1)
#' ms1data <- recalibrate.addMS1data(spec, mode, 15)
#' rcdata <- rbind(rcdata, ms1data)
#' rcdata$recalfield <- rcdata$dppm
#' rcCurve <- recalibrate.loess(rcdata)
#' # define a spectrum and recalibrate it
#' s <- matrix(c(100,150,200,88.8887,95.0005,222.2223), ncol=2)
#' colnames(s) <- c("mz", "int")
#' recalS <- recalibrateSingleSpec(s, rcCurve)
#' 
#' Alternative: define an custom recalibrator function with different parameters
#' recalibrate.MyOwnLoess <- function(rcdata)
#' {
#' 	return(loess(recalfield ~ mzFound, data=rcdata, family=c("symmetric"),
#' 					degree = 2, span=0.4))
#' }
#' # This can then be specified in the RMassBank settings file:
#' # recalibrateMS1: common
#' # recalibrator:
#' #    MS1: recalibrate.loess
#' #    MS2: recalibrate.MyOwnLoess")
#' # [...]
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
recalibrate <- function()
{
	return(NA)
}

#' @export
recalibrate.loess <- function(rcdata)
{
  span <- 0.25
  # ex XCMS (permission by Steffen): heuristically decide on loess vs linear
  mingroups <- nrow(rcdata[!is.na(rcdata$mzFound),])
  if(mingroups < 4)
  {
    warning("recalibrate.loess: Not enough data points, omitting recalibration")
    return(recalibrate.identity(rcdata))
  } else if (mingroups*span < 4) {
    span <- 4/mingroups
    warning("recalibrate.loess: Span too small, resetting to ", round(span, 2))
  }
	return(loess(recalfield ~ mzFound, data=rcdata, family=c("symmetric"),
					degree = 1, span=0.25, surface="direct" ))
}

#' @export 
recalibrate.identity <- function(rcdata)
{
	return(lm(recalfield ~ 0, data=rcdata))
}

#' @export 
recalibrate.mean <- function(rcdata)
{
  return(lm(recalfield ~ 1, data=rcdata))
}

#' @export 
recalibrate.linear <- function(rcdata)
{
  return(lm(recalfield ~ mzFound, data=rcdata))
}

#' Standard progress bar hook.
#' 
#' This function provides a standard implementation for the progress bar in RMassBank.
#' 
#' RMassBank calls the progress bar function in the following three ways:
#' \code{pb <- progressBarHook(object=NULL, value=0, min=0, max=LEN)}
#' to create a new progress bar.
#' \code{pb <- progressBarHook(object=pb, value= VAL)}
#' to set the progress bar to a new value (between the set \code{min} and \code{max})
#' \code{progressBarHook(object=pb, close=TRUE)}
#' to close the progress bar. (The actual calls are performed with \code{do.call}, 
#' e.g. 
#' \code{progressbar <- "progressBarHook"
#' pb <- do.call(progressbar, list(object=pb, value= nProg))
#' }. See the source code for details.)
#' 
#' To substitute the standard progress bar for an alternative implementation (e.g. for
#' use in a GUI), the developer can write his own function which behaves in the same way
#' as \code{progressBarHook}, i.e. takes the same parameters and can be called in the 
#' same way. 
#'  
#' @param object An identifier representing an instance of a progress bar. 
#' @param value The new value to assign to the progress indicator
#' @param min The minimal value of the progress indicator
#' @param max The maximal value of the progress indicator
#' @param close If \code{TRUE}, the progress bar is closed.
#' @return Returns a progress bar instance identifier (i.e. an identifier
#' 		which can be used as \code{object} in subsequent calls.) 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
progressBarHook <- function(object = NULL, value = 0, min = 0, max = 100, close = FALSE)
{
	if(is.null(object))
	{
		object <- txtProgressBar(min, max, value, style=3, file=stderr())
	}
	if(close)
		close(object)
	else
	{
		setTxtProgressBar(object, value)
		return(object)
	}
}
