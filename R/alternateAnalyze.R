fragDataIndexing <- function(fragData){
  index <- vector()
  fragDataWholeMass <- floor(fragData$mass*100)
  
  k <- 0
  
  for(i in 1:length(fragDataWholeMass)){
    while(k <= fragDataWholeMass[i]){
      index[k] <- i
      k <- k + 1
    }
  }

  return(index)
}

is.sub.formula <- function(formula, targetAtomList){
  atomList <- formulastring.to.list(formula)
  if(!all(names(atomList) %in% names(targetAtomList))){
    return(FALSE)
  }
  
  for(atom in names(atomList)){
    if(atomList[[atom]] > targetAtomList[[atom]]){
      return(FALSE)
    }
  }
  return(TRUE)
}

newStep2WorkFlow <- function(w, mode="pH", 
                             confirmMode=FALSE, progressbar = "progressBarHook", 
                             settings = getOption("RMassBank"), analyzeMethod="formula", fragdataFile = NA){
  ##Load the fragment data (locally or from RMassBank package)
  fragData <- read.csv(fragdataFile,colClasses = c("character","numeric"))
  
  ##Progress bar
  nLen <- length(w@files)
  nProg <- 0
  log_info("msmsWorkflow: Step 2. First analysis pre recalibration")
  pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
  
  ##Index the fragment data (for time reasons, "which" is very slow for large matrices)
  fragDataIndex <- fragDataIndexing(fragData)
  
  
  
  w@analyzedSpecs <- lapply(w@specs, function(spec) {
    #print(spec$id)
    s <- analyzeMsMs.optimized(spec, mode=mode, detail=TRUE, run="preliminary",
                     filterSettings = settings$filterSettings,
                     spectraList = settings$spectraList, method = analyzeMethod, fragData = fragData, 
                     fragDataIndex = fragDataIndex)
    # Progress:
    nProg <<- nProg + 1
    pb <- do.call(progressbar, list(object=pb, value= nProg))
    
    return(s)
  })
  
  for(f in w@files)
    w@analyzedSpecs[[basename(as.character(f))]]$name <- basename(as.character(f))
  
  suppressWarnings(do.call(progressbar, list(object=pb, close=TRUE)))
  return(w)
}


analyzeMsMs.optimized <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
                        filterSettings = getOption("RMassBank")$filterSettings,
                        spectraList = getOption("RMassBank")$spectraList, method="formula", fragData,
                        fragDataIndex)
{
  .checkMbSettings()
  
  if(msmsPeaks$foundOK == FALSE)
    return(list(foundOK = FALSE, id=msmsPeaks$id))
  
  if(method=="formula")
  {
    return(analyzeMsMs.formula.optimized(msmsPeaks, mode, detail, run, filterSettings,
                               spectraList, fragData, fragDataIndex))
  }
  else if(method == "intensity")
  {
    return(analyzeMsMs.intensity(msmsPeaks, mode, detail, run, filterSettings))
  }
}

analyzeMsMs.formula.optimized <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
                                filterSettings = getOption("RMassBank")$filterSettings,
                                spectraList = getOption("RMassBank")$spectraList, fragData, fragDataIndex)
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
		  adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
		  if(adductProperties$charge > 0) cut <- 1e4
		  if(adductProperties$charge < 0) cut <- 0
		}
		cutRatio <- filterSettings$prelimCutRatio
	} else{
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
		warning(paste0("The spectra count of the substance ", msmsPeaks$id, " (", nrow(msmsPeaks$childHeaders), " spectra) doesn't match the provided spectra list (", 
						length(spectraList), " spectra).")
	)
    return(list(specOK=FALSE))
    
  }
  
  # here follows the Rcdk analysis
  #------------------------------------
  parentPeaks <- data.frame(mzFound=msmsPeaks$mz$mzCenter, 
                            formula=msmsPeaks$formula,
                            dppm=0,
                            x1=0,x2=0,x3=0)
  
	# get the adduct additions
	adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
	allowed_additions <- adductProperties$addition
	mode.charge <- adductProperties$charge
  
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
  
  parentAtomList <- formulastring.to.list(parent_formula)
  
  rAtoms <- which(c(!grepl("Br",parent_formula),!grepl("N",parent_formula), !grepl("S",parent_formula),!grepl("Cl",parent_formula)))
  
  if(length(rAtoms)){
    newfragData <- fragData[-which(apply(occurrenceMatrix[,rAtoms,drop=FALSE], 1, any)),]
    newfragDataIndex <- fragDataIndexing(newfragData)
  } else{
    newfragData <- fragData
    newfragDataIndex <- fragDataIndex
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
    shot_lo <- shot[(shot$int < cut) | (shot$int < max(shot$int) * cut_ratio),]
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
    
    
    # check whether formula can be fitted into precalculated fragment data
    usePrecalcData <- is.sub.formula(parent_formula,list("C"=50,"H"=70,"N"=50,"O"=50,"S"=10,"Cl"=10,"Br"=10))
    
    
    if(usePrecalcData){
      
		# if yes, split the peaks into 2 sets: Those which can be identified using the
		# fragment data, and those which can't (by maximum fragment Data m/z - 0.5)

		precalcIndex <- which(shot[,mzColname] < (max(newfragData$mass)-0.5))
		smallShots <- shot[precalcIndex,]
		bigShots <- shot[-precalcIndex,]

		# for smaller peaks use the precalculated fragment data
		if(nrow(smallShots)){
		peakmatrixSmall <- lapply(smallShots[,mzColname],function(mass){
			mass.calc <- mass + mode.charge * .emass
			maxminMass <- ppm(mass.calc, ppmlimit, l=TRUE)
			# retrieve only possibly relevant rows of the fragment Data (use the index)
			indices <- newfragDataIndex[floor(maxminMass[2]*100)]:newfragDataIndex[ceiling(maxminMass[1]*100)]
			
			fragmentedFragmentData <- newfragData[indices,]
			
			# narrow it down further using the ppm
			fragmentedFragmentData <- fragmentedFragmentData[which(fragmentedFragmentData$mass > maxminMass[2] & fragmentedFragmentData $mass < maxminMass[1]),]
			
			# return nothing if the narrowed down data is empty at this point
			if(!nrow(fragmentedFragmentData)){
				return(t(c(mzFound=as.numeric(as.character(mass)),
						formula=NA, mzCalc=NA)))
			}
			
			
			# find out if the narrowed down fragments have a sub-formula of the parentformula
			# return the indexes of relevant formulas
			fragmentedFragmentData <- fragmentedFragmentData[which(sapply(fragmentedFragmentData$formula, function(currentFormula) is.sub.formula(currentFormula,parentAtomList))),]
			
			# return nothing if the narrowed down data is empty at this point
			if(!nrow(fragmentedFragmentData)){
			return(t(c(mzFound=as.numeric(as.character(mass)),
						formula=NA, mzCalc=NA)))
			}
			
			return(t(sapply(1:nrow(fragmentedFragmentData), function(f)
			{
			mzCalc <- fragmentedFragmentData$mass[f] - mode.charge * .emass 
			c(mzFound=as.numeric(as.character(mass)),
				formula=fragmentedFragmentData$formula[f], 
				mzCalc=mzCalc)
			})))
			
        })
      } else{
        peakmatrixSmall <- list()
      }
      
      # for bigger peaks use generate.formula from rcdk
      if(nrow(bigShots)){
        system.time(peakmatrixBig <- lapply(bigShots[,mzColname], function(mass){
          # Circumvent bug in rcdk: correct the mass for the charge first, then calculate uncharged formulae
          # finally back-correct calculated masses for the charge
          mass.calc <- mass + mode.charge * .emass
          peakformula <- suppressWarnings(generate.formula(mass.calc, ppm(mass.calc, ppmlimit, p=TRUE), 
                                                   limits, charge=0))
          if(length(peakformula)==0){
            return(t(c(mzFound=as.numeric(as.character(mass)),
                       formula=NA, mzCalc=NA)))
          }else{
            return(t(sapply(peakformula, function(f)
            {
              mzCalc <- f@mass - mode.charge * .emass 
              c(mzFound=as.numeric(as.character(mass)),
                formula=f@string, 
                mzCalc=mzCalc)
            })))
          }
        }))
      } else{
        peakmatrixBig <- list()
      }
      
      peakmatrix <- c(peakmatrixSmall,peakmatrixBig)
      
    } else{
			peakmatrix <- lapply(shot[,mzColname], function(mass){
			# Circumvent bug in rcdk: correct the mass for the charge first, then calculate uncharged formulae
			# finally back-correct calculated masses for the charge
			mass.calc <- mass + mode.charge * .emass
			peakformula <- suppressWarnings(generate.formula(mass.calc, ppm(mass.calc, ppmlimit, p=TRUE), 
												   limits, charge=0))
			if(!is.list(peakformula) || length(peakformula)==0){
				return(t(c(mzFound=as.numeric(as.character(mass)),
						formula=NA, mzCalc=NA)))
			}else{
				return(t(sapply(peakformula, function(f)
				{
					mzCalc <- f@mass - mode.charge * .emass 
					c(mzFound=mass,
					formula=f@string, 
					mzCalc=mzCalc)
				})))
			}
        })
    }
    
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
    bestPpm <- aggregate(as.data.frame(childPeaksInt$dppm), list(childPeaksInt$mzFound),
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
  
  mzranges <- t(sapply(shots, function(p) {
    if(!is.null(p$childRaw)){
      return(range(p$childRaw[,mzColname]))
    } else {
      return(c(NA,NA))
    }
  }))
  
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

findPeaksInFracData <- function(mass, fracData, fracDataIndex){
  
}
