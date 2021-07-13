#' Checks for isotopes in a \code{msmsWorkspace}
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param intensity_cutoff The cutoff (as an absolute intensity value) under which isotopic peaks shouldn't be checked for or accepted as valid.
#'        Please note: The cutoff is not hard in the sense that it interacts with the intensity_precision parameter.
#' @param intensity_precision The difference that is accepted between the calculated and observed intensity of a possible isotopic peak. Further details down below.
#' @param conflict Either "isotopic"(Peak formulas are always chosen if they fit the requirements for an isotopic peak)
#' 		  or "strict"(Peaks are only marked as isotopic when there hasn't been a formula assigned before.)
#' @param isolationWindow Half of the width of the isolation window in Da 
#' @param evalMode Currently no function yet, but planned. Currently must be "complete"
#' @param plotSpectrum A boolean specifiying whether the spectrumshould be plotted
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @details text describing parameter inputs in more detail.
#' \itemize{
#'  \item{\code{intensity_precision}}{This parameter determines how strict the intensity values should adhere to the calculated intensity in relation to the parent peak.
#'  Options for this parameter are \code{"none"}, where the intensity is irrelevant, \code{"low"}, which has an error margin of 70\% and \code{"high"}, where the error margin 
#'  is set to 35\%. The recommended setting is \code{"low"}, but can be changed to adjust to the intensity precision of the mass spectrometer.}
#  \item{\code{evalMode}}{This parameter sets what should be done after the isotopic check. The option "add" adds failpeaks if they are isotopic peaks of previously matched peaks.
#  "check" checks matched peaks with formulas for isotopes and removes them if no isotopic peaks have been found for all formulas. The formula is also
#  adjusted, if the one with matching isotopes doesn't have the lowest dppm. "complete" does both.}
#' }
#' @return The \code{msmsWorkspace} with annotated isolation peaks
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
checkIsotopes <- function(w, mode = "pH", intensity_cutoff = 0, intensity_precision = "none", conflict = "strict", 
							isolationWindow = 2, evalMode = "complete", plotSpectrum = TRUE, settings = getOption("RMassBank")){

	# Load library and data
	requireNamespace("enviPat",quietly=TRUE)
	
	data("isotopes", package="enviPat")
	
	if(!(intensity_precision %in% c("none","low","high"))){
		stop('intensity_precision must be specified as either "none", "low" or "high"')
	}
	
	switch(intensity_precision, 
		none={
			precisionVal <- Inf
		},
		low={
			precisionVal <- 0.7
		},
		high={
			precisionVal <- 0.35
		}
	)
	
	# Load filtersettings
	filterSettings = settings$filterSettings
	
	# get the adduct additions
	adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
	allowed_additions <- adductProperties$addition
	mode.charge <- adductProperties$charge
	
	# "default" isotopes (i.e. those with the highest abundance)
	defIsotopes <- c("107Ag", "27Al", "40Ar", "75As", "197Au", "11B", "138Ba", "9Be", "209Bi",
	"79Br", "12C", "40Ca", "114Cd", "140Ce", "35Cl", "59Co", "52Cr", "133Cs", "63Cu", "164Dy",
	"166Er", "153Eu", "19F", "56Fe", "69Ga", "158Gd", "74Ge", "1H", "4He", "180Hf", "202Hg",
	"165Ho", "127I", "115In", "193Ir", "39K", "84Kr", "139La", "7Li", "175Lu", "24Mg", "55Mn",
	"98Mo", "14N", "23Na", "93Nb", "142Nd", "20Ne", "58Ni", "16O", "192Os", "31P", "231Pa",
	"208Pb", "106Pd", "141Pr", "195Pt", "85Rb", "187Re", "103Rh", "102Ru", "32S", "121Sb",
	"45Sc", "80Se", "28Si", "152Sm", "120Sn", "86Sr", "88Sr", "181Ta", "159Tb", "130Te",
	"232Th", "48Ti", "205Tl", "169Tm", "238U", "51V", "184W", "132Xe", "89Y", "174Yb",
	"64Zn", "90Zr")

	
	# Get the ppm limit from the settings
	ppmlimit <- filterSettings$ppmFine
	
	# Get the cutoff from the settings
	cut <- filterSettings$fineCut
	cut_ratio <- filterSettings$fineCutRatio

	# Extract matched and unmatched peaks from the aggregated peaks
	matchedPeaks <- peaksMatched(w)
	unmatchedPeaks <- peaksUnmatched(w)
	
	wEnv <- environment()
	
	# lapply over all runs
	lapply(w@spectra, function(spec){
		
		# Find parent formula and cpdID
		parent_formula <- add.formula(spec@formula, allowed_additions)
		id <- as.numeric(spec@id)
		specNum <- 0
		specEnv <- environment()
		# lapply over all extracted MS2 spectra
		lapply(spec@children, function(msmsdata){
			
			specEnv$specNum <- specEnv$specNum + 1
			# Extract currently relevant peaks
			currentMPeaks <- matchedPeaks[(matchedPeaks$cpdID == id) & (matchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			currentUPeaks <- unmatchedPeaks[(unmatchedPeaks$cpdID == id) & (unmatchedPeaks$scan == msmsdata@acquisitionNum),,drop=FALSE]
			
			if(nrow(currentMPeaks)){
				rownames(currentMPeaks) <- 1:nrow(currentMPeaks)
			} else {
				log_info(paste0("Compound ", id, " in spectrum #", specEnv$specNum," does not have matched peaks, so no isotopes can be found"))
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
				}
				return(0)
			}
			
			if(nrow(currentUPeaks)){
				rownames(currentUPeaks) <- 1:nrow(currentUPeaks)
			}
			
			# Generate isotopic patterns of the matched peaks
			# sort out possible isotopic peaks according to
			# isolationwindow and intensity
			isoMPatterns <- lapply(1:nrow(currentMPeaks), function(index){
                                .findPattern(currentMPeaks[index,,drop=FALSE], defIsotopes = defIsotopes, intensity_cutoff = intensity_cutoff, 
									precisionVal = precisionVal, ppmlimit = ppmlimit, isolationWindow = isolationWindow)
                            })
			
			# Name the isopatterns
			names(isoMPatterns) <- currentMPeaks$formula
			
			
			# Which isotope patterns still have theoretical intensities above the cutoff?
			peaksToCheck <- which(as.logical(sapply(isoMPatterns,nrow)))
			
            # If there are no peaks left, then abort for this spectrum
			if(!length(peaksToCheck)){
				log_info(paste0("The already annotated peaks of compound ", id, " in spectrum #", specEnv$specNum," are not intense enough to search for isotopic peaks"))
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
				}
				return(0)
			}
			
			# Now, look for isotopic patterns in unmatched peaks with all specified parameters

			# Which peaks have no formula annotation within the dppm as of now?
			peaksNoAnnotation <- currentUPeaks[which(is.na(currentUPeaks$dppm) | abs(currentUPeaks$dppmBest) > ppmlimit),]
			
            # What is the mean of the dppm of the currently "OK" peaks?
            # (Used for calculating the score parameter for intensities)
            dppmMean <- mean(abs(currentMPeaks[currentMPeaks$filterOK,]$dppm))
            
			# If there are any peaks without annotation:
			if(nrow(peaksNoAnnotation)){
				UPList <- .findMatchingPeaks(peaksNoAnnotation,isoMPatterns[peaksToCheck], dppmMean)
			} else{
				# If there are no peaks, fake an empty dataset
				UPList <- list(matrix(character(0),0,29))
			}
            
			# If conflict is "strict", end here (And plot, maybe)
			if(conflict == "strict"){
                # Generate matrix of peaks that should be added
                additionMatrix <- .peakReasoner(list(matrix(character(0),0,29)) , UPList, currentMPeaks)
				if(plotSpectrum){
					plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
					if(nrow(additionMatrix)){
						points(additionMatrix$mzFound, additionMatrix$intensity,type="h", col="green", lwd=3)
					}
				}
                # If there is something in the matrix
                if(nrow(additionMatrix)){
                    # Add all the peaks that could be annotated as isotopic
                    wEnv$w@aggregated[additionMatrix$index,] <- additionMatrix
                }
				return(0)
			}
			
			
			# Now check the matched peaks for the patterns and put all peaks in a matrix
			MPList <- .findMatchingPeaks(currentMPeaks[currentMPeaks$filterOK,], isoMPatterns[peaksToCheck], dppmMean)
			
			# Generate matrix of peaks that should be corrected
			correctionMatrix <- .peakReasoner(MPList, UPList, currentMPeaks)
            
            
            # If there is something in the matrix
            if(nrow(correctionMatrix)){
                # Add all the peaks that could be annotated as isotopic
                wEnv$w@aggregated[correctionMatrix$index,] <- correctionMatrix
            }
			
            # If the newly annotated peaks should be plotted, plot them
            if(plotSpectrum){
                plot(currentMPeaks$mzFound, currentMPeaks$intensity,type="h", main=paste(id,findName(id)), col="black", xlab="m/z", ylab="intensity", lwd=3)
                if(nrow(correctionMatrix)){
                    points(correctionMatrix$mzFound, correctionMatrix$intensity,type="h", col="green", lwd=3)
                }
            }
            return(0)
		})
	})
	return(w)
}

# Pattern finding and evaluation of these patterns inside the checkIsotopes function harmed readability and complicated debugging
# So modularize this function
.findPattern <- function(aggregateRow, defIsotopes, intensity_cutoff, precisionVal, ppmlimit, isolationWindow){
	
	# Find pattern and mass
	outp <- capture.output(pattern <- as.data.frame(enviPat::isopattern(aggregateRow[,"formula"], isotopes = isotopes)[[1]]))
    rm("outp")
	mass <- as.numeric(aggregateRow[,"mzCalc"])
	
	# Find index of monoisotopic molecule and
	# normalize abundance that monoisotopic molecule has always "1"
	mainIndex <- which.min(abs(pattern[,"m/z"]-mass))
	pattern[,"abundance"] <- round(pattern[,"abundance"] * (100/pattern[mainIndex,"abundance"]),3)
	
    # Add the formula of every isotope to every row in the pattern
	pattern <- .annotateFormulaToEnviPatTable(pattern, defIsotopes)
	
	# Sort pattern by abundance
	pattern <- pattern[order(pattern[,"abundance"],decreasing = T),][-mainIndex,]
	
	# Look up which patterns have their m/z inside the isolation window +- 0.5
	# and delete all others 
	keepMzIndex <- which(pattern[,"m/z"] < mass + isolationWindow + 0.2 & pattern[,"m/z"] > mass - isolationWindow - 0.2)
	pattern <- pattern[keepMzIndex,,drop=FALSE]
	if(nrow(pattern) == 0){
		return(pattern)
	}
	
	# Calculate the expected intensities according to the abundance
	# See which expected isotope peaks have an intensity higher than the cutoff
	
	# Find the absolute intensity ranges
	intensities <- apply(pattern, 1, function(patternRow){
		as.numeric(aggregateRow[,"intensity"]) * as.numeric(patternRow["abundance"])/100
	})
	
	roundedInt <- round(intensities,digits=-2)
	keepIntIndex <- which((roundedInt + roundedInt * precisionVal) >= intensity_cutoff)
	pattern <- pattern[keepIntIndex,,drop=FALSE]
	if(nrow(pattern)){
		pattern$minintensity <- roundedInt[keepIntIndex] - roundedInt[keepIntIndex] * precisionVal
		pattern$maxintensity <- roundedInt[keepIntIndex] + roundedInt[keepIntIndex] * precisionVal
        pattern$expectedintensity <- roundedInt[keepIntIndex]
		pattern$monoisotopicFormula <- aggregateRow[,"formula"]
        pattern$monoisotopicIndex <- aggregateRow[,"index"]
        pattern$monoisotopicMz <- aggregateRow[,"mzFound"]
		# Calculate the expected mz range of the isotope peak
		mzCols <- t(sapply(pattern[,"m/z"], ppm, dppm=ppmlimit,l=T))
		colnames(mzCols) <- c("maxMZ","minMZ")
		pattern <- cbind(pattern,mzCols)
	}
	
	return(pattern)
}

.findMatchingPeaks <- function(aggregatedPeakMatrix, isoPatterns, dppmMean){
	# Only unique mzs (peaks can be annotated multiple times and will therefore be included multiple times)
	# It doesn't matter which one, so take the first
	aggregatedPeakMatrix <- aggregatedPeakMatrix[sapply(unique(aggregatedPeakMatrix$mzFound), function(mz){
		which(aggregatedPeakMatrix$mzFound == mz)[1]
	}),]
	
	# Iterate over all supplied patterns
	lapply(isoPatterns, function(pattern){

		do.call(rbind,lapply(1:nrow(pattern), function(index){
            patternRow <- pattern[index,,drop=FALSE]
			# Find peaks that fit the specified intensities and m/z
			pIndex <- which(aggregatedPeakMatrix[,"mzFound"] < patternRow[,"maxMZ"]
									& aggregatedPeakMatrix[,"mzFound"] > patternRow[,"minMZ"]
									& aggregatedPeakMatrix[,"intensity"] < patternRow[,"maxintensity"]
									& aggregatedPeakMatrix[,"intensity"] > patternRow[,"minintensity"]
			)
			
			# Note these Peaks
			candidates <- aggregatedPeakMatrix[pIndex,]
			# If there are any: Change parameters in the aggregated matrix
			if(nrow(candidates)){
                # General parameters that need to be changed
				candidates$dppm <- (candidates$mzFound /  as.numeric(patternRow[,"m/z"]) - 1) * 1e6
				candidates$mzCalc <- patternRow[,"m/z"]
				candidates$formula <- patternRow[,"formula"]
				candidates$matchedReanalysis <- NA
				candidates$filterOK <- TRUE
				candidates$good <- TRUE
				candidates$dppmBest <- candidates$dppm
				candidates$formulaCount <- 1
                # New parameters (are used to make peak reasoning easier, we would
                # need to find them out later anyways, so do it now when it's easy)
                
                # 1) The monoisotopic formula of the now added peak
                candidates$monoisotopicFormula <- patternRow[,"monoisotopicFormula"]
                
                # 2) The calculated expected intensity (added for debug reasons)
                candidates$intCalc <- patternRow[,"expectedintensity"]
                
                # 3) Produce "scores" in the range of dppm and scales nicely
                #candidates$dint <- (max(as.numeric(candidates$intensity),as.numeric(patternRow["expectedintensity"]))/min(as.numeric(candidates$intensity),as.numeric(patternRow["expectedintensity"])))^2
                maxInt <- max(as.numeric(candidates$intensity),as.numeric(patternRow[,"expectedintensity"]))
                minInt <- min(as.numeric(candidates$intensity),as.numeric(patternRow[,"expectedintensity"]))

                candidates$dint  <- (maxInt/minInt) * dppmMean
                candidates$monoisotopicIndex <- patternRow[,"monoisotopicIndex"]
                
                candidates$monoisotopicMz <- as.numeric(patternRow[,"monoisotopicMz"])
			} else{
                candidates$monoisotopicFormula <- character(0)
                candidates$intCalc <- numeric(0)
                candidates$dint <- numeric(0)
                candidates$monoisotopicIndex <- integer(0)
                candidates$monoisotopicMz <- numeric(0)
            }
			return(candidates)
		}))
	})
}

.peakReasoner <- function(MPList, UPList, currentMPeaks){
	
    # Unlist the possible isotope peak lists for the previously matched and previously unmatched peaks
	adjustmentMatrix <- rbind(do.call(rbind, MPList),do.call(rbind, UPList))
	
    # Remove peaks with a much too high dint
    dintRemoval <- which(adjustmentMatrix$dint>50)
    if(as.logical(length(dintRemoval))){
        adjustmentMatrix <-  adjustmentMatrix[-which(adjustmentMatrix$dint>50),]
    }
    # No isotopes found, abort
    if(!as.logical(nrow(adjustmentMatrix))){
		return(adjustmentMatrix)
	}
    
    # Now we need to remove rows from this matrix, so that 4 conditions are fulfilled:
    # a) Values in the mzFound can not occur in monoisotopicMz (because an isotope can't be monoisotopic)
    # b) No monoisotopicMz can have more than one index assigned (because this would mean a peak has 2 formulas)
    # c) All unique mzFound can only be included once
    # d) The number of indices must be minimal (it is highly likely that the monoisotopic peak with the most 
    # isotope peaks is correctly annotated)

    uniqueMonoMz <- unique(adjustmentMatrix$monoisotopicMz)
    uniqueMzFound <- unique(adjustmentMatrix$mzFound)
    
    checkAllMzFound <- function(){
        all(uniqueMzFound %in% adjustmentMatrix$mzFound) && all(table(adjustmentMatrix$mzFound) == 1)
    }
    
    # Condition a) Go through all monoisotopicMz
    # and remove those that happen to be in mzFound
    removalIndex <- which(adjustmentMatrix$monoisotopicMz %in% uniqueMzFound)
    if(length(removalIndex)){
        adjustmentMatrix <- adjustmentMatrix[-removalIndex,]
    }
    
    
    # Condition b) and d)
    for(monoMz in uniqueMonoMz){
        monoMzRows <- which(adjustmentMatrix$monoisotopicMz == monoMz)
        
        # If there are different indices, use the one that appears most often, else there is no problem
        if(!all(adjustmentMatrix$monoisotopicIndex[monoMzRows] == adjustmentMatrix$monoisotopicIndex[monoMzRows][1])){
            # Store all indices that occur equally as often and most often
            maxIndex <- vector()
            maxlength <- 0
            for(index in unique(adjustmentMatrix$monoisotopicIndex[monoMzRows])){
                currlength <- length(which(adjustmentMatrix$monoisotopicIndex[monoMzRows] == index))
                if(currlength > maxlength){
                    maxlength <- currlength
                    maxIndex <- index
                }
                if(currlength == maxlength){
                    maxlength <- currlength
                    if(!(index %in% maxIndex)){
                        maxIndex <-  c(maxIndex,index)
                    }
                }
            }
            
            # One index appears most often, so use that one
            if(length(maxIndex) == 1){
                monoMzRows <- monoMzRows[-which(adjustmentMatrix$monoisotopicIndex[monoMzRows] == maxIndex)]
                adjustmentMatrix <- adjustmentMatrix[-monoMzRows,]
            } else { # Sometimes different indices appear equally as often and most often, use scoring metric
                bestscore <- 1000
                bestIndex <- 0
                # Go through all possible indices
                for(conflictIndex in maxIndex){
                    # Calculate Score and note the best score
                    score <- mean(apply(abs(adjustmentMatrix[which(adjustmentMatrix$monoisotopicIndex == conflictIndex),c("dint","dppm")]),1,sum))
                    if(score < bestscore){
                        bestscore <- score
                        bestIndex <- conflictIndex
                    }
                }
                # Throw out all indices that don't have the best average score
                monoMzRows <- monoMzRows[-which(adjustmentMatrix$monoisotopicIndex[monoMzRows] == bestIndex)]
                adjustmentMatrix <- adjustmentMatrix[-monoMzRows,]
            }
        }
    }
    
    # Everything should be fine, but if it is not
    # some peaks have the same monoisotopicMz but not the same Formula,
    # i.e. one monoisotopic peak has 2 isotope formulas that fit the
    # same peak
    if(length(unique(adjustmentMatrix$mzFound)) != length(adjustmentMatrix$mzFound)){
        occurrences <- as.data.frame(table(adjustmentMatrix$mzFound))
        for(value in occurrences[which(occurrences[,2] > 1),1]){
            rows <- which(adjustmentMatrix$mzFound == value)
            scores <- apply(adjustmentMatrix[rows,c("dint","dppm")],1,sum)
            notMinRows <- rows[-which.min(scores)]
            adjustmentMatrix <- adjustmentMatrix[-notMinRows,]
        }
    }
    
    # Now: Check for every peak if the monoisotopic peak is already the peak that is passed
    # else change the peak that gets passed
    problemPeaks <- which(!sapply(adjustmentMatrix$monoisotopicIndex, function(index) currentMPeaks$filterOK[which(currentMPeaks$index == index)]))
    problemMz <- adjustmentMatrix$monoisotopicMz[problemPeaks]
    correctIndex <- adjustmentMatrix$monoisotopicIndex[problemPeaks]
    
    adjustmentMatrix <- adjustmentMatrix[,-(25:29)]
    
    if(as.logical(length(problemMz))){
        # Go through every mz and change the peak to the correct one
        for(p in 1:length(problemMz)){
            addPeak <- currentMPeaks[which(currentMPeaks$index == correctIndex[p]),,drop=FALSE]
            addPeak$filterOK <- TRUE
            removePeak <- currentMPeaks[which((currentMPeaks$mzFound == problemMz[p]) & currentMPeaks$filterOK),,drop=FALSE]
            removePeak$filterOK <- FALSE
            adjustmentMatrix <- rbind(adjustmentMatrix,addPeak,removePeak)
        }
    }
    
    return(adjustmentMatrix)
}


# Further modularization of formula addition
# Is essentially one task, so this makes it more understandable
.annotateFormulaToEnviPatTable <- function(pattern, defIsotopes){
    # Find all isotope atom names (only in colnames of patterns, sadly)
	isoCols <- colnames(pattern)[3:ncol(pattern)]
	
	# Order the formulas to be in Hill system:
	# First the Cs, then the Hs, then lexicographical
	# First: Split isotope names so that there only are atoms
	# Find position of C and H
	splitNames <- gsub('^([0-9]+)(.+)$', '\\2', isoCols)
	CHPos <- c(which(splitNames == "C"),which(splitNames == "H"))
	
	# Account for special case: No Cs or Hs
	if(length(CHPos)){
		# If there are Cs and Hs, overwrite the positions so they always get sorted to the front
		splitNames[CHPos] <- ""
	}
	
	# Default isotopes are always first in the colnames, so no need to order internally between isotopes
	atomOrder <- order(splitNames)
	
	# If there are Cs and Hs, the order must be preserved
	if(length(CHPos)){
		atomOrder[1:length(CHPos)] <- CHPos
	}
	# else it is already ordered
	
	# Mark new names for formula creation.
	newnames <- unname(sapply(isoCols, function(currentCol){
		if(currentCol %in% defIsotopes){
			# "Default" isotope (no isotope mass in formula, e.g. "C")
			return(gsub('^(.{0})([0-9]+)(.+)$', '\\3', currentCol))
		}else{
			# Other isotopes (isotope mass in formula, e.g. "[13]C")
			return(gsub('^(.{0})([0-9]+)(.+)$', '\\1[\\2]\\3', currentCol))
		}
	}))
	
	# Generate the formula for every pattern
	pattern$formula <- apply(pattern,1,function(p){
		paste0(sapply(atomOrder+2,function(isoIndex){
			if(p[isoIndex] == 0){ 
				return("")
			}
			if(p[isoIndex] == 1){
				return(newnames[isoIndex-2])
			}
			return(paste0(newnames[isoIndex-2],p[isoIndex]))
		}),collapse="")
	})
    
    return(pattern)
}
