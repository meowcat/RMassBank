# Based on the code from MSnbase
# * contributed by Guangchuang Yu <guangchuangyu@gmail.com>
# * Modified by Sebastian Gibb <mail@sebastiangibb.de>
# * adapted into RMassBank by Michele Stravs <stravsmi@eawag.ch>
## setMethod("writeMgfData",
##           signature = signature("Spectrum1"),
##           function(object,
##                    con = "spectrum.mgf",
##                    COM = NULL,
##                    TITLE = NULL) {
##             writeMgfDataFile1(list(object), con = con, COM = COM, TITLE = TITLE,
##                              verbose = FALSE)
##           })
## 
  
setMethod("writeMgfData",
		signature=signature("RmbSpectraSet"),
		function(object,
				con="spectrum.mgf",
				COM = NULL,
				TITLE = NULL,
				exactPrecursor=FALSE) {
			writeMgfSpectraSet(object,
					con = con, COM = COM, TITLE = TITLE,
					verbose = FALSE, exactPrecursor=exactPrecursor)
		})

setMethod("writeMgfData",
		signature=signature("RmbSpectrum2List"),
		function(object,
				con="spectrum.mgf",
				COM = NULL,
				TITLE = NULL) {
			writeMgfRmbSpectrum2List(object,
					con = con, COM = COM, TITLE = TITLE,
					verbose = FALSE)
		})


writeMgfSpectraSet <- function(object,
					con = con, COM = COM, TITLE = TITLE,
					verbose = FALSE, exactPrecursor = FALSE)
{
	if (class(con) == "character" && file.exists(con)) {
		log_info("Overwriting ", con, "!")
		unlink(con)
	}
	
	con <- file(description = con, open = "at")
	on.exit(close(con))
	
	if (is.null(COM)) {
		COM <- paste0("COM=", object@id, " ", object@name, " ", object@mz, " ", object@formula,
				"exported by MSnbase/RMassBank on ", date())
	}
	cat(COM, file = con, sep = "")
	MSnbase:::writeMgfContent(object@parent, TITLE = NULL, con = con)
	for(chi in as.list(object@children))
	{
		if(exactPrecursor)
			chi@precursorMz <- object@mz
		MSnbase:::writeMgfContent(chi, TITLE=NULL, con=con)
	}	
}

writeMgfRmbSpectrum2List <- function(object,
		con = con, COM = COM, TITLE = TITLE,
		verbose = FALSE)
{
	if (class(con) == "character" && file.exists(con)) {
		log_info("Overwriting ", con, "!")
		unlink(con)
	}
	
	con <- file(description = con, open = "at")
	on.exit(close(con))
	
	if (is.null(COM)) {
		COM <- paste0("COM=", "RmbSpectrum2List",
				" exported by MSnbase/RMassBank on ", date())
	}
	cat(COM, file = con, sep = "")
	for(chi in as.list(object))
	{
		MSnbase:::writeMgfContent(chi, TITLE=NULL, con=con)
	}	
}

