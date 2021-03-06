% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/webAccess.R
\name{getCompTox}
\alias{getCompTox}
\title{Retrieves DTXSID (if it exists) from EPA Comptox Dashboard}
\usage{
getCompTox(query)
}
\arguments{
\item{query}{The InChIKey of the compound.}
}
\value{
Returns the DTXSID.
}
\description{
Retrieves DTXSID (if it exists) from EPA Comptox Dashboard
}
\examples{

\dontrun{
# getCompTox("MKXZASYAUGDDCJ-NJAFHUGGSA-N")
}

}
\author{
Adelene Lai <adelene.lai@uni.lu>
}
