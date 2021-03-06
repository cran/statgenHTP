% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createTimePoints.R
\name{summary.TP}
\alias{summary.TP}
\title{Summary function for TP objects}
\usage{
\method{summary}{TP}(object, ...)
}
\arguments{
\item{object}{An object of class TP.}

\item{...}{Ignored.}
}
\value{
No return value, a summary is printed.
}
\description{
Function for creating a short summary of the contents of a \code{TP} object.
The summary consists of the name of the experiment, the number of time
points, the first and last time point and the genotypes defined as checks.
}
\examples{
## Create a TP object containing the data from the Phenovator.
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Replicate",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("check1", "check2",
                                               "check3","check4"))
## Create a summary.
summary(phenoTP)

}
\seealso{
Other functions for data preparation: 
\code{\link{as.data.frame.TP}()},
\code{\link{createTimePoints}()},
\code{\link{getTimePoints}()},
\code{\link{plot.TP}()},
\code{\link{removeTimePoints}()}
}
\concept{functions for data preparation}
