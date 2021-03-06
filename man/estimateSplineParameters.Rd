% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSpline.R
\name{estimateSplineParameters}
\alias{estimateSplineParameters}
\title{Extract estimates from fitted splines.}
\usage{
estimateSplineParameters(
  HTPSpline,
  estimate = c("predictions", "derivatives", "derivatives2"),
  what = c("min", "max", "mean", "AUC", "p"),
  AUCScale = c("min", "hour", "day"),
  timeMin = NULL,
  timeMax = NULL,
  genotypes = NULL,
  plotIds = NULL
)
}
\arguments{
\item{HTPSpline}{An object of class HTPSpline, the output of the
\code{\link{fitSpline}} function.}

\item{estimate}{The P-Spline component for which the estimate should be
extracted, the predictions, the first derivatives or the second derivatives
("derivatives2")}

\item{what}{The type of estimate that should be extracted. Either minimum
("min"), maximum ("max"), mean, area under the curve ("AUC") or a percentile.
Percentiles should be given as p + percentile. E.g. for the 10th percentile
specify what = "p10"}

\item{AUCScale}{The area under the curve is dependent on the scale used on
the x-axis. By default the area is computed assuming a scale in minutes. This
can be changed to either hours or days.}

\item{timeMin}{The lower bound of the time interval from which the
estimates should be extracted. If \code{NULL} the smallest time value for
which the splines were fitted is used.}

\item{timeMax}{The upper bound of the time interval from which the
estimates should be extracted. If \code{NULL} the largest time value for
which the splines were fitted is used.}

\item{genotypes}{A character vector indicating the genotypes for which
estimates should be extracted. If \code{NULL}, estimates will be extracted
for all genotypes for which splines where fitted.}

\item{plotIds}{A character vector indicating the plotIds for which
estimates should be extracted. If \code{NULL}, estimates will be extracted
for all plotIds for which splines where fitted.}
}
\value{
A data.frame containing the estimated parameters.
}
\description{
Function for extracting parameter estimates from fitted splines on a
specified interval.
}
\examples{
## Run the function to fit P-splines on a subset of genotypes.
subGeno <- c("G160", "G151")
fit.spline <- fitSpline(inDat = spatCorrectedVator,
                        trait = "EffpsII_corr",
                        genotypes = subGeno,
                        knots = 50)

## Estimate the maximum value of the predictions at the beginning of the time course.
paramVator <- estimateSplineParameters(HTPSpline = fit.spline,
                                       estimate = "predictions",
                                       what = "max",
                                       timeMin = 1527784620,
                                       timeMax = 1528500000,
                                       genotypes = subGeno)
head(paramVator)

}
\seealso{
Other functions for fitting splines: 
\code{\link{fitSpline}()},
\code{\link{plot.HTPSpline}()}
}
\concept{functions for fitting splines}
