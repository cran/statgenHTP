#' Detect outliers for series of observations
#'
#' Function for detecting strange time courses. The function uses the estimates
#' for the spline coefficients per time course (typically per plant).
#' Correlations between those coefficient vectors are calculated to identify
#' outlying time courses, i.e., plants. An outlying time course will have low
#' correlation to the majority of time courses. To support the analysis by
#' correlations, a principal component analysis is done on the plant
#' (time course) by spline coefficient matrix. A PCA plot of the plant scores
#' will show the outlying plants.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param predDat A data.frame with predicted data from a fitted spline.
#' @param coefDat A data.frame with coefficients from a fitted spline.
#' @param trait A character string indicating the trait for which to detect
#' outliers.
#' @param genotypes A character vector indicating the genotypes for which to
#' detect outliers. If \code{NULL}, outlier detection will be done for all
#' genotypes.
#' @param geno.decomp A character vector indicating the variables to use to
#' group the genotypic variance in the model.
#' @param thrCor A numerical value used as threshold for determining outliers
#' based on correlation between plots.
#' @param thrPca A numerical value used as threshold for determining outliers
#' based on angles (in degrees) between PCA scores.
#'
#' @return An object of class \code{serieOut}, a \code{data.frame} with outlying
#' series of observations.
#'
#' @examples
#' \donttest{
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-splines on a subset of genotypes
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the data.frames with predicted values and P-Spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses.
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30)
#'
#' ## The `outVator` can be visualized for selected genotypes.
#' plot(outVator, genotypes = "G151")
#' }
#'
#' @family functions for detecting outliers for series of observations
#'
#' @importFrom utils combn
#' @export
detectSerieOut <- function(corrDat,
                           predDat,
                           coefDat,
                           trait,
                           genotypes = NULL,
                           geno.decomp = NULL,
                           thrCor = 0.9,
                           thrPca = 30) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(corrDat, "data.frame")) {
    stop("corrDat should be a data.frame.\n")
  }
  corrCols <- c("plotId", "genotype", trait, geno.decomp)
  if (!all(hasName(x = corrDat, name = corrCols))) {
    stop("corrDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!inherits(predDat, "data.frame")) {
    stop("predDat should be a data.frame.\n")
  }
  predCols <- c("plotId", "genotype", "pred.value")
  if (!all(hasName(x = predDat, name = predCols))) {
    stop("predDat should at least contain the following columns: ",
         paste(predCols, collapse = ", "))
  }
  if (!inherits(coefDat, "data.frame")) {
    stop("coefDat should be a data.frame.\n")
  }
  coefCols <- c("plotId", "genotype", "type", "obj.coefficients")
  if (!all(hasName(x = coefDat, name = coefCols))) {
    stop("coefDat should at least contain the following columns: ",
         paste(coefCols, collapse = ", "))
  }
  if (is.null(genotypes)) {
    genotypes <- unique(as.character(predDat[["genotype"]]))
  } else {
    if (!is.character(genotypes)) {
      stop("genotypes should be a character vector.\n")
    }
    if (!all(genotypes %in% predDat[["genotype"]])) {
      stop("all genotypes should be in predDat")
    }
  }
  if (!all(genotypes %in% coefDat[["genotype"]])) {
    stop("all genotypes should be in coefDat")
  }
  if (!all(genotypes %in% corrDat[["genotype"]])) {
    stop("all genotypes should be in corrDat")
  }
  if (!is.numeric(thrCor) || any(thrCor < 0) || any(thrCor > 1)) {
    stop("thrCor should be a numerical vector with values between 0 and 1.\n")
  }
  if (!is.numeric(thrPca) || any(thrPca < 0)) {
    stop("thrPca should be a numerical vector with positive values.\n")
  }
  genoPlotId <- sapply(X = genotypes, FUN = function(genotype) {
    length(unique(corrDat[corrDat[["genotype"]] == genotype, "plotId"]))
  })
  genoPlotIdLim <- names(genoPlotId[genoPlotId < 3])
  if (length(genoPlotIdLim) > 0) {
    warning("The following genotypes have less than 3 plotIds and are skipped ",
            "in the outlier detection:\n",
            paste(genoPlotIdLim, collapse = ", "), "\n", call. = FALSE)
    genotypes <- genotypes[!genotypes %in% genoPlotIdLim]
    if (length(genotypes) == 0) {
      stop("No genotypes left for performing outlier detection.\n")
    }
  }
  ## Restrict corrDat to genotypes.
  corrDat <- corrDat[corrDat[["genotype"]] %in% genotypes, ]
  ## Get number of values for geno.decomp.
  if (!is.null(geno.decomp)) {
    genoDecompLevs <- unique(corrDat[[geno.decomp]])
    nGenoDecomp <- length(genoDecompLevs)
    if (length(thrCor) == 1) {
      thrCor <- setNames(rep(thrCor, times = nGenoDecomp), genoDecompLevs)
    } else if (is.null(names(thrCor)) ||
               !all(genoDecompLevs %in% names(thrCor))) {
      stop("thrCor should be a named vector, with names matching the levels ",
           "in geno.decomp.\n")
    }
    if (length(thrPca) == 1) {
      thrPca <- setNames(rep(thrPca, times = nGenoDecomp), genoDecompLevs)
    } else if (is.null(names(thrPca)) ||
               !all(genoDecompLevs %in% names(thrPca))) {
      stop("thrPca should be a named vector, with names matching the levels ",
           "in geno.decomp.\n")
    }
  } else {
    if (length(thrCor) != 1) {
      stop("thrCor should be a vector of length 1.\n")
    }
    if (length(thrPca) != 1) {
      stop("thrPca should be a vector of length 1.\n")
    }
  }
  ## Restrict predDat and coefDat to genotypes.
  ## Get corrected and predicted data per genotype.
  ## Merge geno.decomp to coefDat.
  if (!is.null(geno.decomp)) {
    predDat <- merge(predDat, unique(corrDat[c("plotId", geno.decomp)]))
  }
  ## Restrict to corrected plots that are also in predictions.
  ## Some plots are removed while predicting the splines.
  corrDatPred <- corrDat[corrDat[["plotId"]] %in% predDat[["plotId"]], ]
  NAplants <- tapply(corrDatPred[[trait]], droplevels(corrDatPred[["plotId"]]),
                     FUN = function(x) {
                       all(is.na(x))
                     })
  corrDatPred <- corrDatPred[corrDatPred[["plotId"]] %in%
                               names(NAplants)[!NAplants], ]
  genoDats <- split(x = corrDatPred,
                    f = corrDatPred[c("genotype", geno.decomp)], drop = TRUE)
  predDat <- predDat[predDat[["plotId"]] %in%
                       names(NAplants)[!NAplants], ]
  genoPreds <- split(x = predDat,
                     f = predDat[c("genotype", geno.decomp)], drop = TRUE)
  ## Reshape the spline coefficients per plant x geno.
  ## First restrict to selected genotypes.
  coefDat <- coefDat[coefDat[["genotype"]] %in% genotypes, ]
  ## Merge geno.decomp to coefDat.
  if (!is.null(geno.decomp)) {
    coefDat <- merge(coefDat, unique(corrDatPred[c("plotId", geno.decomp)]))
  }
  plantDats <- lapply(X = split(x = coefDat,
                                f = coefDat[c("genotype", geno.decomp)],
                                drop = TRUE),
                      FUN = function(dat) {
                        plantDat <- reshape2::acast(dat,
                                                    formula = type ~ plotId,
                                                    value.var = "obj.coefficients")
                        ## Remove intercept.
                        plantDat <- plantDat[-1, ]
                        ## Remove plants with only NA.
                        NAplants <- apply(X = plantDat, MARGIN = 2,
                                          FUN = function(plant) {
                                            all(is.na(plant))
                                          })
                        plantDat <- plantDat[, !NAplants]
                        ## Add geno.decomp as attribute for later use.
                        if (!is.null(geno.decomp)) {
                          attr(plantDat, which = "genoDecomp") <-
                            as.character(unique(dat[[geno.decomp]]))
                        }
                        return(plantDat)
                      })
  ## Compute correlation matrix.
  cormats <- lapply(X = plantDats, FUN = function(plantDat) {
    if (!is.null(dim(plantDat))) {
      ## if there are plants, estimate the correlation...
      cormat <- cor(plantDat)
      diag(cormat) <- NA
      attr(cormat, which = "genoDecomp") <- attr(plantDat, which = "genoDecomp")
      if (!is.null(geno.decomp)) {
        attr(cormat, which = "thrCor") <-
          thrCor[attr(plantDat, which = "genoDecomp")]
      } else {
        attr(cormat, which = "thrCor") <- thrCor
      }
      return(cormat)
    } else {
      ## ... if not, return a null matrix
      return(NULL)
    }
  })
  plantPcas <- lapply(X = plantDats, FUN = function(plantDat) {
    ## Perform a PCA on the spline coefficients per genotype.
    ## Run the PCA.
    if (!is.null(dim(plantDat))) {
      plantPca <- prcomp(plantDat, center = TRUE, scale. = TRUE)
      attr(plantPca, which = "genoDecomp") <- attr(plantDat, which = "genoDecomp")
      ## if there are plants, perform the PCA...
      return(plantPca)
    } else {
      ## ... if not, return a null object
      return(NULL)
    }
  })
  annotatePlantsCor <- lapply(X = names(cormats), FUN = function(geno) {
    if (!is.null(cormats[[geno]])) {
      ## Compute mean based on all but the worst correlation per plant.
      ## This prevents one bad correlation causing all plants to be outliers.
      meanCor <- apply(cormats[[geno]], MARGIN = 1, FUN = function(plant) {
        mean(plant[-which(rank(plant) == 1)], na.rm = TRUE)
      })
      thrCorPlant <- attr(cormats[[geno]], which = "thrCor")
    } else {
      meanCor <- NULL
    }
    if (any(meanCor < thrCorPlant)) {
      ## Create data.frame with info on plants with average correlation
      ## below threshold.
      annPlotsCorr <- meanCor[meanCor < thrCorPlant]
      return(data.frame(plotId = names(annPlotsCorr), reason = "mean corr",
                        value = annPlotsCorr))
    } else {
      return(NULL)
    }
  })
  cormats <- lapply(X = cormats, FUN = function(cormat) {
    ## Set lower part of cormat to NA for plotting.
    cormat[lower.tri(cormat)] <- NA
    ## Melt to format used by ggplot.
    meltedCormat <- reshape2::melt(cormat, na.rm = TRUE)
    attr(meltedCormat, which = "thrCor") <- attr(cormat, which = "thrCor")
    return(meltedCormat)
  })
  annotatePlantsPcaAngle <- lapply(X = names(plantPcas), FUN = function(geno) {
    ## Calculate the pairwise difference of coordinates on the 2nd axis and
    ## annotate plant with average diff larger than threshold.
    plantPcaPlot <- factoextra::fviz_pca_var(plantPcas[[geno]])
    PcVecs <- as.matrix(plantPcaPlot$data[, 2:3])
    PcAngles <- matrix(nrow = nrow(PcVecs), ncol = nrow(PcVecs),
                       dimnames = list(rownames(PcVecs), rownames(PcVecs)))
    PcAngles[lower.tri(PcAngles)] <-
      combn(x = rownames(PcVecs), m = 2, FUN = function(plants) {
        angle(PcVecs[plants, ])
      }, simplify = TRUE)
    PcAngles[upper.tri(PcAngles)] <- t(PcAngles)[upper.tri(PcAngles)]
    meanPcAngles <- rowMeans(PcAngles, na.rm = TRUE)
    if (!is.null(geno.decomp)) {
      thrPcaPlant <- thrPca[attr(plantPcas[[geno]], which = "genoDecomp")]
    } else {
      thrPcaPlant <- thrPca
    }
    if (any(meanPcAngles >= thrPcaPlant)) {
      ## Create data.frame with info on plants with average difference.
      ## above threshold.
      annPlotsPcAngle <- meanPcAngles[meanPcAngles >= thrPcaPlant]
      return(data.frame(plotId = names(annPlotsPcAngle), reason = "angle",
                        value = annPlotsPcAngle))
    } else {
      return(NULL)
    }
  })
  ## Create full data.frame with annotated plants.
  annotatePlants <- do.call(rbind, c(annotatePlantsCor, annotatePlantsPcaAngle))
  if (!is.null(annotatePlants)) {
    ## Merge genotype and geno.decomp to annotated plants.
    annotatePlants <- merge(unique(corrDatPred[c("genotype", geno.decomp,
                                                 "plotId")]),
                            annotatePlants)
    ## Order by genotype, geno.decomp and plotId.
    if (!is.null(geno.decomp)) {
      annOrd <- order(annotatePlants[["genotype"]],
                      annotatePlants[[geno.decomp]], annotatePlants[["plotId"]])
    } else {
      annOrd <- order(annotatePlants[["genotype"]], annotatePlants[["plotId"]])
    }
    annotatePlants <- droplevels(annotatePlants[annOrd, ])
  } else {
    annotatePlants <- data.frame()
  }
  plotInfo <- unique(corrDatPred[c("genotype", geno.decomp)])
  class(annotatePlants) <- c("serieOut", class(annotatePlants))
  attr(x = annotatePlants, which = "thrCor") <- thrCor
  attr(x = annotatePlants, which = "thrPca") <- thrPca
  attr(x = annotatePlants, which = "trait") <- trait
  attr(x = annotatePlants, which = "geno.decomp") <- geno.decomp
  attr(x = annotatePlants, which = "plotInfo") <- plotInfo
  attr(x = annotatePlants, which = "cormats") <- cormats
  attr(x = annotatePlants, which = "plantPcas") <- plantPcas
  attr(x = annotatePlants, which = "genoPreds") <- genoPreds
  attr(x = annotatePlants, which = "genoDats") <- genoDats
  return(annotatePlants)
}

#' Plot outliers for series of observations
#'
#' Plot the fitted spline, correlation matrix and PCA biplot for each of the
#' genotypes. Outlying series of observations are shown as filled dots in the
#' fitted spline plot, other observations are shown as open dots.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class \code{serieOut}.
#' @param ... Ignored.
#' @param genotypes A character vector indicating which genotypes should be
#' plotted.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint
#' in the labels on the x-axis?
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector
#' indicating the column containing the numerical time to use.
#'
#' @return A list of ggplot objects is invisibly returned.
#'
#' @examples
#' \donttest{
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-splines on a subset of genotypes
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the data.frames with predicted values and P-Spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses.
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30)
#'
#' ## The `outVator` can be visualized for selected genotypes.
#' plot(outVator, genotypes = "G151")
#' }
#'
#' @family functions for detecting outliers for series of observations
#'
#' @export
plot.serieOut <- function(x,
                          ...,
                          genotypes = NULL,
                          useTimeNumber = FALSE,
                          timeNumber = NULL,
                          title = NULL,
                          output = TRUE) {
  if (useTimeNumber && (is.null(timeNumber) || !is.character(timeNumber) ||
                        length(timeNumber) > 1)) {
    stop("timeNumber should be a character string of length 1.\n")
  }
  thrCor <- attr(x = x, which = "thrCor")
  thrPca <- attr(x = x, which = "thrPca")
  trait <- attr(x = x, which = "trait")
  geno.decomp <- attr(x = x, which = "geno.decomp")
  plotInfo <- attr(x = x, which = "plotInfo")
  if (!is.null(genotypes) && (!is.character(genotypes) ||
                              !all(genotypes %in% plotInfo[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes used for ",
         "outlier detection.\n")
  }
  if (is.null(genotypes)) {
    genotypes <- interaction(plotInfo[c("genotype", geno.decomp)])
  } else {
    genotypes <- interaction(plotInfo[plotInfo[["genotype"]] %in% genotypes,
                                      c("genotype", geno.decomp)])
  }
  genotypes <- as.character(genotypes)
  cormats <- attr(x = x, which = "cormats")[genotypes]
  plantPcas <- attr(x = x, which = "plantPcas")[genotypes]
  genoPreds <- attr(x = x, which = "genoPreds")[genotypes]
  genoDats <- attr(x = x, which = "genoDats")[genotypes]
  ## Get minimum correlation. Cannot be higher than 0.8.
  minCor <- min(c(unlist(cormats, use.names = FALSE), 0.8), na.rm = TRUE)
  ## Compute the number of breaks for the time scale based on all plants.
  ## If there are less than 3 time points use the number of time points.
  ## Otherwise use 3.
  useTimePoint <- hasName(x = genoPreds[[1]], name = "timePoint") &&
    !useTimeNumber
  timeVar <- if (useTimePoint) "timePoint" else "timeNumber"
  timeVar2 <- if (useTimeNumber) timeNumber else timeVar
  ## Create plots.
  if (!output) {
    ## When calling arrangeGrob a blank page is opened if no plotting
    ## device is currently open.
    ## Adding a call to an empty pdf file prevents this empty plot from
    ## being opened.
    pdf(file = NULL)
    on.exit(dev.off(), add = TRUE)
  }
  p <- lapply(X = genotypes, FUN = function(genotype) {
    ## Create shape for plotting.
    ## Defaults to open circle.
    plotIds <- unique(genoPreds[[genotype]][["plotId"]])
    plotShapes <- setNames(rep(1, times = length(plotIds)), plotIds)
    ## Annotated plants get a closed circle.
    plotShapes[names(plotShapes) %in% x[["plotId"]]] <- 19
    ## Plot of time course per genotype: corrected data + spline per plant.
    kinetic <- ggplot2::ggplot(genoDats[[genotype]],
                               ggplot2::aes_string(x = timeVar2, y = trait,
                                                   color = "plotId")) +
      ggplot2::geom_point(ggplot2::aes_string(shape = "plotId"), size = 2,
                          na.rm = TRUE) +
      ggplot2::geom_line(data = genoPreds[[genotype]],
                         ggplot2::aes_string(x = timeVar,
                                             y = "pred.value"),
                         size = 0.5, na.rm = TRUE) +
      ggplot2::scale_shape_manual(values = plotShapes) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_text(size = 13))
    if (useTimePoint) {
      ## Format the time scale to Month + day.
      kinetic <- kinetic +
        ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                  labels = scales::date_format("%B %d"))
    }
    ## Correlation plot.
    thrCorGeno <- attr(cormats[[genotype]], which = "thrCor")
    correl <- ggplot2::ggplot(data = cormats[[genotype]],
                              ggplot2::aes_string("Var2", "Var1",
                                                  fill = "value")) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradientn(colors = c("red", "white", "blue"),
                                    values = scales::rescale(c(minCor,
                                                               thrCorGeno, 1)),
                                    limits = c(minCor, 1),
                                    name = "Pearson\nCorrelation") +
      ## Move y-axis to the right.
      ggplot2::scale_y_discrete(position = "right") +
      ## Use coord fixed to create a square shaped output.
      ggplot2::coord_fixed() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                     panel.border = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_line(color = "grey92"),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.ticks = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                         size = 12, hjust = 1),
                     legend.position = "left") +
      ggplot2::labs(title = "Correl of coef", x = NULL, y = NULL)
    ## PCA biplot.
    pcaplot <- factoextra::fviz_pca_var(plantPcas[[genotype]])
    ## Arrange plots.
    lay <- rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(2, 3))
    ## grid arrange always plots results.
    titleGeno <- paste("Geno", genotype,
                       if (!is.null(title)) paste("-", title), "\n")
    if (genotype == genotypes[1] && output) {
      ## Arrange grob always needs an open device.
      ## This creates a blank first page when first plotting.
      ## By opening a new page for the first plot and then using
      ## newpage = FALSE in the actual plot this blank page is overwritten.
      grid::grid.newpage()
    }
    pGeno <-
      gridExtra::arrangeGrob(kinetic, correl, pcaplot,
                             layout_matrix = lay,
                             top = grid::textGrob(label = titleGeno,
                                                  gp = grid::gpar(fontsize = 15,
                                                                  fontface = 2)))
    if (output) {
      gridExtra::grid.arrange(pGeno, newpage = (genotype != genotypes[1]))
    }
    return(pGeno)
  })
  invisible(p)
}

#' Replace outliers for series of observations by NA
#'
#' Function for replacing outliers for series of observations in the data by NA.
#' The input can either be a data.frame, specified in \code{dat}, or the output
#' of the \code{fitSpline} function, specified in \code{fitSpline}. Exactly one
#' of these should be provided as input for the function.
#'
#' @param dat A \code{data.frame}.
#' @param fitSpline An object of class \code{HTPSpline}, the output of the
#' \code{\link{fitSpline}} function.
#' @param serieOut A data.frame with at least the column plotId with
#' values corresponding to those in dat/fitSpline.
#' @param trait The trait that should be replaced by NA. Can be ignored when
#' using the output of \code{detectSerieOut} as input.
#'
#' @return Depending on the input either a \code{data.frame} or an object of
#' class \code{HTPSpline} for which the outliers specified in \code{serieOut}
#' are replace by NA.
#'
#' @examples
#' ## Run the function to fit P-spline on a subset of genotypes.
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the tables of predicted values and P-spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30)
#'
#' ## Replace the outliers by NA in the corrected data.
#' spatCorrectedVatorOut <- removeSerieOut(dat = spatCorrectedVator,
#'                                         serieOut = outVator)
#'
#' @family functions for detecting outliers for series of observations
#'
#' @export
removeSerieOut <- function(dat = NULL,
                           fitSpline = NULL,
                           serieOut,
                           trait = attr(x = serieOut,
                                        which = "trait")) {
  ## Check that one of dat and fitSpline are specified.
  if ((is.null(dat) && is.null(fitSpline)) || (
    !is.null(dat) && !is.null(fitSpline))) {
    stop("Specify exactly one of dat and fitSpline as inputs.\n")
  }
  if (!is.null(dat)) {
    if (!inherits(dat, "data.frame")) {
      stop("dat should be a data.frame.\n")
    }
    if (!hasName(dat, "plotId")) {
      stop("dat should at least contain the column plotId.\n")
    }
  }
  if (!is.null(fitSpline)) {
    if (!inherits(fitSpline, "HTPSpline")) {
      stop("dat should be an object of class HTPSpline.\n")
    }
  }
  if (!inherits(serieOut, "data.frame")) {
    stop("serieOut should be a data.frame.\n")
  }
  if (nrow(serieOut) > 0) {
    if (!hasName(serieOut, "plotId")) {
      stop("serieOut should at least contain the column plotId.\n")
    }
    if (!is.null(dat)) {
      ## Remove plots that are in serieOut.
      dat[dat[["plotId"]] %in% serieOut[["plotId"]], trait] <- NA
    } else if (!is.null(fitSpline)) {
      fitSpline$coefDat[fitSpline$coefDat[["plotId"]] %in%
                          serieOut[["plotId"]], "obj.coefficients"] <- NA
      fitSpline$predDat[fitSpline$predDat[["plotId"]] %in%
                          serieOut[["plotId"]], c("pred.value", "deriv")] <- NA
      modDat <- attr(x = fitSpline, which = "modDat")
      modDat[modDat[["plotId"]] %in% serieOut[["plotId"]], trait] <- NA
      attr(x = fitSpline, which = "modDat") <- modDat
    }
  }
  res <- if (!is.null(dat)) dat else fitSpline
  return(res)
}