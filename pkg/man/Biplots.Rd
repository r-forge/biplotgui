\name{Biplots}
\alias{Biplots}
\title{Interactive Biplots in R}
\description{
  Initialises a GUI with which to construct and interact with biplots.
}
\usage{%
  Biplots(Data, groups = rep(1, nrow(Data)),  
          PointLabels = rownames(Data),  
          AxisLabels = colnames(Data), excel = NULL,  
          ExcelGroupsCol = 0)
}
\arguments{
  \item{Data}{ A matrix or data frame of numerical data. Its \eqn{n}{n} samples (observations) will be represented as points in the biplots; its \eqn{p}{p} variables will be represented as calibrated biplot axes. }
  \item{groups}{ A vector or factor of length \eqn{n}{n} specifying the group membership of the samples. By default, all samples are taken to be from a single group. The group labels are taken from this argument, deprecated to 14 characters each.}
  \item{PointLabels}{ A vector of length \eqn{n}{n} specifying the point labels. By default, the point labels are taken to be the row names of \code{Data}. }
  \item{AxisLabels}{ A vector of length \eqn{p}{p} specifying the axis labels. By default, the axis labels are taken to be the column names of \code{Data}, deprecated to 14 characters each.}
  \item{excel}{ Deprecated as from version 0.0-4.1.}
  \item{ExcelGroupsCol}{ Deprecated as from version 0.0-4.1.}
}
\details{
  \code{Biplots} is the sole function of the \pkg{BiplotGUI} package. The function initialises the GUI for a particular data set. All further options are available from within the GUI. The GUI features themselves are documented in a separate manual available from the Help menu of the GUI. The manual is also included as a vignette to the package.
  
  Due to the removal of the \pkg{xlsReadWrite} package from CRAN, the direct import of data from Excel 1997-2003 files has been deprecated as from \pkg{BiplotGUI} 0.0-4.1. As an alternative mechanism, consider the \pkg{RODBC} package.
  
At present, \pkg{BiplotGUI} is intended to be run under Windows. In Windows, it runs marginally better in SDI mode, rather than MDI mode. 
}
\references{ 
  Gower JC, Hand DJ (1996). \emph{Biplots}. Monographs on Statistics and Applied Probability. Chapman \& Hall, London, UK.
}
\author{
Author and maintainer: Anthony la Grange <amlg at sun.ac.za>, <\url{http://biplotgui.r-forge.r-project.org/}>
}
\examples{
data(Countries)
\dontrun{Biplots(Data = Countries)}

data(AntiqueFurniture)
\dontrun{%
Biplots(Data = AntiqueFurniture[, -1],  
        groups = AntiqueFurniture[, 1])}

data(FighterAircraft)
\dontrun{Biplots(Data = FighterAircraft)}
}
\keyword{dplot}
\keyword{hplot}
\keyword{iplot}
\keyword{dynamic}
\keyword{multivariate}
