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
  \item{excel}{ If the data are to be extracted from an Excel 1997--2003 file, the path to and name of the file. If provided, \code{Data} is ignored. Requires the \pkg{xlsReadWrite} package. Only available under Windows.}
  \item{ExcelGroupsCol}{ If the argument \code{excel} is used and there are multiple groups, the column number of the worksheet column in which the group specifications are given. If sample names are provided in the first column of the worksheet, this column is not counted. If provided, \code{groups} is ignored. Requires the \pkg{xlsReadWrite} package. Only available under Windows.}
}
\details{
  \code{Biplots} is the sole function of the \pkg{BiplotGUI} package. The function initialises the GUI for a particular data set. All further options are available from within the GUI. The GUI features themselves are documented in a separate manual available from the Help menu of the GUI. The manual is also included as a vignette to the package.
  
  The import of data from Excel is rudimentary. It requires the \pkg{xlsReadWrite} package, and is only available under Windows. If provided, \code{Biplots} simply passes the file name and path, \code{excel}, to the \code{read.xls} function of the \pkg{xlsReadWrite} package. Therefore the data should be in the first worksheet of the Excel file. Variable names must be provided in the first row. These are taken to be the axis labels. Sample names may optionally be provided in the first column. These should be textual, not numerical. If provided, the sample names are taken to be the point labels. When read from an Excel file, the sample names must be unique. If sample names are provided, cell A1 should be left blank. Excel 2007 \file{.xlsx} files are not supported.
  
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
