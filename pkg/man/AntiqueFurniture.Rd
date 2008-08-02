\name{AntiqueFurniture}
\alias{AntiqueFurniture}
\docType{data}
\title{Antique furniture data set}
\description{
  Microscopic measurements on tree species used to make Old-Cape furniture in South Africa between 1652 and 1900.
}
\usage{data(AntiqueFurniture)}
\format{
  A data frame with 37 observations on the following 7 variables.
  \describe{
    \item{\code{Species}}{The tree species. A factor with levels \code{Obul} \code{Oken} \code{Opor}. }
    \item{\code{VesD}}{Tangential vessel diameter in micrometres. A numeric vector.  }
    \item{\code{VesL}}{Vessel element length in micrometres. A numeric vector. }
    \item{\code{FibL}}{Fibre length in micrometres. A numeric vector.}
    \item{\code{RayH}}{Ray height in micrometres. A numeric vector.}
    \item{\code{RayW}}{Ray width in micrometres. A numeric vector.}
    \item{\code{NumVes}}{The number of vessels per millimetre squared. A numeric vector.}
  }
}
\details{
  During the period 1652--1900, wood from both the indiginous \emph{Ocotea bullata} (`Stinkwood') and the imported \emph{Ocotea perosa} (`Imbuia') were used to make Old-Cape furniture in South Africa. The data set contains mean measurements made on such wood, together with a third species, \emph{Ocotea kenyensis}. Twenty samples of \emph{Ocotea bullata} (\code{Obul}), 10 samples of \emph{Ocotea perosa} (\code{Opor}) and 7 samples of \emph{Ocotea kenyensis} (\code{Oken}) were inspected microscopically, with eight variables measured on each sample 50 times over. The data are the mean measurements over the 50 repetitions.
}
\source{
  Swart JPJ (1985), \emph{Lauracea}. Unpublished Masters thesis, Stellenbosch University, Stellenbosch, South Africa. 
}
\references{
  Burden M, Gardner S, Le Roux NJ, Swart JPJ (2001). ``Ou-Kaapse meubels en stinkhoutidentifikasie: Moontlikhede met kanoniese veranderlike-analise en bistippings.'' \emph{South African Journal of Cultural History}, \bold{15}, 50--73.
  
  Le Roux NJ, Gardner S (2005). ``Analysing your multivariate data as a pictorial: A case for applying biplot methodology?'' \emph{International Statistical Review}, \bold{73}(3), 365--387
}
\examples{
data(AntiqueFurniture)
# Paste the following into an interactive R session (without the ##):
## Biplots(Data = AntiqueFurniture[, -1], groups = AntiqueFurniture[, 1])
}
\keyword{datasets}
