\name{Countries}
\alias{Countries}
\docType{data}
\title{Countries data set}
\description{
  Eight variables measured on the countries with the 15 largest economies by purchasing price parity (PPP) gross domesitic product (GDP) in 2007.
}
\usage{data(Countries)}
\format{
A data frame with 15 observations on the following 8 variables.
  \describe{
    \item{\code{GDP}}{Purchasing price parity (PPP) gross domesitc product (GDP) per capita in 2007 US dollars. A numeric vector.}
    \item{\code{HIV.Aids}}{HIV/Aids prevalence as a percentage of the population. A numeric vector.}
    \item{\code{Life.exp.}}{Life expectancy in years. A numeric vector.}
    \item{\code{Mil.}}{Military spending as a percentage of GDP. A numeric vector.}
    \item{\code{Oil.cons.}}{Oil consumption in barrels per annum per capita. A numeric vector.}
    \item{\code{Pop.}}{Population in millions. A numeric vector.}
    \item{\code{Tel.}}{Number of fixed line telephones per 1000 people. A numeric vector.}
    \item{\code{Unempl.}}{Percentage unemployed. A numeric vector.}
  }}
\details{
  The data have been derived largely from the 2007 CIA World Factbook, and are intended for illustrative purposes only.
}
\source{
  Agency C (2007). \emph{The World Factbook: 2007, CIA's 2006}, Potomac Books, Washington, DC, USA.
}
\examples{
data(Countries)
\dontrun{Biplots(Data = Countries)}
}
\keyword{datasets}
