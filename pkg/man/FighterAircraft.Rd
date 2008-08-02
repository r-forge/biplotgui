\name{FighterAircraft}
\alias{FighterAircraft}
\docType{data}
\title{ Fighter aircraft data set }
\description{
  Eight variables measured on 21 types of US fighter aircraft.
}
\usage{data(FighterAircraft)}
\format{
  A data frame with 21 observations on the following 4 variables.
  \describe{
    \item{\code{SPR}}{Specific power, proportional to power per unit weight. A numeric vector.}
    \item{\code{RGF}}{Flight range factor. A numeric vector.}
    \item{\code{PLF}}{Payload as a fraction of gross weight. A numeric vector.}
    \item{\code{SLF}}{Sustained load factor. A numeric vector.}
  }
}
\details{
  Measurements of four variables on 21 of 22 types of US fighter aircraft extracted by Cook and Weisberg (1982) from a report by Stanley and Miller (1979).
}
\source{
  Stanley W, Miller M (1979). ``Measuring technological change in jet fighter aircraft.'' \emph{Technical Report R-2249-AF}, RAND Corporation, Santa Monica, CA, USA.
}
\references{
  Cook RD, Weisberg S (1982). \emph{Residuals and influence in regression}. Monographs on Statistics and Applied Probability. Chapman \& Hall, London, UK.
  
  Gower JC, Hand DJ (1996). \emph{Biplots}. Monographs on Statistics and Applied Probability. Chapman \& Hall, London, UK.
}
\examples{
data(FighterAircraft)
# Paste the following into an interactive R session (without the ##):
## Biplots(Data = FighterAircraft)
}
\keyword{datasets}
