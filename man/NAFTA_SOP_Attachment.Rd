\name{NAFTA_SOP_Attachment}
\alias{NAFTA_SOP_Attachment}
\docType{data}
\title{
Example datasets from Attachment 1 to the NAFTA SOP published 2015
}
\description{
Data taken from from Attachment 1 of the SOP.
}
\usage{NAFTA_SOP_Attachment}
\format{
  A list (NAFTA_SOP_Attachment) containing 16 datasets suitable
  for the evaluation with \code{\link{nafta}}
}
\source{
  NAFTA (2011) Guidance for evaluating and calculating degradation kinetics
  in environmental media. NAFTA Technical Working Group on Pesticides
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/guidance-evaluating-and-calculating-degradation}
  accessed 2019-02-22

  US EPA (2015) Standard Operating Procedure for Using the NAFTA Guidance to
  Calculate Representative Half-life Values and Characterizing Pesticide
  Degradation
  \url{https://www.epa.gov/pesticide-science-and-assessing-pesticide-risks/standard-operating-procedure-using-nafta-guidance}
}
\examples{
  nafta_att_p5a <- nafta(NAFTA_SOP_Attachment[["p5a"]], cores = 1)
  print(nafta_att_p5a)
  plot(nafta_att_p5a)
}
\keyword{datasets}
