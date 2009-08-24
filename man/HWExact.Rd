\name{HWExact}
\Rdversion{1.1}
\alias{HWExact}
\title{Function to calculate Hardy-Weinburg exact p-values
}
\description{
 This function calculates Hardy-Weinburg (Fisher's) exact p-values for GWAS SNP data.
}
\usage{
HWExact(GenotypeCounts)
}

\arguments{
  \item{GenotypeCounts}{A dataframe of genotype counts, with columns called \code{nAA}, \code{nAa} and \code{naa}, one row for each SNP.}
}

\details{
}

\value{A vector of exact p-values.}

\references{ Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76} 

\author{Ian Painter}

\note{ This function uses a C function \code{SNPHWE.c} written by Jan Wigginton as described in the above reference.}



\seealso{
}
\examples{

pA<- runif(1)
pAA<- pA^2
pAa<- 2*pA*(1-pA)
paa<- (1-pA)^2

counts<- rmultinom(1000, 3000, c(pAA, pAa, paa) )
genotypes<- data.frame(nAA = counts[1,], nAa = counts[2,], naa = counts[3,])
hwPvalues<- HWExact(genotypes)

}

\keyword{ htest }