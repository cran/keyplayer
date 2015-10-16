\name{keyplayer-package}
\alias{keyplayer}
\docType{package}
\title{
Locating Key Players in Social Networks
}
\description{
\pkg{keyplayer} is developed in three aspects. First, it introduces four centrality measures built on the recent research: \code{\link{mreach.degree}}, \code{\link{mreach.closeness}} (An and Liu, 2015), \code{\link{fragment}} (Borgatti, 2006), and \code{\link{diffusion}} (Banerjee et.al., 2013). Second, it extends the degree, closeness, betweenness, and eigenvector centralities as well as the above four newly included measures from a node-level to a group-level. \code{\link{kpcent}} encompasses all the eight measures and is able to compute the group centrality when measure and grouping method are specified. Third, it implements a greedy search algorithm to locate key players in a social network. \code{\link{kpset}} reports the players who form the most central set and its centrality score, when measure and target group size are spcified. Please cite the package as follows. \cr

  An, Weihua and Yu-Hsin Liu (2015). ``keyplayer: Locating Key Players in Social Networks." R package version 1.0.1.
}

\details{
\tabular{ll}{
Package: \tab keyplayer\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2015-10-15\cr
Imports: \tab sna, matpow\cr
License: \tab GPL-3\cr
}
}

\author{
Weihua An \email{weihuaan@indiana.edu};
Yu-Hsin Liu \email{yuhsliu@indiana.edu}
}
\references{
An, Weihua. (2015). "Multilevel Meta Network Analysis with Application to Studying Network Dynamics of Network Interventions." \emph{Social Networks} 43: 48-56.\cr

An, Weihua and Yu-Hsin Liu (2015). "keyplayer: An R Package for Locating Key Players in Social Networks."
\emph{Working Paper}, Indiana Univeristy.\cr

Banerjee, A., A. Chandrasekhar, E. Duflo, and M. Jackson (2013): "Diffusion of Microfinance," \emph{Science}, Vol. 341. p.363\cr

Borgatti, Stephen P. (2006). "Identifying Sets of Key Players in a Network." \emph{Computational, Mathematical and Organizational Theory}, 12(1):21-34.\cr

Butts, Carter T. (2014). sna: Tools for Social Network Analysis. R package
  version 2.3-2. \url{http://CRAN.R-project.org/package=sna}\cr

Matloff, Norm and Jack Norman. (2014). matpow: matrix powers. R package
  version 0.1.1. \url{http://CRAN.R-project.org/package=matpow}\cr
}

\keyword{ network key players }
\keyword{ group centrality }


