% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idid.R
\name{idid_wald2}
\alias{idid_wald2}
\title{Two-sample Summary-Data Instrumented DID Wald estimator}
\usage{
idid_wald2(mu.y, mu.d, se.y, se.d, tv, zv)
}
\arguments{
\item{mu.y}{A 4-dimensional vector of outcome means, for every (t,z) combination}

\item{mu.d}{A 4-dimensional vector of exposure rates, for every (t,z) combination}

\item{se.y}{A 4-dimensional vector of the standard error of \code{mu.y}, for every (t,z) combination}

\item{se.d}{A 4-dimensional vector of the standard error of \code{mu.d}, for every (t,z) combination}

\item{tv}{A 4-dimensional vector for time}

\item{zv}{A 4-dimensional vector for the level of IV for DID}
}
\value{
\describe{
\item{beta.est}{Two-sample summary-data instrumented DID Wald estimate.}
\item{beta.se}{Standard error of \code{beta.est}.}
\item{f.statistic}{F statistic from the first-stage regression. It is preferably larger than 10.}
}
}
\description{
Two-sample Summary-Data Instrumented DID Wald estimator
}
\examples{
dfa<-data_gen(1e5) # outcome dataset
dfb<-data_gen(1e5) # exposure dataset
dfa<-dfa[dfa$x1<=0&dfa$x2<=0,]
dfb<-dfb[dfb$x1<=0&dfb$x2<=0,]
res<-data.frame(mu.y=0,mu.d=0, se.y=0, se.d=0, t=c(0,1,0,1),z=c(0,0,1,1))
for(i in 1:4){
t<-res$t[i]
z<-res$z[i]
res$mu.y[i]<-mean(dfa$y[dfa$t==t & dfa$z==z])
res$mu.d[i]<-mean(dfb$d[dfb$t==t & dfb$z==z])
res$se.y[i]<-sqrt(var(dfa$y[dfa$t==t & dfa$z==z])/length(dfa$y[dfa$t==t & dfa$z==z]))
res$se.d[i]<-sqrt(var(dfb$d[dfb$t==t & dfb$z==z])/length(dfb$d[dfb$t==t & dfb$z==z]))
}
idid_wald2(mu.y=res$mu.y,mu.d=res$mu.d,se.y=res$se.y,se.d=res$se.d,tv=res$t,zv=res$z)

}
\references{
Ting Ye, Ashkan Ertefaie, James Flory, Sean Hennessy, and Dylan S. Small (2021). Instrumented Difference-in-Differences.
}
