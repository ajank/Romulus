% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitRomulus.R
\name{fitRomulus}
\alias{fitRomulus}
\title{Fit the Romulus model.}
\usage{
fitRomulus(cuts1, cuts2, anno, priors, bins1, bins2, nbound = NA,
  PriorLik = NULL, addIntercept = T, mixingDelta = 0.5, maxIter = 100,
  maxPostProbDiff = 0.001)
}
\arguments{
\item{cuts1}{integer matrix of forward strand DNase I cuts, with a row for each candidate binding site. The columns should correspond to the genomic locations upstream and within the candidate binding site.}

\item{cuts2}{integer matrix of reverse strand DNase I cuts, with a row for each candidate binding site. The columns should correspond to the genomic locations within the candidate binding site and downstream.}

\item{anno}{data frame with annotations of the candidate binding sites. The numeric columns to be used are specified in \code{priors}.}

\item{priors}{character vector or list of character vectors specifying the names of columns from \code{anno} to be considered in the logistic regression. If a list, each item specifies the column names for each bound state separately, otherwise the same column names will be used for all the bound states.}

\item{bins1}{integer vector or matrix specifying how the columns of \code{cuts1} are grouped into bins. If a matrix, each row specifies the grouping for each bound state separately, otherwise the same grouping will be used for all the bound states. The numbers specifying the grouping must form the set of the first N natural numbers (\eqn{1}, \eqn{\dots}, \eqn{N}) for some \eqn{N}.}

\item{bins2}{integer vector or matrix specifying how the columns of \code{cuts2} are grouped into bins, as above.}

\item{nbound}{optional integer, specifying the number of bound states. If not provided, the value will be guessed from the length of \code{priors} or the number of rows in \code{bins1} and \code{bins2}.}

\item{PriorLik}{optional numeric matrix, with a row for each candidate binding site and a column for each bound state, containing the initial prior likelihoods. If not provided, the default initialization procedure will be used.}

\item{addIntercept}{logical. Should an additional intercept term be included in the model?}

\item{maxIter}{integer. Maximal number of Expectation-Maximization iterations to perform.}

\item{maxPostProbDiff}{numeric. The Expectation-Maximization procedure will be terminated when the absolute differences in posterior probabilities between the iterations will become smaller than this value.}
}
\value{
A list with the elements
\item{nstates}{number of states (including the unbound state).}
\item{nbound}{number of bound states.}
\item{bins1}{matrix specifying how the columns of \code{cuts1} are grouped into bins, with a row for each state.}
\item{bins2}{matrix specifying how the columns of \code{cuts2} are grouped into bins, with a row for each state.}
\item{binsizes1}{list specifying, for each state, how many columns of \code{cuts1} fall into each bin.}
\item{binsizes2}{list specifying, for each state, how many columns of \code{cuts2} fall into each bin.}
\item{Beta}{list of estimated logistic regression coefficients \eqn{\beta_j^{(k)}} for each state.}
\item{NegBinom}{matrix of estimated negative binomial parameters for each state, with columns \code{"NegBinomR1"} \code{"NegBinomR2"}, \code{"NegBinomP1"} and \code{"NegBinomP2"}, representing \eqn{r^{+(k)}}, \eqn{r^{-(k)}}, \eqn{p^{+(k)}} and \eqn{p^{-(k)}}, respectively.}
\item{Lambda1}{list of estimated multinomial parameters \eqn{\lambda_b^{+(k)}} for each state.}
\item{Lambda2}{list of estimated multinomial parameters \eqn{\lambda_b^{-(k)}} for each state.}
\item{LambdaReg1}{matrix of multinomial parameter estimates for forward strand DNase I cuts, smoothed by replacing the fixed-value bins by a piecewise linear function, with a row for each state.}
\item{LambdaReg2}{matrix of multinomial parameter estimates for reverse strand DNase I cuts, as above.}
\item{PriorProb}{matrix of prior probabilities calculated from the logistic regression component, with a row for each candidate binding site and a column for each state.}
\item{LogLikelihood}{matrix of log-likelihoods calculated from the negative binomial and multinomial components, with a row for each candidate binding site and a column for each state.}
\item{PostProb}{matrix of posterior probabilities, calculated from the complete model, with a row for each candidate binding site and a column for each state.}
}
\description{
Fit the Romulus model.
}
\details{
Fit the Romulus model using an Expectation-Maximization approach.
}
\examples{
# Clip the DNase-seq data for NRSF at 99.9\% quantile
thresh <- max(quantile(NRSF.cuts, 0.999), 1L)
cuts <- pmin(NRSF.cuts, thresh)

# Forward strand cuts only upstream and within the candidate binding site
cuts1 <- cuts[, 1:(ncol(cuts) / 2 - NRSF.margin)]
# Reverse strand cuts only within the candidate binding site and downstream
cuts2 <- cuts[, ncol(cuts) / 2 + (NRSF.margin + 1):(ncol(cuts) / 2)]

# Take 20 bp bins outside the candidate binding site and 1 bp bins within it
NRSF.width <- (ncol(cuts) - 4 * NRSF.margin) / 2
bins1 <- c(rep(1:10, each = 20), 10 + 1:NRSF.width)
bins2 <- c(1:NRSF.width, NRSF.width + rep(1:10, each = 20))

# Fit the Romulus model
r.fit <- fitRomulus(cuts1, cuts2, NRSF.anno, list(c("score")), bins1, bins2)


# Benchmarking of Romulus and CENTIPEDE

\dontrun{
library(ROCR)
library(CENTIPEDE)
c.fit <- fitCentipede(Xlist = list(DNase = as.matrix(NRSF.cuts)),
  Y = cbind(1, NRSF.anno$score))

r.pred <- prediction(1 - r.fit$PostProb[, r.fit$nstates],
  as.integer(NRSF.anno$signalValue > 0))
r.perf <- performance(r.pred, measure = "tpr", x.measure = "fpr")
r.auc <- performance(r.pred, measure = "auc")

c.pred <- prediction(c.fit$PostPr, as.integer(NRSF.anno$signalValue > 0))
c.perf <- performance(c.pred, measure = "tpr", x.measure = "fpr")
c.auc <- performance(c.pred, measure = "auc")

plot(r.perf, col = "red",
  main = "NRSF binding predictions benchmarked using ChIP-seq data")
lines(c.perf@x.values[[1]], c.perf@y.values[[1]], col = "blue")
legend("bottomright", col = c("red", "blue"), lty = 1,
  legend = c(sprintf("Romulus, AUC = \%0.4f", r.auc@y.values[[1]]),
  sprintf("CENTIPEDE, AUC = \%0.4f", c.auc@y.values[[1]])))
}

}
\references{
Jankowski, A., Tiuryn, J. and Prabhakar, S. (2016).  Romulus: robust multi-state identification of transcription factor binding sites from DNase-seq data.  \emph{Bioinformatics} \bold{32}, 2419--2426. doi: 10.1093/bioinformatics/btw209.
}
