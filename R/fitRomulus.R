#' Fit the Romulus model.
#'
#' @param cuts1 integer matrix of forward strand DNase I cuts, with a row for each candidate binding site. The columns should correspond to the genomic locations upstream and within the candidate binding site.
#' @param cuts2 integer matrix of reverse strand DNase I cuts, with a row for each candidate binding site. The columns should correspond to the genomic locations within the candidate binding site and downstream.
#' @param anno data frame with annotations of the candidate binding sites. The numeric columns to be used are specified in \code{priors}.
#' @param priors character vector or list of character vectors specifying the names of columns from \code{anno} to be considered in the logistic regression. If a list, each item specifies the column names for each bound state separately, otherwise the same column names will be used for all the bound states.
#' @param bins1 integer vector or matrix specifying how the columns of \code{cuts1} are grouped into bins. If a matrix, each row specifies the grouping for each bound state separately, otherwise the same grouping will be used for all the bound states. The numbers specifying the grouping must form the set of the first N natural numbers (\eqn{1}, \eqn{\dots}, \eqn{N}) for some \eqn{N}.
#' @param bins2 integer vector or matrix specifying how the columns of \code{cuts2} are grouped into bins, as above.
#' @param nbound optional integer, specifying the number of bound states. If not provided, the value will be guessed from the length of \code{priors} or the number of rows in \code{bins1} and \code{bins2}.
#' @param PriorLik optional numeric matrix, with a row for each candidate binding site and a column for each bound state, containing the initial prior likelihoods. If not provided, the default initialization procedure will be used.
#' @param addIntercept logical. Should an additional intercept term be included in the model?
#' @param maxIter integer. Maximal number of Expectation-Maximization iterations to perform.
#' @param maxPostProbDiff numeric. The Expectation-Maximization procedure will be terminated when the absolute differences in posterior probabilities between the iterations will become smaller than this value.
#'
#' @details
#' Fit the Romulus model using an Expectation-Maximization approach.
#'
#' @return A list with the elements
#' \item{nstates}{number of states (including the unbound state).}
#' \item{nbound}{number of bound states.}
#' \item{bins1}{matrix specifying how the columns of \code{cuts1} are grouped into bins, with a row for each state.}
#' \item{bins2}{matrix specifying how the columns of \code{cuts2} are grouped into bins, with a row for each state.}
#' \item{binsizes1}{list specifying, for each state, how many columns of \code{cuts1} fall into each bin.}
#' \item{binsizes2}{list specifying, for each state, how many columns of \code{cuts2} fall into each bin.}
#' \item{Beta}{list of estimated logistic regression coefficients \eqn{\beta_j^{(k)}} for each state.}
#' \item{NegBinom}{matrix of estimated negative binomial parameters for each state, with columns \code{"NegBinomR1"} \code{"NegBinomR2"}, \code{"NegBinomP1"} and \code{"NegBinomP2"}, representing \eqn{r^{+(k)}}, \eqn{r^{-(k)}}, \eqn{p^{+(k)}} and \eqn{p^{-(k)}}, respectively.}
#' \item{Lambda1}{list of estimated multinomial parameters \eqn{\lambda_b^{+(k)}} for each state.}
#' \item{Lambda2}{list of estimated multinomial parameters \eqn{\lambda_b^{-(k)}} for each state.}
#' \item{LambdaReg1}{matrix of multinomial parameter estimates for forward strand DNase I cuts, smoothed by replacing the fixed-value bins by a piecewise linear function, with a row for each state.}
#' \item{LambdaReg2}{matrix of multinomial parameter estimates for reverse strand DNase I cuts, as above.}
#' \item{PriorProb}{matrix of prior probabilities calculated from the logistic regression component, with a row for each candidate binding site and a column for each state.}
#' \item{LogLikelihood}{matrix of log-likelihoods calculated from the negative binomial and multinomial components, with a row for each candidate binding site and a column for each state.}
#' \item{PostProb}{matrix of posterior probabilities, calculated from the complete model, with a row for each candidate binding site and a column for each state.}
#'
#' @references Jankowski, A., Tiuryn, J. and Prabhakar, S. (2016) Romulus: Robust multi-state identification of transcription factor binding sites from DNase-seq data. Bioinformatics. doi: 10.1093/bioinformatics/btw209
#'
#' @examples
#' # Clip the DNase-seq data for NRSF at 99.9% quantile
#' thresh <- max(quantile(NRSF.cuts, 0.999), 1L)
#' cuts <- pmin(NRSF.cuts, thresh)
#' 
#' # Forward strand cuts only upstream and within the candidate binding site
#' cuts1 <- cuts[, 1:(ncol(cuts) / 2 - NRSF.margin)]
#' # Reverse strand cuts only within the candidate binding site and downstream
#' cuts2 <- cuts[, ncol(cuts) / 2 + (NRSF.margin + 1):(ncol(cuts) / 2)]
#' 
#' # Take 20 bp bins outside the candidate binding site and 1 bp bins within it
#' NRSF.width <- (ncol(cuts) - 4 * NRSF.margin) / 2
#' bins1 <- c(rep(1:10, each = 20), 10 + 1:NRSF.width)
#' bins2 <- c(1:NRSF.width, NRSF.width + rep(1:10, each = 20))
#' 
#' # Fit the Romulus model
#' r.fit <- fitRomulus(cuts1, cuts2, NRSF.anno, list(c("score")), bins1, bins2)
#' 
#' 
#' # Benchmarking of Romulus and CENTIPEDE
#' 
#' \dontrun{
#' library(ROCR)
#' library(CENTIPEDE)
#' c.fit <- fitCentipede(Xlist = list(DNase = as.matrix(NRSF.cuts)),
#'   Y = cbind(1, NRSF.anno$score))
#' 
#' r.pred <- prediction(1 - r.fit$PostProb[, r.fit$nstates],
#'   as.integer(NRSF.anno$signalValue > 0))
#' r.perf <- performance(r.pred, measure = "tpr", x.measure = "fpr")
#' r.auc <- performance(r.pred, measure = "auc")
#' 
#' c.pred <- prediction(c.fit$PostPr, as.integer(NRSF.anno$signalValue > 0))
#' c.perf <- performance(c.pred, measure = "tpr", x.measure = "fpr")
#' c.auc <- performance(c.pred, measure = "auc")
#' 
#' plot(r.perf, col = "red",
#'   main = "NRSF binding predictions benchmarked using ChIP-seq data")
#' lines(c.perf@@x.values[[1]], c.perf@@y.values[[1]], col = "blue")
#' legend("bottomright", col = c("red", "blue"), lty = 1, 
#'   legend = c(sprintf("Romulus, AUC = %0.4f", r.auc@@y.values[[1]]),
#'   sprintf("CENTIPEDE, AUC = %0.4f", c.auc@@y.values[[1]])))
#' }

fitRomulus <- function(cuts1, cuts2, anno, priors, bins1, bins2, nbound = NA,
  PriorLik = NULL, addIntercept = T, maxIter = 100, maxPostProbDiff = 0.001)
{
  stopifnot(is.matrix(cuts1))
  stopifnot(is.matrix(cuts2))
  stopifnot(ncol(cuts1) == ncol(bins1))
  stopifnot(ncol(cuts2) == ncol(bins2))
  stopifnot(nrow(cuts1) == nrow(anno))
  stopifnot(nrow(cuts2) == nrow(anno))
  nr <- nrow(anno)

  if (is.na(nbound))
  {
    if (is.list(priors))
      nbound <- length(priors)
    else if (is.matrix(bins1))
      nbound <- nrow(bins1)
    else if (is.matrix(bins2))
      nbound <- nrow(bins2)
    else
      nbound <- 1
  }

  nstates <- nbound + 1L
  cat(paste0("Fitting Romulus model with ", nstates, " states, i.e. ", nbound, " bound state", ifelse(nbound > 1, "s", ""), " and 1 unbound state.\n"))

  if (!is.list(priors))
    priors <- replicate(nbound, priors, simplify = F)
  if (!is.matrix(bins1))
    bins1 <- t(replicate(nbound, bins1))
  if (!is.matrix(bins2))
    bins2 <- t(replicate(nbound, bins2))

  stopifnot(length(priors) == nbound)
  stopifnot(nrow(bins1) == nbound)
  stopifnot(nrow(bins2) == nbound)

  # adding the unbound state with uniform cut distribution
  bins1 <- rbind(bins1, 1L)
  bins2 <- rbind(bins2, 1L)

  # test if bin numbers are from the sequence 1...numbins[k] in each row
  numbins1 <- apply(bins1, 1, function(v) length(unique(v)))
  for (k in 1:nstates)
    if (!all(bins1[k, ] %in% 1:numbins1[k]))
      stop(paste0("bins1[", k, ", ] must form the set of the first N natural numbers (1, ..., N) for some N"))
  numbins2 <- apply(bins2, 1, function(v) length(unique(v)))
  for (k in 1:nstates)
    if (!all(bins2[k, ] %in% 1:numbins2[k]))
      stop(paste0("bins2[", k, ", ] must form the set of the first N natural numbers (1, ..., N) for some N"))

  result <- list()
  result$nstates <- nstates
  result$nbound <- nbound
  result$bins1 <- bins1
  result$bins2 <- bins2

  for (k in 1:nbound)
    if (!all(priors[[k]] %in% colnames(anno)))
      stop(paste0("'priors' refers to column names not present in 'anno'"))
  cols <- lapply(1:nbound, function(k) colnames(anno) %in% priors[[k]])

  cat("Priors and numbers of parameters for forward+reverse strand DNase I footprints:\n")
  for (k in 1:nbound)
    cat(paste0("state ", k, ": priors ", ifelse(addIntercept, '"(Intercept)", ', ''),
      paste0('"', colnames(anno)[cols[[k]]], '"', collapse = ", "), ",\n  ",
      "2+2 negative binomial parameters, ", numbins1[k], "+", numbins2[k], " multinomial parameters\n"))
  cat(paste0("state ", nstates, ": no priors,\n  ",
    "2+2 negative binomial parameters, ", numbins1[nstates], "+", numbins2[nstates], " multinomial parameters\n"))

  lgammaCuts1 <- lgamma(cuts1 + 1)
  lgammaCuts2 <- lgamma(cuts2 + 1)
  sumLgammaCuts1 <- apply(lgammaCuts1, 1, sum)
  sumLgammaCuts2 <- apply(lgammaCuts2, 1, sum)
  sumcuts1 <- apply(cuts1, 1, sum)
  sumcuts2 <- apply(cuts2, 1, sum)
  lgammaSumcuts1 <- lgamma(sumcuts1 + 1)
  lgammaSumcuts2 <- lgamma(sumcuts2 + 1)

  bincuts1 <- lapply(1:nstates, function(k) sapply(1:numbins1[k], function(bin) apply(cuts1[, bins1[k, ] == bin, drop = F], 1, sum)))
  bincuts2 <- lapply(1:nstates, function(k) sapply(1:numbins2[k], function(bin) apply(cuts2[, bins2[k, ] == bin, drop = F], 1, sum)))
  binsizes1 <- lapply(1:nstates, function(k) sapply(1:numbins1[k], function(bin) sum(bins1[k, ] == bin)))
  binsizes2 <- lapply(1:nstates, function(k) sapply(1:numbins2[k], function(bin) sum(bins2[k, ] == bin)))
  result$binsizes1 <- binsizes1
  result$binsizes2 <- binsizes2

  if (is.null(PriorLik))
  {
    PriorLik <- matrix(0.01, nr, nbound)
    # initial estimate of bound states
    scsel <- sumcuts1 + sumcuts2 > quantile(sumcuts1 + sumcuts2, probs = 0.9)

    for (uc in unique(cols))
    {
      st <- which(sapply(cols, identical, uc))
      if (length(st) < 2)
        # if only one state has this set of priors
        PriorLik[scsel, st] <- 100 # apply the initial estimate
      else
      {
        # if priors are the same for multiple states, differentiate the states by kmeans clustering
        cat(paste0("Applying k-means clustering to differentiate tkmeahe initial prior probabilities of states: ",
          paste0(st, collapse = ", "), "\n"))
        scselcuts <- cbind(bincuts1[[st[1]]][scsel, ], bincuts2[[st[1]]][scsel, ])
        scselclust <- kmeans(scselcuts, length(st), nstart = 100)$cluster
        for (cl in seq_along(st))
          PriorLik[scsel, st[cl]][scselclust == cl] <- 100 # apply the initial estimate
        rm(scselcuts, scselclust)
      }
    }
    rm(scsel)
  }
  stopifnot(nrow(PriorLik) == nrow(anno))
  stopifnot(ncol(PriorLik) == nbound)
  PriorLogLik <- log(PriorLik)

  Beta <- list()
  for (k in 1:nbound)
  {
    if (addIntercept)
      thislm <- lm(PriorLogLik[, k] ~ ., data = as.data.frame(anno[, cols[[k]], drop = F]))
    else
      thislm <- lm(PriorLogLik[, k] ~ . - 1, data = as.data.frame(anno[, cols[[k]], drop = F]))

    Beta[[k]] <- thislm$coefficients
    # if some of the coefficients could not be fitted, then they are NAs; replace them with 0s
    Beta[[k]][is.na(Beta[[k]])] <- 0

  #  PriorLogLik[, k] <- extanno[, extcols[[k]], drop = F] %*% Beta[[k]] # equivalent to predict(thislm)
  }
  cat("\nInitial logistic regression parameters:\n")
  print(Beta)

  colsel <- do.call(pmax, cols) == 1
  if (addIntercept)
  {
    extanno <- cbind(1, as.matrix(anno[, colsel]))
    extcols <- lapply(cols, function(v) c(T, v[colsel]))
  }
  else
  {
    extanno <- as.matrix(anno[, colsel])
    extcols <- lapply(cols, function(v) v[colsel])
  }

  PriorLogLik <- sapply(1:nbound, function(k) extanno[, extcols[[k]], drop = F] %*% Beta[[k]])
  PriorLik <- exp(PriorLogLik)
  PriorProbUnbound <- 1 / (1 + apply(PriorLik, 1, sum))
  PriorProb <- cbind(PriorLik * PriorProbUnbound, PriorProbUnbound, deparse.level = 0)

  # only for the first iteration
  PostProb <- PriorProb
  OldPostProb <- NA


  for (iter in 1:maxIter)
  {
    cat(paste("\nIteration: ", iter, "\n", sep = ""))

    # now the maximization step

    # LogLC1: part of complete likelihood with Beta

    #print("LogLC1")

    fnLogLC1 <- function(BetaUnlisted)
    {
      thisBeta <- relist(BetaUnlisted, Beta)
      thisPriorLogLik <- sapply(1:nbound, function(k) extanno[, extcols[[k]], drop = F] %*% thisBeta[[k]])
      thisPriorLik <- exp(thisPriorLogLik)

      return(sum(PostProb[, -nstates] * thisPriorLogLik) - nbound * sum(log(1 + apply(thisPriorLik, 1, sum))))
    }

    grLogLC1 <- function(BetaUnlisted)
    {
      thisBeta <- relist(BetaUnlisted, Beta)
      thisPriorLogLik <- sapply(1:nbound, function(k) extanno[, extcols[[k]], drop = F] %*% thisBeta[[k]])
      thisPriorLik <- exp(thisPriorLogLik)

      grLogLC1k <- function(k)
      {
        s1 <- crossprod(PostProb[, k], extanno[, extcols[[k]], drop = F])
        s2 <- - crossprod(thisPriorLik[, k] / (1 + apply(thisPriorLik, 1, sum)), extanno[, extcols[[k]], drop = F])
        return(s1 + s2)
      }

      return(do.call(c, lapply(1:nbound, grLogLC1k)))
    }

    opt <- optim(unlist(Beta), fnLogLC1, grLogLC1, method = "BFGS", control = list(fnscale = -1))
    Beta <- relist(opt$par, Beta)
    LogLC1 <- opt$value
    cat("\nUpdated logistic regression parameters:\n")
    print(Beta)
    #print(LogLC1)

    # LogLC2: part of complete likelihood with NegBinomP and NegBinomR

    #print("LogLC2")

    sumPostProb <- apply(PostProb, 2, sum)
    sumPostProbSumcuts1 <- apply(PostProb * sumcuts1, 2, sum)
    sumPostProbSumcuts2 <- apply(PostProb * sumcuts2, 2, sum)
    sumPostProbLgammaSumcuts1 <- apply(PostProb * lgammaSumcuts1, 2, sum)
    sumPostProbLgammaSumcuts2 <- apply(PostProb * lgammaSumcuts2, 2, sum)

    fnLogLC21 <- function(thisNegBinomR, k)
    {
      if (!(thisNegBinomR > 0)) return(NA)
      thisNegBinomP <- thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts1[k])

      s1 <- sum(PostProb[, k] * lgamma(thisNegBinomR + sumcuts1))
      s2 <- - sumPostProb[k] * lgamma(thisNegBinomR)
      s3 <- - sumPostProbLgammaSumcuts1[k]
      s4 <- sumPostProb[k] * thisNegBinomR * log(thisNegBinomP)
      s5 <- sumPostProbSumcuts1[k] * log(1 - thisNegBinomP)

      return(s1 + s2 + s3 + s4 + s5)
    }

    fnLogLC22 <- function(thisNegBinomR, k)
    {
      if (!(thisNegBinomR > 0)) return(NA)
      thisNegBinomP <- thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts2[k])

      s1 <- sum(PostProb[, k] * lgamma(thisNegBinomR + sumcuts2))
      s2 <- - sumPostProb[k] * lgamma(thisNegBinomR)
      s3 <- - sumPostProbLgammaSumcuts2[k]
      s4 <- sumPostProb[k] * thisNegBinomR * log(thisNegBinomP)
      s5 <- sumPostProbSumcuts2[k] * log(1 - thisNegBinomP)

      return(s1 + s2 + s3 + s4 + s5)
    }

    grLogLC21 <- function(thisNegBinomR, k)
    {
      if (!(thisNegBinomR > 0)) return(NA)
      thisNegBinomP <- thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts1[k])

      s1 <- sum(PostProb[, k] * digamma(thisNegBinomR + sumcuts1))
      s2 <- - sumPostProb[k] * digamma(thisNegBinomR)
      s4 <- sumPostProb[k] * (1 - thisNegBinomP + log(thisNegBinomP))
      s5 <- - sumPostProbSumcuts1[k] * thisNegBinomP / thisNegBinomR

      return(s1 + s2 + s4 + s5)
    }

    grLogLC22 <- function(thisNegBinomR, k)
    {
      if (!(thisNegBinomR > 0)) return(NA)
      thisNegBinomP <- thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts2[k])

      s1 <- sum(PostProb[, k] * digamma(thisNegBinomR + sumcuts2))
      s2 <- - sumPostProb[k] * digamma(thisNegBinomR)
      s4 <- sumPostProb[k] * (1 - thisNegBinomP + log(thisNegBinomP))
      s5 <- - sumPostProbSumcuts2[k] * thisNegBinomP / thisNegBinomR

      return(s1 + s2 + s4 + s5)
    }

    NegBinomR1 <- rep(1, nstates)
    NegBinomR2 <- rep(1, nstates)
    LogLC2 <- 0
    for (k in 1:nstates)
    {
      opt1 <- optim(NegBinomR1[k], fnLogLC21, grLogLC21, k = k, method = "BFGS", control = list(fnscale = -1))
      opt2 <- optim(NegBinomR2[k], fnLogLC22, grLogLC22, k = k, method = "BFGS", control = list(fnscale = -1))
      NegBinomR1[k] <- opt1$par
      NegBinomR2[k] <- opt2$par
      LogLC2 <- LogLC2 + opt1$value + opt2$value
    }
    NegBinomP1 <- NegBinomR1 * sumPostProb / (NegBinomR1 * sumPostProb + sumPostProbSumcuts1)
    NegBinomP2 <- NegBinomR2 * sumPostProb / (NegBinomR2 * sumPostProb + sumPostProbSumcuts2)

    cat("Updated negative binomial parameters for DNase I footprint:\n")
    print(cbind(NegBinomR1, NegBinomR2, NegBinomP1, NegBinomP2))
    #print(LogLC2)
    #LogLC2 <- sum(sapply(1:nstates, function(k) sum(PostProb[, k] * dnbinom(sumcuts, NegBinomR[k], NegBinomP[k], log = T))))
    #print(LogLC2)

    # LogLC3: part of complete likelihood with Lambda[[i]]

    regularize <- function(values, bins)
    {
      rl <- rle(bins)
      i <- 0

      for (j in seq_len(length(rl$values) - 1))
      {
        i <- i + rl$lengths[j]
        il <- i - (rl$lengths[j] - 1) %/% 2
        ir <- i + (rl$lengths[j + 1] + 1) %/% 2
        values[il:ir] <- seq(from = values[il], to = values[ir], length.out = ir - il + 1)
      }

      return(values)
    }

    #print("LogLC3")
    Lambda1 <- lapply(1:nstates, function(k)
    {
      lambda <- as.vector(crossprod(PostProb[, k], bincuts1[[k]])) / binsizes1[[k]]
      lambda <- lambda + sum(lambda * binsizes1[[k]]) / ncol(bins1) # estimator shrinkage
      lambda / sum(lambda * binsizes1[[k]])
    })
    #print(Lambda1)

    LambdaReg1 <- t(sapply(1:nstates, function(k) Lambda1[[k]][bins1[k, ]]))
    for (k in 1:nstates)
      LambdaReg1[k, ] <- regularize(LambdaReg1[k, ], bins1[k, ])

    Lambda2 <- lapply(1:nstates, function(k)
    {
      lambda <- as.vector(crossprod(PostProb[, k], bincuts2[[k]])) / binsizes2[[k]]
      lambda <- lambda + sum(lambda * binsizes2[[k]]) / ncol(bins2) # estimator shrinkage
      lambda / sum(lambda * binsizes2[[k]])
    })
    #print(Lambda2)
    cat("\nUpdated multinomial parameters for DNase I footprint (not shown).\n")

    LambdaReg2 <- t(sapply(1:nstates, function(k) Lambda2[[k]][bins2[k, ]]))
    for (k in 1:nstates)
      LambdaReg2[k, ] <- regularize(LambdaReg2[k, ], bins2[k, ])

    #sapply(1:nstates, function(k) print(sum(Lambda[[k]][bins[k, ]]) - 1))

    LogLC3 <- sum(sapply(1:nstates, function(k)
    {
      s1 <- sum(PostProb[, k] * lgammaSumcuts1) + sum(PostProb[, k] %*% cuts1 * log(LambdaReg1[k, ])) - sum(PostProb[, k] * sumLgammaCuts1)
      s2 <- sum(PostProb[, k] * lgammaSumcuts2) + sum(PostProb[, k] %*% cuts2 * log(LambdaReg2[k, ])) - sum(PostProb[, k] * sumLgammaCuts2)
      s1 + s2
    }))
    #print(LogLC3)

    #LogLC3alt <- 0
    #for (i in seq_len(nr))
    #{
    #  LogLC3alt <- LogLC3alt + PostProb[, nstates][i] * dmultinom(cuts[i, ], prob = Lambda[[nstates]][bins[1, ]], log = T)
    #  LogLC3alt <- LogLC3alt + PostProb[i] * dmultinom(cuts[i, ], prob = Lambda[[1]][bins[1, ]], log = T)
    #}
    #print(LogLC3alt)

    LogLC <- LogLC1 + LogLC2 + LogLC3
    #print("LogLC total")
    #print(LogLC)


    # now the maximization step

    PriorLogLik <- sapply(1:nbound, function(k) extanno[, extcols[[k]], drop = F] %*% Beta[[k]])
    PriorLik <- exp(PriorLogLik)
    PriorProbUnbound <- 1 / (1 + apply(PriorLik, 1, sum))
    PriorProb <- cbind(PriorLik * PriorProbUnbound, PriorProbUnbound, deparse.level = 0)

    cat("\nAverage prior probabilites:\n")
    m <- apply(PriorProb, 2, mean)
    names(m) <- paste0("V", seq_along(m))
    print(m)

    LogCondProb11 <- sapply(1:nstates, function(k) dnbinom(sumcuts1, NegBinomR1[k], NegBinomP1[k], log = T))
    LogCondProb12 <- sapply(1:nstates, function(k) dnbinom(sumcuts2, NegBinomR2[k], NegBinomP2[k], log = T))
    LogCondProb21 <- sapply(1:nstates, function(k) lgammaSumcuts1 + cuts1 %*% log(LambdaReg1[k, ]) - sumLgammaCuts1)
    LogCondProb22 <- sapply(1:nstates, function(k) lgammaSumcuts2 + cuts2 %*% log(LambdaReg2[k, ]) - sumLgammaCuts2)

    # here we use the Bayes theorem
    LogLikelihood <- LogCondProb11 + LogCondProb12 + LogCondProb21 + LogCondProb22
    # avoiding numerical overflow
    JointProb <- PriorProb * exp(LogLikelihood - apply(LogLikelihood, 1, max))
    PostProb <- JointProb / apply(JointProb, 1, sum)

    cat("\nAverage posterior probabilites:\n")
    m <- apply(PostProb, 2, mean)
    names(m) <- paste0("V", seq_along(m))
    print(m)

    if (iter > 1)
    {
      maxdiff <- max(abs(OldPostProb - PostProb))
      cat(paste0("\nMaximum difference in posterior probability: ", format(maxdiff), "\n"))
      if (maxdiff < maxPostProbDiff) break
    }
    OldPostProb <- PostProb

    result$Beta <- Beta
    result$NegBinom <- cbind(NegBinomR1, NegBinomR2, NegBinomP1, NegBinomP2)
    result$Lambda1 <- Lambda1
    result$Lambda2 <- Lambda2
    result$LambdaReg1 <- LambdaReg1
    result$LambdaReg2 <- LambdaReg2
    cat("\n")
  }

  result$PriorProb <- PriorProb
  result$LogLikelihood <- LogLikelihood
  result$PostProb <- PostProb
  return(result)
}
