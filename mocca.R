options(echo = T)
args = commandArgs(trailingOnly = T)
prefix = args[1]
print(prefix)
file.prefix = args[2]
print(file.prefix)
dnase.prefix = args[3]
print(dnase.prefix)

library(ROCR)

margin = 200

pred = list()
predr = list()
perf = list()
perfr = list()
auc = list()
prior = list()
loglik = list()
prob = list()
prof = list()
params = list()

{

pdf(paste(file.prefix, ".pdf", sep = ""))
load(paste(dnase.prefix, "/", prefix, ".Rdata", sep = ""))
#cuts = cuts[1:10000, ]
#anno = anno[1:10000, ]

cuts = cuts[, !is.na(cuts[1, ])]
stopifnot(ncol(cuts) %% 2 == 0)
len = ncol(cuts) / 2 - 2 * margin
print(dim(cuts))
print(summary(anno))

thresh = max(quantile(cuts, 0.999), 1)
print(thresh)
cuts[cuts > thresh] = thresh

cuts1 = cuts[, 1:(ncol(cuts) / 2 - margin)]
cuts2 = cuts[, ncol(cuts) / 2 + (margin + 1):(ncol(cuts) / 2)]
rm(cuts)

lgammaCuts1 = lgamma(cuts1 + 1)
lgammaCuts2 = lgamma(cuts2 + 1)
sumLgammaCuts1 = apply(lgammaCuts1, 1, sum)
sumLgammaCuts2 = apply(lgammaCuts2, 1, sum)
sumcuts1 = apply(cuts1, 1, sum)
sumcuts2 = apply(cuts2, 1, sum)
lgammaSumcuts1 = lgamma(sumcuts1 + 1)
lgammaSumcuts2 = lgamma(sumcuts2 + 1)

signalcols = names(anno)[grep("signalValue", names(anno))]
anno$signalValue = 1

compare.proper = function(this_pred, this_anno, col)
{
#  return(cor(pred, anno, method = "spearman"))

  id = paste(prefix, col, K)
  pred[[id]] <<- prediction(this_pred, ifelse(this_anno > 0, 1, 0))
  predr[[id]] <<- prediction(rank(this_pred, ties.method = "random"), ifelse(this_anno > 0, 1, 0))
  perf[[id]] <<- performance(pred[[id]], measure = "tpr", x.measure = "fpr")
  perfr[[id]] <<- performance(predr[[id]], measure = "rec")
  auc[[id]] <<- performance(pred[[id]], measure = "auc")

  plot(perf[[id]], col = rainbow(10), main = id, sub = paste(length(this_anno), "sites, of which", sum(this_anno), "positives"))
  return(auc[[id]]@y.values[[1]])
}

compare = function(pred, col)
{
  print(c(compare.proper(pred, anno[, col], col), 42))
}


#for (K in 0:3)
K = length(grep("score", names(anno))) - 1
{
  # FIXME should accept list of colnames, list of boolean values, list of formulas?
  Priors = list()
  for (k in 1:(K + 1))
    Priors[[k]] = unique(c("score", names(anno)[grep("score", names(anno))][k]))

  bins1 = rbind(c(rep(1:10, each = 20), 10 + 1:len))
  for (k in seq_len(K))
    bins1 = rbind(bins1, c(rep(1:10, each = 20), 10 + 1:len))
  bins1 = rbind(bins1, rep(1, 200 + len))

  bins2 = rbind(c(1:len, len + rep(1:10, each = 20)))
  for (k in seq_len(K))
    bins2 = rbind(bins2, c(1:len, len + rep(1:10, each = 20)))
  bins2 = rbind(bins2, rep(1, 200 + len))


### here main MOCCA code

# assume we have K
zero = K + 2
nr = nrow(anno)
stopifnot(nrow(cuts1) == nr)
stopifnot(nrow(cuts2) == nr)

id = paste(prefix, K)
params[[id]] = list()
params[[id]]$K = K
params[[id]]$zero = zero
params[[id]]$bins1 = bins1
params[[id]]$bins2 = bins2

# genomic sites go in rows
# subsequent genomic positions go in columns
# cut bins go in columns
# binding modes go in rows

#bins = c(rep(11, 20), rep(1, 45), rep(2, 35), rep(3, len), rep(4, 35), rep(5, 45), rep(11, 20),
#  rep(11, 20), rep(6, 45), rep(7, 35), rep(8, len), rep(9, 35), rep(10, 45), rep(11, 20))
# assume we have bins
#bins = c(rep(10, 20), rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20), rep(5, len), rep(6, 20), rep(7, 20), rep(8, 20), rep(9, 20), rep(10, 20),
#  rep(10, 20), rep(9, 20), rep(8, 20), rep(7, 20), rep(6, 20), rep(5, len), rep(4, 20), rep(3, 20), rep(2, 20), rep(1, 20), rep(10, 20))
stopifnot(ncol(bins1) == ncol(cuts1))
stopifnot(ncol(bins2) == ncol(cuts2))
stopifnot(nrow(bins1) == zero)
stopifnot(nrow(bins2) == zero)

# test if bin numbers are from the sequence 1...numbins[k] in each row
numbins1 = apply(bins1, 1, function(v) length(unique(v)))
numbins2 = apply(bins2, 1, function(v) length(unique(v)))
for (k in 1:zero)
{
  stopifnot(bins1[k, ] %in% 1:numbins1[k])
  stopifnot(bins2[k, ] %in% 1:numbins2[k])
}

bincuts1 = lapply(1:zero, function(k) sapply(1:numbins1[k], function(bin) apply(cuts1[, bins1[k, ] == bin, drop = F], 1, sum)))
bincuts2 = lapply(1:zero, function(k) sapply(1:numbins2[k], function(bin) apply(cuts2[, bins2[k, ] == bin, drop = F], 1, sum)))
binsizes1 = lapply(1:zero, function(k) sapply(1:numbins1[k], function(bin) sum(bins1[k, ] == bin)))
binsizes2 = lapply(1:zero, function(k) sapply(1:numbins2[k], function(bin) sum(bins2[k, ] == bin)))
print(binsizes1)
print(binsizes2)
params[[id]]$binsizes1 = binsizes1
params[[id]]$binsizes2 = binsizes2

PriorLik = matrix(0.01, nr, K + 1)
# FIXME initialization?
scsel = sumcuts1 + sumcuts2 > quantile(sumcuts1 + sumcuts2, probs = 0.9)
PriorLik[scsel, 1] = 100 # guess for monomer binding
for (k in 2:(K + 1))
  PriorLik[scsel & anno[, names(anno)[grep("score", names(anno))][k]] < quantile(anno[, names(anno)[grep("score", names(anno))][k]], probs = 0.1), k] = 100
print(summary(PriorLik))

PriorLogLik = log(PriorLik)

# assume we have Priors
add.intercept = T

cols = lapply(1:(zero - 1), function(k) colnames(anno) %in% Priors[[k]])

Beta = list()
for (k in 1:(K + 1))
{
  if (add.intercept)
    thislm = lm(PriorLogLik[, k] ~ ., data = as.data.frame(anno[, cols[[k]], drop = F]))
  else
    thislm = lm(PriorLogLik[, k] ~ . - 1, data = as.data.frame(anno[, cols[[k]], drop = F]))

  Beta[[k]] = thislm$coefficients
#  PriorLogLik[, k] = extanno[, extcols[[k]], drop = F] %*% Beta[[k]] # equivalent to predict(thislm)
#  print(Beta)
}

colsel = do.call(pmax, cols) == 1
if (add.intercept)
{
  extanno = cbind(1, as.matrix(anno[, colsel]))
  extcols = lapply(cols, function(v) c(T, v[colsel]))
}
else
{
  extanno = as.matrix(anno[, colsel])
  extcols = lapply(cols, function(v) v[colsel])
}

PriorLogLik = sapply(1:(zero - 1), function(k) extanno[, extcols[[k]], drop = F] %*% Beta[[k]])
PriorLik = exp(PriorLogLik)
PriorProbUnbound = 1 / (1 + apply(PriorLik, 1, sum))
PriorProb = cbind(PriorLik * PriorProbUnbound, PriorProbUnbound, deparse.level = 0)

# only for the first iteration
PostProb = PriorProb
print(summary(PostProb[, 1]))
OldPostProb = NA


for (iter in 1:100)
{
cat(paste("\nIteration: ", iter, "\n", sep = ""))

# now the maximization step

#tLogLC = sum(PostProb * PriorLogLik) - sum(log(1 + PriorLik))
#print(tLogLC)

# LogLC1: part of complete likelihood with Beta

print("LogLC1")

fnLogLC1 = function(BetaUnlisted)
{
  thisBeta = relist(BetaUnlisted, Beta)
  thisPriorLogLik = sapply(1:(zero - 1), function(k) extanno[, extcols[[k]], drop = F] %*% thisBeta[[k]])
  thisPriorLik = exp(thisPriorLogLik)

  return(sum(PostProb[, -zero] * thisPriorLogLik) - (K + 1) * sum(log(1 + apply(thisPriorLik, 1, sum))))
}

grLogLC1 = function(BetaUnlisted)
{
  thisBeta = relist(BetaUnlisted, Beta)
  thisPriorLogLik = sapply(1:(zero - 1), function(k) extanno[, extcols[[k]], drop = F] %*% thisBeta[[k]])
  thisPriorLik = exp(thisPriorLogLik)

  grLogLC1k = function(k)
  {
    s1 = crossprod(PostProb[, k], extanno[, extcols[[k]], drop = F])
    s2 = - crossprod(thisPriorLik[, k] / (1 + apply(thisPriorLik, 1, sum)), extanno[, extcols[[k]], drop = F])
    return(s1 + s2)
  }

  return(do.call(c, lapply(1:(zero - 1), grLogLC1k)))
}

opt = optim(unlist(Beta), fnLogLC1, grLogLC1, method = "BFGS", control = list(fnscale = -1))
Beta = relist(opt$par, Beta)
LogLC1 = opt$value
print(Beta)
print(LogLC1)

# LogLC2: part of complete likelihood with NegBinomP and NegBinomR

print("LogLC2")

sumPostProb = apply(PostProb, 2, sum)
sumPostProbSumcuts1 = apply(PostProb * sumcuts1, 2, sum)
sumPostProbSumcuts2 = apply(PostProb * sumcuts2, 2, sum)
sumPostProbLgammaSumcuts1 = apply(PostProb * lgammaSumcuts1, 2, sum)
sumPostProbLgammaSumcuts2 = apply(PostProb * lgammaSumcuts2, 2, sum)

fnLogLC21 = function(thisNegBinomR, k)
{
  if (!(thisNegBinomR > 0)) return(NA)
  thisNegBinomP = thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts1[k])

  s1 = sum(PostProb[, k] * lgamma(thisNegBinomR + sumcuts1))
  s2 = - sumPostProb[k] * lgamma(thisNegBinomR)
  s3 = - sumPostProbLgammaSumcuts1[k]
  s4 = sumPostProb[k] * thisNegBinomR * log(thisNegBinomP)
  s5 = sumPostProbSumcuts1[k] * log(1 - thisNegBinomP)

  return(s1 + s2 + s3 + s4 + s5)
}

fnLogLC22 = function(thisNegBinomR, k)
{
  if (!(thisNegBinomR > 0)) return(NA)
  thisNegBinomP = thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts2[k])

  s1 = sum(PostProb[, k] * lgamma(thisNegBinomR + sumcuts2))
  s2 = - sumPostProb[k] * lgamma(thisNegBinomR)
  s3 = - sumPostProbLgammaSumcuts2[k]
  s4 = sumPostProb[k] * thisNegBinomR * log(thisNegBinomP)
  s5 = sumPostProbSumcuts2[k] * log(1 - thisNegBinomP)

  return(s1 + s2 + s3 + s4 + s5)
}

grLogLC21 = function(thisNegBinomR, k)
{
  if (!(thisNegBinomR > 0)) return(NA)
  thisNegBinomP = thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts1[k])

  s1 = sum(PostProb[, k] * digamma(thisNegBinomR + sumcuts1))
  s2 = - sumPostProb[k] * digamma(thisNegBinomR)
  s4 = sumPostProb[k] * (1 - thisNegBinomP + log(thisNegBinomP))
  s5 = - sumPostProbSumcuts1[k] * thisNegBinomP / thisNegBinomR

  return(s1 + s2 + s4 + s5)
}

grLogLC22 = function(thisNegBinomR, k)
{
  if (!(thisNegBinomR > 0)) return(NA)
  thisNegBinomP = thisNegBinomR * sumPostProb[k] / (thisNegBinomR * sumPostProb[k] + sumPostProbSumcuts2[k])

  s1 = sum(PostProb[, k] * digamma(thisNegBinomR + sumcuts2))
  s2 = - sumPostProb[k] * digamma(thisNegBinomR)
  s4 = sumPostProb[k] * (1 - thisNegBinomP + log(thisNegBinomP))
  s5 = - sumPostProbSumcuts2[k] * thisNegBinomP / thisNegBinomR

  return(s1 + s2 + s4 + s5)
}

NegBinomR1 = rep(1, zero)
NegBinomR2 = rep(1, zero)
LogLC2 = 0
for (k in 1:zero)
{
  opt1 = optim(NegBinomR1[k], fnLogLC21, grLogLC21, k = k, method = "BFGS", control = list(fnscale = -1))
  opt2 = optim(NegBinomR2[k], fnLogLC22, grLogLC22, k = k, method = "BFGS", control = list(fnscale = -1))
  NegBinomR1[k] = opt1$par
  NegBinomR2[k] = opt2$par
  LogLC2 = LogLC2 + opt1$value + opt2$value
}
NegBinomP1 = NegBinomR1 * sumPostProb / (NegBinomR1 * sumPostProb + sumPostProbSumcuts1)
NegBinomP2 = NegBinomR2 * sumPostProb / (NegBinomR2 * sumPostProb + sumPostProbSumcuts2)

print(cbind(NegBinomR1, NegBinomR2, NegBinomP1, NegBinomP2))
print(LogLC2)
#LogLC2 = sum(sapply(1:zero, function(k) sum(PostProb[, k] * dnbinom(sumcuts, NegBinomR[k], NegBinomP[k], log = T))))
#print(LogLC2)

# LogLC3: part of complete likelihood with Lambda[[i]]

regularize = function(values, bins)
{
  rl = rle(bins)
  i = 0

  for (j in seq_len(length(rl$values) - 1))
  {
    i = i + rl$lengths[j]
    il = i - (rl$lengths[j] - 1) %/% 2
    ir = i + (rl$lengths[j + 1] + 1) %/% 2
    values[il:ir] = seq(from = values[il], to = values[ir], length.out = ir - il + 1)
  }

  return(values)
}

print("LogLC3")
Lambda1 = lapply(1:zero, function(k)
{
  lambda = as.vector(crossprod(PostProb[, k], bincuts1[[k]])) / binsizes1[[k]]
  lambda = lambda + sum(lambda * binsizes1[[k]]) / ncol(bins1) # FIXME estimator shrinkage
  lambda / sum(lambda * binsizes1[[k]])
})
print(Lambda1)

LambdaReg1 = t(sapply(1:zero, function(k) Lambda1[[k]][bins1[k, ]]))
for (k in 1:zero)
  LambdaReg1[k, ] = regularize(LambdaReg1[k, ], bins1[k, ])

Lambda2 = lapply(1:zero, function(k)
{
  lambda = as.vector(crossprod(PostProb[, k], bincuts2[[k]])) / binsizes2[[k]]
  lambda = lambda + sum(lambda * binsizes2[[k]]) / ncol(bins2) # FIXME estimator shrinkage
  lambda / sum(lambda * binsizes2[[k]])
})
print(Lambda2)

LambdaReg2 = t(sapply(1:zero, function(k) Lambda2[[k]][bins2[k, ]]))
for (k in 1:zero)
  LambdaReg2[k, ] = regularize(LambdaReg2[k, ], bins2[k, ])

#sapply(1:zero, function(k) print(sum(Lambda[[k]][bins[k, ]]) - 1))

LogLC3 = sum(sapply(1:zero, function(k)
{
  s1 = sum(PostProb[, k] * lgammaSumcuts1) + sum(PostProb[, k] %*% cuts1 * log(LambdaReg1[k, ])) - sum(PostProb[, k] * sumLgammaCuts1)
  s2 = sum(PostProb[, k] * lgammaSumcuts2) + sum(PostProb[, k] %*% cuts2 * log(LambdaReg2[k, ])) - sum(PostProb[, k] * sumLgammaCuts2)
  s1 + s2
}))
print(LogLC3)

#LogLC3alt = 0
#for (i in seq_len(nr))
#{
#  LogLC3alt = LogLC3alt + PostProb[, zero][i] * dmultinom(cuts[i, ], prob = Lambda[[zero]][bins[1, ]], log = T)
#  LogLC3alt = LogLC3alt + PostProb[i] * dmultinom(cuts[i, ], prob = Lambda[[1]][bins[1, ]], log = T)
#}
#print(LogLC3alt)

LogLC = LogLC1 + LogLC2 + LogLC3
print("LogLC total")
print(LogLC)


# now the maximization step

PriorLogLik = sapply(1:(zero - 1), function(k) extanno[, extcols[[k]], drop = F] %*% Beta[[k]])
PriorLik = exp(PriorLogLik)
PriorProbUnbound = 1 / (1 + apply(PriorLik, 1, sum))
PriorProb = cbind(PriorLik * PriorProbUnbound, PriorProbUnbound, deparse.level = 0)

LogCondProb11 = sapply(1:zero, function(k) dnbinom(sumcuts1, NegBinomR1[k], NegBinomP1[k], log = T))
LogCondProb12 = sapply(1:zero, function(k) dnbinom(sumcuts2, NegBinomR2[k], NegBinomP2[k], log = T))
LogCondProb21 = sapply(1:zero, function(k) lgammaSumcuts1 + cuts1 %*% log(LambdaReg1[k, ]) - sumLgammaCuts1)
LogCondProb22 = sapply(1:zero, function(k) lgammaSumcuts2 + cuts2 %*% log(LambdaReg2[k, ]) - sumLgammaCuts2)

# here we use the Bayes theorem
LogLikelihood = LogCondProb11 + LogCondProb12 + LogCondProb21 + LogCondProb22
# avoiding numerical overflow
JointProb = PriorProb * exp(LogLikelihood - apply(LogLikelihood, 1, max))
PostProb = JointProb / apply(JointProb, 1, sum)
print(summary(PostProb))

if (iter > 1)
{
  maxdiff = max(abs(OldPostProb - PostProb))
  print(maxdiff)
  if (maxdiff < 0.001) break
}
OldPostProb = PostProb

params[[id]]$Beta = Beta
params[[id]]$NegBinom = cbind(NegBinomR1, NegBinomR2, NegBinomP1, NegBinomP2)
params[[id]]$Lambda1 = Lambda1
params[[id]]$Lambda2 = Lambda2
params[[id]]$LambdaReg1 = LambdaReg1
params[[id]]$LambdaReg2 = LambdaReg2
}

### end of MOCCA code

for (col in signalcols)
{
  cat("\n")
  print(paste(col, K))
  id = paste(prefix, col, K)

#sum 1000
#  sel = order(anno$signalValue.AR.vehicle.siCTRL, decreasing = T)

  print(c(sum(anno[, col] > 0), 42))
#  anno$signalValue = ifelse(anno[, col] > 0, 1, 0)

  compare(1 - PostProb[, zero], col)
  prior[[id]] = PriorProb
  loglik[[id]] = LogLikelihood
  prob[[id]] = PostProb
  prof[[id]] = cbind(crossprod(PostProb, cuts1) / apply(PostProb, 2, sum), crossprod(PostProb, cuts2) / apply(PostProb, 2, sum)) # FIXME divide only once
}

}
}

print(warnings())
save(pred, predr, perf, perfr, auc, prior, prob, loglik, prof, params, file = paste(file.prefix, ".Rdata", sep = ""))
