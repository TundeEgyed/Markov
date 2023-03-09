#minta generálása
possibleValues = c("A", "B", "C", "D")
numberOfSamples = 1
sampleLength = 50
sample = matrix(ncol = numberOfSamples, nrow = sampleLength)

for (i in 1:numberOfSamples) {
  for (j in 1:2) {
    sample[j, i] = sample(possibleValues, size = 1, replace = TRUE)
  }
}

for (i in 1:numberOfSamples) {
  for (j in 3:sampleLength) {
    sample[j, i] = sample(c(rep(sample[j - 2, i], 2), rep(sample[j - 1, i], 3), possibleValues), size = 1, replace = TRUE)
  }
}

#maximum likelihood módszer
shift1 <- function(col) {
  return(rbind(head(col, -1), col[0:-1]))
}

shift2 <- function(col) {
  return(rbind(paste0(head(col, -2), head(col[-1], -1)), col[-1:-2]))
}

shift3 <- function(col) {
  return(rbind(paste0(head(col, -3), head(col[-1], -2), head(col[-1:-2], -1)), col[-1:-3]))
}

shift3 <- function(col) {
  return(rbind(paste0(head(col, -3), head(col[-1], -2), head(col[-1:-2], -1)), col[-1:-3]))
}

possibleTransitions = vector()
for (i in 1:length(possibleValues)){
  for (j in 1:length(possibleValues)){
    possibleTransitions = c(possibleTransitions, paste0(possibleValues[i], possibleValues[j]))
  }
}

calculateFrequencies <- function(order, sample) {
  frequencies <- apply(sample, 2, paste0("shift", order))
  frequencies <- data.frame(t(matrix(frequencies, nrow = 2)))
  names(frequencies) <- c("From", "To")

  return(table(frequencies))
}

calculateTransitionMatrix <- function(order, sample) {
  frequencies = calculateFrequencies(order, sample)
  print(frequencies)
  transitionMatrix = frequencies / rowSums(frequencies)
  
  return(transitionMatrix)
}

a = apply(frequencies, 1, sum)
plot(attr(a,"names"), a, type = "h")

calculateLoglikelihood <- function(order, sample) {
  loglikelihood = 0
  # loglikelihood = log(1 / length(possibleValues))
  frequencies = calculateFrequencies(order, sample)
  transitionMatrix = frequencies / rowSums(frequencies)
  for (i in 1:dim(transitionMatrix)[1]) {
    for (j in 1:length(possibleValues)) {
      logTransition = log(transitionMatrix)[i, j]
      logTransition[logTransition == -Inf] = 0
      loglikelihood = loglikelihood + logTransition * frequencies[i, j]
    }
  }

  return(loglikelihood)
}

loglikelihood = vector(length = 3)
for (i in 1:3) {
  loglikelihood[i] = calculateLoglikelihood(i, sample)
}
plot(loglikelihood, main = "Log likelihoods", xlab = "Rend", ylab = "Log likelihood", pch = 16, type = 'o', col = "blue", cex = 1.5)

likelihoodRatioTest <- function(k, m, sample) {
  T = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample))
  degreesOfFreedom = (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  print(T)
  print(degreesOfFreedom)
  print(pchisq(T, degreesOfFreedom))
  p = 1 - pchisq(T, degreesOfFreedom)
  
  return(p)
}

likelihoodPValue = matrix(ncol = 2, nrow = 0)
for (i in 1:(3 - 1)) {
  for (j in (i + 1):3) {
    likelihoodPValue = rbind(likelihoodPValue, c(paste(i, j), likelihoodRatioTest(i, j, sample)))
  }
}
likelihoodPValue

# Akaike
calculateAkaike <- function(k, m, sample) {
  AIC = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  
  return(AIC)
}

AIC = vector(length = 3)
for (i in 1:3) {
  AIC[i] = calculateAkaike(i, 3, sample)
}
plot(AIC, main = "AIC with test model m = 3", pch = 16, type = 'o', col = "blue")

# Bayes
calculateBayes <- function(k, m, sample) {
  BIC = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1) * log(sampleLength * numberOfSamples)
  
  return(BIC)
}

BIC = vector(length = 3)
for (i in 1:3) {
  BIC[i] = calculateBayes(i, 3, sample)
}
plot(BIC, main = "BIC with test model m = 3", pch = 16, type = 'o', col = "blue")

# cross validation
calculateAverageRank <- function(order, sample, sample2) {
  averageRank = 0
  transitionMatrix = calculateTransitionMatrix(order, sample)
  print(transitionMatrix)
  frequencies = calculateFrequencies(order, sample2)
  for (k in 1:length(possibleValues)^order) {
    for (i in 1:length(possibleValues)) {
      for (j in 1:length(possibleValues)) {
        if (transitionMatrix[k, i] == sort(transitionMatrix[k,], decreasing = TRUE)[j]) {
          averageRank = averageRank + frequencies[k, i] * j
        }
      }
    }
  }
  averageRank = averageRank / sum(frequencies)
  
  return(averageRank)
}

sample1Length = 50
sample1 = sample[1:sample1Length,]
sample2 = sample[(sample1Length + 1):length(sample),]
dim(sample1) = c(sample1Length, numberOfSamples)
dim(sample2) = c(length(sample) - sample1Length, numberOfSamples)
calculateAverageRank(2, sample1, sample2)

#két-lépéses visszatérés
shiftTwoStep1 <- function(col) {
  return(rbind(head(col, -1), col[0:-1]))
}

shiftTwoStep2 <- function(col) {
  return(rbind(paste0(head(col, -2), head(col[-1], -1)), paste0(head(col[-1], -1), col[-1:-2])))
}

shiftTwoStep3 <- function(col) {
  return(rbind(paste0(head(col, -3), head(col[-1], -2), head(col[-1:-2], -1)), paste0(head(col[-1], -2), head(col[-1:-2], -1), col[-1:-3])))
}

calculateTransitions <- function(order){
  transitions <- apply(sample, 2, paste0("shiftTwoStep", order))
  transitions <- data.frame(t(matrix(transitions, nrow = 2)))
  
  return(transitions)
}

calculateFrequencies2 <- function(order){
  transitions <- calculateTransitions(order)
  names(transitions) <- c("From", "To")
  frequencies <- table(transitions)
  
  return(frequencies)
}

calculateTransitionMatrix2 <- function(order) {
  frequencies = calculateFrequencies2(order)
  transitionMatrix = frequencies / rowSums(frequencies)
  
  return(transitionMatrix)
}

calculateTwoStep <- function(order) {
  transitionMatrix = calculateTransitionMatrix2(order)
  normalisedEigenVector = eigen(t(transitionMatrix))$vectors[,1] / sum(eigen(t(transitionMatrix))$vectors[,1])
  
  transitionMatrixSquared = transitionMatrix %*% transitionMatrix
  possibleTransitions = unique(calculateTransitions(order)[,1])
  A = vector()
  for (i in sort(possibleTransitions, decreasing = FALSE)) {
    A[i] = sum(transitionMatrixSquared[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 1, 1) == substr(i, 1, 1)]])
  }
  
  return(A %*% normalisedEigenVector)
}

#relativeFrequencies
relativeFrequencies = (function() {
  possibleTransitions = unique(calculateTransitions(2)[,1])
  relativeFrequencies = 0
  frequencies = calculateFrequencies2(2)
  for (i in possibleTransitions) {
    relativeFrequencies = relativeFrequencies + sum((frequencies[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 2, 2) == substr(i, 1, 1)]]))
  }
  
  return(relativeFrequencies / sum(frequencies))
})()

#entropy
calculateEntropyRate <- function(order) {
  transitionMatrix = calculateTransitionMatrix2(order)
  normalisedEigenVector = eigen(t(transitionMatrix))$vectors[,1] / sum(eigen(t(transitionMatrix))$vectors[,1])
  
  A = vector()
  for (i in 1:dim(transitionMatrix)[1]) {
    logTransitionMatrix = log(transitionMatrix)[i]
    logTransitionMatrix[logTransitionMatrix == -Inf] = 0
    A[i] = sum((transitionMatrix * logTransitionMatrix)[i]) * normalisedEigenVector[i]
  }
  
  return(-sum(A))
}
