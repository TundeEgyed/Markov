#minta generálása
possibleValues = c("A", "B", "C", "D")
numberOfSamples = 1
sampleLength = 500
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

calculateFrequencies <- function(order, sample) {
  frequencies <- apply(sample, 2, paste0("shift", order))
  frequencies <- data.frame(t(matrix(frequencies, nrow = 2)))
  names(frequencies) <- c("From", "To")

  return(table(frequencies))
}

calculateTransitionMatrix <- function(order, sample) {
  frequencies = calculateFrequencies(order, sample)
  transitionMatrix = frequencies / rowSums(frequencies)
  
  return(transitionMatrix)
}

calculateAverageRank <- function(order, sample) {
  averageRank = 0
  transitionMatrix = calculateTransitionMatrix(order, sample)
  frequencies = calculateFrequencies(order, sample)
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

calculateAverageRank(1, ts(sample, start = 1, end = 250))

calculateLoglikelihood <- function(order) {
  loglikelihood = 0
  # loglikelihood = log(1 / length(possibleValues))
  frequencies = calculateFrequencies(order)
  transitionMatrix = frequencies / rowSums(frequencies)
  print(transitionMatrix)
  for (i in 1:dim(transitionMatrix)[1]) {
    for (j in 1:length(possibleValues)) {
      loglikelihood = loglikelihood + log(transitionMatrix[i, j] ^ frequencies[i, j])
    }
  }

  return(loglikelihood)
}

loglikelihood = vector(length = 3)
for (i in 1:3) {
  loglikelihood[i] = calculateLoglikelihood(i)
}
plot(loglikelihood, main = "Log likelihoods", pch = 16, type = 'o', col = "blue", cex = 1.5)

likelihoodRatioTest <- function(k, m) {
  (T = -2 * (calculateLoglikelihood(k) - calculateLoglikelihood(m)))
  degreesOfFreedom = (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  p = 1 - pchisq(T, degreesOfFreedom)
  
  return(p)
}

likelihoodPValue = matrix(ncol = 2, nrow = 0)
for (i in 1:(3 - 1)) {
  for (j in (i + 1):3) {
    likelihoodPValue = rbind(likelihoodPValue, c(paste(i, j), likelihoodRatioTest(i, j)))
  }
}
likelihoodPValue

# Akaike
calculateAkaike <- function(k, m) {
  AIC = -2 * (calculateLoglikelihood(k) - calculateLoglikelihood(m)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  
  return(AIC)
}

AIC = vector(length = 3)
for (i in 1:3) {
  AIC[i] = calculateAkaike(i, 3)
}
plot(AIC, main = "AIC with test model m = 3", pch = 16, type = 'o', col = "blue")

# Bayes
calculateBayes <- function(k, m) {
  BIC = -2 * (calculateLoglikelihood(k) - calculateLoglikelihood(m)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1) * log(sampleLength * numberOfSamples)
  
  return(BIC)
}

BIC = vector(length = 3)
for (i in 1:3) {
  BIC[i] = calculateBayes(i, 3)
}
plot(BIC, main = "BIC with test model m = 3", pch = 16, type = 'o', col = "blue")

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

calculateTwoStep <- function(order) {
  transitions <- apply(sample, 2, paste0("shiftTwoStep", order))
  transitions <- data.frame(t(matrix(transitions, nrow = 2)))
  names(transitions) <- c("From", "To")
  frequencies <- table(transitions)
  transitionMatrix <- frequencies / rowSums(frequencies)
  
  NormalisedEigenVector = eigen(t(transitionMatrix))$vectors[,1] / sum(eigen(t(transitionMatrix))$vectors[,1])
  
  transitionMatrixSquared = transitionMatrix %*% transitionMatrix
  A = vector(length = length(possibleValues) ^ order)
  for (i in 1:length(possibleValues) ^ order) {
    A[i] = sum(transitionMatrixSquared[i, seq(i %% length(possibleValues), length(possibleValues) ^ order, by = length(possibleValues))])
  }
  return(A %*% NormalisedEigenVector)
}

transitions <- apply(sample, 2, shiftTwoStep2)
transitions <- data.frame(t(matrix(transitions, nrow = 2)))
names(transitions) <- c("From", "To")
frequencies <- table(transitions)
frequenciesTwoStep = 0
for (i in 1:length(possibleValues)^2) {
  frequenciesTwoStep = frequenciesTwoStep + frequencies[i, ceiling(i / length(possibleValues)) + ((i %% length(possibleValues)) + length(possibleValues) * (i %% length(possibleValues) == 0) - 1) * length(possibleValues)]
}
relFrequenciesTwoStep = frequenciesTwoStep / sum(frequencies)
relFrequenciesTwoStep
