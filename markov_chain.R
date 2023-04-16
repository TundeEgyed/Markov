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
shift <- function(order, sample) {
  a = t(head(sample, -order))
  if (order != 1) {
    for (i in 1:(order - 1)) {
      a = paste0(a, head(sample[-1:-i], -(order - i)))
    }
  }
  
  return(rbind(a, sample[0:-order]))
}

calculateTransitions <- function(order, sample) {
  transitions = shift(order, sample)
  transitions = data.frame(t(matrix(transitions, nrow = 2)))
  
  return(transitions)
}

calculateFrequencies <- function(order, sample) {
  frequencies = calculateTransitions(order, sample)
  names(frequencies) <- c("From", "To")
  
  return(table(frequencies))
}

calculateTransitionMatrix <- function(order, sample) {
  frequencies = calculateFrequencies(order, sample)
  transitionMatrix = frequencies / rowSums(frequencies)
  
  return(transitionMatrix)
}

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

loglikelihood = vector(length = 5)
for (i in 1:5) {
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
for (i in 1:(5 - 1)) {
  for (j in (i + 1):5) {
    likelihoodPValue = rbind(likelihoodPValue, c(paste(i, j), likelihoodRatioTest(i, j, sample)))
  }
}
likelihoodPValue

# Akaike
calculateAkaike <- function(k, m, sample) {
  AIC = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  
  return(AIC)
}

m = 5
AIC = vector()
for (i in 1:m) {
  AIC[i] = calculateAkaike(i, m, sample)
}
plot(AIC, main = paste("AIC with test model m = ", m), pch = 16, type = 'o', col = "blue")

# Bayes
calculateBayes <- function(k, m, sample) {
  BIC = - 2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1) * log(sampleLength * numberOfSamples)
  
  return(BIC)
}

BIC = vector()
for (i in 1:m) {
  BIC[i] = calculateBayes(i, m, sample)
}
plot(BIC, main = paste("BIC with test model m = ", m), pch = 16, type = 'o', col = "blue")

# cross validation
calculateAverageRank <- function(order, trainingSet, validationSet) {
  possibleTransitionsTraning = sort(unique(calculateTransitions(order, trainingSet)[,1]))
  possibleTransitionsValid = sort(unique(calculateTransitions(order, validationSet)[,1]))
  rankMatrix = matrix(nrow = length(possibleTransitionsValid), ncol = length(possibleValues))
  transitionMatrix = calculateTransitionMatrix(order, trainingSet)
  for (i in 1:length(possibleTransitionsValid)) {
    for (j in 1:length(possibleTransitionsTraning)) {
      if (possibleTransitionsValid[i] == possibleTransitionsTraning[j]) {
        for (k in 1:length(possibleValues)) {
          for (l in 1:length(possibleValues))
            if (transitionMatrix[j, k] == sort(transitionMatrix[j,], decreasing = TRUE)[l]) {
              rankMatrix[i, k] = length(possibleValues) - l + 1
              break
            }
        }
      }
    }
  }
  rankMatrix[is.na(rankMatrix)] = length(possibleValues)
  frequencies = calculateFrequencies(order, validationSet)
  averageRank = sum(rankMatrix * frequencies) / sum(frequencies)
  
  return(averageRank)
}

trainingSetLength = 500
trainigSet = sample[1:trainingSetLength,]
validationSet = sample[(trainingSetLength + 1):length(sample),]
dim(trainigSet) = c(trainingSetLength, numberOfSamples)
dim(validationSet) = c(length(sample) - trainingSetLength, numberOfSamples)
calculateAverageRank(2, trainigSet, validationSet)

#két-lépéses visszatérés
shift2 <- function(order, sample) {
  a = head(sample, -order)
  if (order != 1) {
    for (i in 1:(order - 1)) {
      a = paste0(a, head(sample[-1:-i], -(order - i)))
    }
  }
  
  return(rbind(t(head(a, -1)), t(a[-1])))
}

calculateTransitions2 <- function(order){
  transitions = shift2(order, sample)
  transitions = data.frame(t(matrix(transitions, nrow = 2)))
  
  return(transitions)
}

calculateFrequencies2 <- function(order){
  transitions <- calculateTransitions2(order)
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
  possibleTransitions = unique(calculateTransitions2(order)[,1])
  twoStepProbability = vector()
  for (i in sort(possibleTransitions, decreasing = FALSE)) {
    twoStepProbability[i] = sum(transitionMatrixSquared[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 1, 1) == substr(i, 1, 1)]])
  }
  
  return(Re(twoStepProbability %*% normalisedEigenVector))
}

#relativeFrequencies
relativeFrequencies = (function() {
  possibleTransitions = unique(calculateTransitions2(2)[,1])
  relativeFrequencies = 0
  frequencies = calculateFrequencies2(2)
  for (i in possibleTransitions) {
    relativeFrequencies = relativeFrequencies + sum((frequencies[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 2, 2) == substr(i, 1, 1)]]))
  }
  
  return(Re(relativeFrequencies / sum(frequencies)))
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
  
  return(-Re(sum(A)))
}
