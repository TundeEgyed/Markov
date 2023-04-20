#Transition matrix
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

#maximum likelihood method
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

likelihoodRatioTest <- function(k, m, sample) {
  T = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample))
  degreesOfFreedom = (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  print(T)
  print(degreesOfFreedom)
  print(pchisq(T, degreesOfFreedom))
  p = 1 - pchisq(T, degreesOfFreedom)
  
  return(p)
}

# Akaike
calculateAkaike <- function(k, m, sample) {
  AIC = -2 * (calculateLoglikelihood(k, sample) - calculateLoglikelihood(m, sample)) - 2 * (length(possibleValues) ^ m - length(possibleValues) ^ k)*(length(possibleValues) - 1)
  
  return(AIC)
}

# cross validation
calculateAverageRank <- function(order, trainingSet, validationSet) {
  possibleTransitionsTraning = sort(unique(calculateTransitions(order, trainingSet)[,1]))
  possibleTransitionsValid = sort(unique(calculateTransitions(order, validationSet)[,1]))
  rankMatrix = matrix(nrow = length(possibleTransitionsValid), ncol = length(possibleValues))
  frequenciesOfTraningSet = calculateFrequencies(order, trainingSet)
  for (i in 1:length(possibleTransitionsValid)) {
    for (j in 1:length(possibleTransitionsTraning)) {
      if (possibleTransitionsValid[i] == possibleTransitionsTraning[j]) {
        for (k in 1:length(possibleValues)) {
          for (l in 1:length(possibleValues))
            if (frequenciesOfTraningSet[j, k] == sort(frequenciesOfTraningSet[j,])[l]) {
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
