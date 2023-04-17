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

#maximum likelihood
loglikelihood = vector(length = 5)
for (i in 1:5) {
  loglikelihood[i] = calculateLoglikelihood(i, sample)
}
plot(loglikelihood, main = "Log likelihoods", xlab = "Rend", ylab = "Log likelihood", pch = 16, type = 'o', col = "blue", cex = 1.5)

likelihoodPValue = matrix(ncol = 2, nrow = 0)
for (i in 1:(5 - 1)) {
  for (j in (i + 1):5) {
    likelihoodPValue = rbind(likelihoodPValue, c(paste(i, j), likelihoodRatioTest(i, j, sample)))
  }
}
likelihoodPValue

#Akaike
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

#cross validation
trainingSetLength = 500
trainigSet = sample[1:trainingSetLength,]
validationSet = sample[(trainingSetLength + 1):length(sample),]
dim(trainigSet) = c(trainingSetLength, numberOfSamples)
dim(validationSet) = c(length(sample) - trainingSetLength, numberOfSamples)
calculateAverageRank(2, trainigSet, validationSet)

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