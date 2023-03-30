sample = read.csv("/home/tunde/Dokumentumok/Önálló projekt/SolarFlare_list_GOES.csv", sep = ",")$GOES.Class.2
# sample = read.csv("SolarFlare_list_GOES.csv", sep = ",")$GOES.Class

dim(sample) = c(length(sample), 1)
possibleValues = unique(sample)
numberOfSamples = 1
sampleLength = length(sample)

#likelihood
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

#AIC
m = 5
AIC = vector()
for (i in 1:m) {
  AIC[i] = calculateAkaike(i, m, sample)
}
plot(AIC, main = paste("AIC with test model m = ", m), pch = 16, type = 'o', col = "blue")

#BIC
BIC = vector()
for (i in 1:m) {
  BIC[i] = calculateBayes(i, m, sample)
}
plot(BIC, main = paste("BIC with test model m = ", m), pch = 16, type = 'o', col = "blue")

#cross validation
trainingSetLength = 5000
trainigSet = sample[1:trainingSetLength,]
validationSet = sample[(trainingSetLength + 1):length(sample),]
calculateAverageRank(1, trainigSet, validationSet)

# két lépéses visszatérés
relativeFrequencies = (function() {
  possibleTransitions = unique(calculateTransitions2(2)[,1])
  relativeFrequencies = 0
  frequencies = calculateFrequencies2(2)
  for (i in possibleTransitions) {
    relativeFrequencies = relativeFrequencies + sum((frequencies[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 2, 2) == substr(i, 1, 1)]]))
  }
  
  return(relativeFrequencies / sum(frequencies))
})()

calculateTwoStep(1)

#entropy
calculateEntropyRate(1)
