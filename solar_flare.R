library(gridExtra)
library(ggplot2)

sample = read.csv("/home/tunde/Dokumentumok/Önálló projekt/SolarFlare_list_GOES.csv", sep = ",")$GOES.Class.3
# sample = read.csv("SolarFlare_list_GOES.csv", sep = ",")$GOES.Class
df = data.frame(sample)
ggplot(df, aes(sample)) + geom_bar(aes(y=..count.., fill=sample))  + labs(title = "Frequencies of Categories", x = "Categories", y = "Frequency") + theme(legend.position = 'none')

dim(sample) = c(length(sample), 1)
possibleValues = unique(sample)
numberOfSamples = 1
sampleLength = length(sample)
maximumOrder = 3

#likelihood
loglikelihood = vector(length = maximumOrder)
for (i in 1:maximumOrder) {
  loglikelihood[i] = calculateLoglikelihood(i, sample)
}
plot(loglikelihood, main = "Log likelihoods", xlab = "Order", ylab = "Log likelihood", pch = 16, type = 'o', col = "blue", cex = 1.5)

likelihoodPValue = matrix(ncol = 2, nrow = 0)
for (i in 1:(maximumOrder - 1)) {
  for (j in (i + 1):maximumOrder) {
    likelihoodPValue = rbind(likelihoodPValue, c(paste(i, j), likelihoodRatioTest(i, j, sample)))
  }
}

#AIC
AIC = vector()
for (i in 1:maximumOrder) {
  AIC[i] = calculateAkaike(i, maximumOrder, sample)
}
plot(AIC, main = paste("AIC with test model m = ", maximumOrder), pch = 16, type = 'o', col = "blue")

#BIC
BIC = vector()
for (i in 1:maximumOrder) {
  BIC[i] = calculateBayes(i, maximumOrder, sample)
}
plot(BIC, main = paste("BIC with test model m = ", maximumOrder), pch = 16, type = 'o', col = "blue")

#cross validation
trainingSetLength = 10000
trainigSet = sample[1:trainingSetLength,]
validationSet = sample[(trainingSetLength + 1):length(sample),]
calculateAverageRank(1, trainigSet, validationSet)

averageRank = vector()
for (i in 1:maximumOrder) {
  averageRank[i] = calculateAverageRank(i, trainigSet, validationSet)
}
plot(averageRank, main = "Cross Validation", pch = 16, type = 'o', col = "blue")


# két lépéses visszatérés
relativeFrequencies = (function() {
  possibleTransitions = unique(calculateTransitions2(2)[,1])
  relativeFrequencies = 0
  frequencies = calculateFrequencies2(2)
  for (i in possibleTransitions) {
    relativeFrequencies = relativeFrequencies + sum((frequencies[possibleTransitions[possibleTransitions == i], possibleTransitions[substr(possibleTransitions, 2, 2) == substr(i, 1, 1)]]))
  }
  
  return(Re(relativeFrequencies / sum(frequencies)))
})()

twoStepReturn = vector()
twoStepReturnDiff = vector()
for (i in 1:maximumOrder) {
  twoStepReturn[i] = calculateTwoStep(i)
  twoStepReturnDiff[i] = abs(twoStepReturn[i] - relativeFrequencies)
}
plot.new()
grid.table(cbind(1:maximumOrder, round(twoStepReturn, 2)), cols = c("Order", "Two-Step Return"))

#entropy
entropyRate = vector()
for (i in 1:maximumOrder) {
  entropyRate[i] = calculateEntropyRate(i)
}
plot.new()
grid.table(cbind(1:maximumOrder, round(entropyRate, 2)), cols = c("Order", "Entropy Rate"))
