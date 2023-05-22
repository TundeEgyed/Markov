#LAMP
shift <- function(mu, sample) {
  a = t(head(sample, -mu))
  return(rbind(a, sample[0:-mu]))
}

calculateTransitions <- function(mu, sample) {
  transitions = shift(mu, sample)
  transitions = data.frame(t(matrix(transitions, nrow = 2)))
  
  return(transitions)
}

calculateFrequencies <- function(mu, sample) {
  frequencies = calculateTransitions(mu, sample)
  names(frequencies) <- c("From", "To")
  
  return(table(frequencies))
}

calculateRealtiveFrequencies <- function(mu, sample) {
  frequencies = calculateFrequencies(mu, sample)
  realtiveFrequencies = frequencies / rowSums(frequencies)
  realtiveFrequencies[realtiveFrequencies == 0] = 0.001
  return(realtiveFrequencies)
}

numberOfIterations = 100
possibleValues = sort(unique(sample))
n = length(possibleValues)
sampleLength = length(sample)

calculateLamp <- function(order) {
  #véletlen indítás
  #psi = runif(order, 0, 1)
  #psi = psi / sum(psi)
  psi = rep(1 / order, order) # egyenletes
  a = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))
  for (i in 1:order) {
    a[,,i] = calculateRealtiveFrequencies(mu, sample)
  }
  for (i in 1:order) {
    ###relativ gyak
    #a[,,i] = calculateRealtiveFrequencies(mu, sample)
    ###véletlen
    a[,,i] = matrix(runif(n * n, 0, 1), nrow = n)
    a[,,i] = a[,,i] / rowSums(a[,,i])
    ###egyenletes
    #a[,,i] = matrix(rep(1 / n, n * n), nrow = n)
  }
  psi2 = vector(length = order)
  a2 = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))
  loglikelihood = rep(0, numberOfIterations)
  
  for (j in 1:numberOfIterations) {
    a2[1 == 1] = 0
    for (mu in 1:order) {
      psi2[mu] = 0
      for (i in (order + 1):length(sample)) {
        denominator = 0
        for (nu in 1:order) {
          denominator = denominator + (psi[nu] * a[sample[i - nu], sample[i], nu])
        }
        psi2[mu] = psi2[mu] + (psi[mu] * a[sample[i - mu], sample[i], mu] / denominator)
        a2[sample[i - mu], sample[i], mu] = a2[sample[i - mu], sample[i], mu] + (psi[mu] * a[sample[i - mu], sample[i], mu] / denominator)
      }
    }
    for (i in (order + 1):length(sample)) {
      p = 0
      for (mu in 1:order) {
        p = p + psi[mu] * a[sample[i - mu], sample[i], mu] 
      }
      loglikelihood[j] = loglikelihood[j] + log(p)
    }
    psi = psi2 / sum(psi2)
    for (mu in 1:order) {
      a[,,mu] = a2[,,mu] / rowSums(a2[,,mu])
    }
  }
  #a = round(a, 3)
  
  #return(psi)
  return(loglikelihood[numberOfIterations])
}

maximumOrder = 5
loglikelihood = vector(length = maximumOrder)
#for (i in 1:maximumOrder) {
#  loglikelihood[i] = calculateLamp(i)
#}
#plot(loglikelihood, main = "Log likelihood értékek", xlab = "Rend", ylab = "Log likelihood", pch = 16, type = 'o', col = "blue", cex = 1.5)


calculateAkaikeLamp <- function(order) {
  AIC = 2 * (order * n * (n - 1) + order - 1) - 2 * calculateLamp(order)
  
  return(AIC)
}
#calculateAkaikeLamp(1) - calculateAkaikeLamp(3)
#calculateAkaikeLamp(2) - calculateAkaikeLamp(3)

calculateBayesLamp <- function(order) {
  BIC = log(sampleLength) * (order * n * (n - 1) + order - 1) - 2 * calculateLamp(order)
  
  return(BIC)
}
