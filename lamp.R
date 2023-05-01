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

order = 1
possibleValues = unique(sample)
n = length(possibleValues)
psi = rep(1 / order, order)
a = array(dim = c(n, n, order))
for (i in 1:order) {
  a[,,i] = matrix(rep(1 / n, n * n), nrow = n)
}

for (i in 1:100) {
  for (j in 1:order) {
    psi[j] = psi[j] * sum(calculateFrequencies(j, sample) * a[,,j])
  }
  psi = psi / sum(psi)
  
  for (j in 1:order) {
    a[,,j] = psi[j] * calculateFrequencies(j, sample) * a[,,j]
    a[,,j] = a[,,j] / rowSums(a[,,j])
  }
}
