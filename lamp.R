#LAMP
order = 2
possibleValues = sort(unique(sample))
n = length(possibleValues)
psi = runif(order, 0, 1)
psi = psi / sum(psi)
#psi = rep(1 / order, order)
a = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))
for (i in 1:order) {
  a[,,i] = matrix(runif(n * n, 0, 1), nrow = n)
  a[,,i] = a[,,i] / rowSums(a[,,i])
}
#for (i in 1:order) {
#  a[,,i] = matrix(rep(1 / n, n * n), nrow = n)
#}
psi2 = vector(length = order)
a2 = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))
loglikelihood = rep(0, 10)

for (j in 1:10) {
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
a = round(a, 3)
