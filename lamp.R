#LAMP
order = 2
possibleValues = sort(unique(sample))
n = length(possibleValues)
psi = rep(1 / order, order)
a = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))
for (i in 1:order) {
  a[,,i] = matrix(rep(1 / n, n * n), nrow = n)
}
psi2 = vector(length = order)
a2 = array(dim = c(n, n, order), dimname = list(possibleValues, possibleValues))

for (j in 1:100) {
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
  psi = psi2 / sum(psi2)
  for (mu in 1:order) {
    a[,,mu] = a2[,,mu] / rowSums(a2[,,mu])
  }
}
