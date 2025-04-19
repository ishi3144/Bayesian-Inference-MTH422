rhathat = function(type, yObs, classifyParam, q, sd, xObs, mu, sigma) {
  nA = 10
  a = numeric(nA)
  for (i in 1:nA) {
    xGenerated = rnorm(length(xObs), mean = mu, sd = sigma)
    a[i] = classify(type, xObs, xGenerated, classifyParam)
  }
  rhat = mean(a)
  
  tau = rgeom(1, q) + 10
  result = 0
  
  for (n in 0:tau) {
    coeff = if (n == 0) 1 else (-1)^n / factorial(n)
    if (n > 10) coeff = coeff / ((1 - q)^(n - 10))
    
    Hn = H(n, (yObs - rhat) / sd)
    prod_term = 1
    
    if (n >= 1) {
      for (j in 1:n) {
        xGenerated = rnorm(length(xObs), mean = mu, sd = sigma)
        s_val = classify(type, xObs, xGenerated, classifyParam)
        prod_term = prod_term * (s_val - rhat) / sd
      }
    }
    
    result = result + coeff * Hn * prod_term
  }
  
  r = result * dnorm((yObs - rhat) / sd) / sd
  
  if (is.nan(r) || is.infinite(r)) r = 0
  return(r)
}
