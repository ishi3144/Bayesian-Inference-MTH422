rhathat = function(type, yObs, classifyParam, q, sd, xObs, mu, sigma) {
  nA = 10
  a = rep(0, nA)
  for (i in 1:nA) {
    xGenerated = rnorm(length(xObs), mean = mu, sd = sigma)
    a[i] = classify(type, xObs, xGenerated, classifyParam)
  }
  a = mean(a)
  
  tau = rgeom(1, q) + 10
  r = 1
  
  if (tau >= 1) {
    for (n in 1:tau) {
      term = (-1)^n / factorial(n)
      if (n > 10) term = term / ((1.0 - q)^(n - 10))
      term = term * H(n, (yObs - a) / sd)
      
      for (i in 1:n) {
        xGenerated = rnorm(length(xObs), mean = mu, sd = sigma)
        y0 = classify(type, xObs, xGenerated, classifyParam)
        term = term * (a - y0) / sd
      }
      
      r = r + term  # ✅ now correctly inside loop
    }
  }
  
  r = r * exp(-0.5 * ((yObs - a) / sd)^2)
  if (is.nan(r) || is.infinite(r)) r = 0
  return(r)
}
