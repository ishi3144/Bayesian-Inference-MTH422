rhathat = function(type,yObs,classifyParam,q,sd,xObs,mu,sigma) {
  nA = 10
  a = rep(0,nA)
  for (i in c(1:nA)) {
    xGenerated = rnorm(length(xObs),mean = mu,sd = sigma)
    a[i] = classify(type,xObs,xGenerated,classifyParam)
  }
  a = mean(a)
  
  tau = rgeom(1,q) + 10
  
  r = 1
  if (tau >= 1) {
    for (n in c(1:tau)) {
      term = (-1)^n / factorial(n)
      if (n > 10) term = term / ((1.0 - q)^(n-10))
      term = term * H(n,(yObs - a) / sd)
      
      for (i in c(1:n)) {
        y0 = 0
        for (k in c(1:1)) {
          xGenerated = rnorm(length(xObs),mean = mu,sd = sigma)
          y0 = y0 + classify(type,xObs,xGenerated,classifyParam)
        }
        y0 = y0 / 1
        
        term = term * (a - y0) / sd
      }
    }
    
    r = r + term
  }
  
  r = r * exp( - 0.5 * ((yObs - a) / sd)^2)
  
  r
}
