vv = function(type,yObs,classifyParam,q,sd,xObs,mu,sigma) {
    nA = 10
    a = rep(0,nA)
    for (i in c(1:nA)) {
        xGenerated = rnorm(length(xObs),mean = mu,sd = sigma)
        a[i] = classify(type,xObs,xGenerated,classifyParam)
    }
    muA = mean(a)
    varA = var(a)

    a = c()
    var = c()
    bruttovarians = c()

    for (aa in seq(from=0,to=1,by=0.01)) {
        tau = 10

        vv = 0
        mm = 0
        
        for (n in c(0:tau)) {
            factor = (-1)^n / factorial(n)
            factor = factor * H(n,(yObs - aa) / sd)

            vv = vv + factor^2 * (varA / (sd^2))^n
            mm = mm + factor * ((aa - muA) / sd)^n
        }
        a = c(a,aa)
        var = c(var,vv)
        mean10 = mm

        
        tau = 40
        mm = 0
        
        for (n in c(0:tau)) {
            factor = (-1)^n / factorial(n)
            factor = factor * H(n,(yObs - aa) / sd)

            mm = mm + factor * ((aa - muA) / sd)^n
        }
        mean20 = mm

        bruttovarians = c(bruttovarians,log((mean10-mean20)^2 + vv) - 0.5 * ((yObs-aa)/sd)^2)
    }


    
    vv = list(a = a,var = var,bruttovarians = bruttovarians,mean,muA = muA,varA = varA)
}

        
