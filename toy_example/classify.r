classify = function(type,xObs,xGenerated,classifyParam=NULL,method="logistic",nfolds=5,shuffle=FALSE) 
{
    n = length(xObs)
    
    if(shuffle) 
    {
        xObs = xObs[sample(n)]
        xGenerated = xGenerated[sample(n)]
    }

    predictions = rep(NA,n)

    idx.test.start = round(seq(1,n+1,length.out=nfolds+1))
    for(i in 1:nfolds) 
    {
        idx.test = idx.test.start[i]:(idx.test.start[i+1]-1)
        n.test = length(idx.test)
        idx.train = (1:n)[-idx.test]
        n.train = length(idx.train)
        x.train = c( xObs[idx.train], xGenerated[idx.train] )
        y.train = c( rep(0,n.train), rep(1,n.train) )
        if(type == 0) 
        {
            x.test = xObs[idx.test]
            y.test = rep(0,n.test)
        } 
        else 
        {
            x.test = xGenerated[idx.test]
            y.test = rep(1,n.test)
        }
        if(method == "logistic") 
        {
            data.train = data.frame(x = x.train, x2 = x.train^2)
            fit = glm(y.train ~ ., family = binomial("logit"), data = data.train)
            data.test = data.frame(x = x.test, x2 = x.test^2)
            predictions[idx.test] = predict(object = fit, type = "response", newdata = data.test)
        } 
        else if(method == "knn") 
        {
            kk = classifyParam$kk
            if(type == 0) 
            {
                res = knn(train = matrix(x.train,ncol=1), test = matrix(x.test,ncol=1), cl = y.train, k = kk, prob=TRUE)
                prob = (res == 0) * (1 - attr(res,"prob")) + (res == 1) * attr(res,"prob")
                predictions[idx.test] = 1.0 - prob
            } 
            else 
            {
                res = knn(train = matrix(x.train,ncol=1), test = matrix(x.test,ncol=1), cl = y.train, k = kk, prob=TRUE)
                prob = (res == 0) * (1 - attr(res,"prob")) + (res == 1) * attr(res,"prob")
                predictions[idx.test] = prob 
            }
        }
    }
    return(mean(predictions))
}
