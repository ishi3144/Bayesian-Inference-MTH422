H = function(n,x) 
{
    value = 0
    if (n == 0) value = 1
    if (n == 1) value = x
    if (n >= 2) 
    {
        HH = rep(0,n+1)
        HH[1] = 1
        HH[2] = x
        for (nn in c(2:n)) 
        {
            HH[nn + 1] = x * HH[nn] - (nn - 1) * HH[nn-1]
        }

        value = HH[n + 1]
    }
    value
}
