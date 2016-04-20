myinla.smarginal <-
function (marginal, log = FALSE, extrapolate = 0, keep.type = FALSE,len=300) 
{
    is.mat = is.matrix(marginal)
    m = INLA:::inla.marginal.fix(marginal)
    r = range(m$x)
    r = r[2] - r[1]
    ans = spline(m$x, log(m$y), xmin = min(m$x) - extrapolate * 
        r, xmax = max(m$x) + extrapolate * r, n = len)
    if (!log) {
        ans$y = exp(ans$y)
    }
    if (is.mat && keep.type) {
        return(cbind(ans$x, ans$y))
    }
    else {
        return(ans)
    }
}
