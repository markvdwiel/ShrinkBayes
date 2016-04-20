linear.spline <-
function (marginal) 
{
    m = INLA:::inla.marginal.fix(marginal)
    r = range(m$x)
    return(list(range = r, fun = approxfun(m$x, log(m$y))))
}
