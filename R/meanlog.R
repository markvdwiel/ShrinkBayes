meanlog <-
function(fit) {marg <- fit$marginals.hyper[[1]]; inla.expectation(function(x) log(x)^1, marg)}

