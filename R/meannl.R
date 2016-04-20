meannl <-
function(fit) {marg <- fit$marginals.hyper$Size; inla.expectation(function(x) x, marg)}

