FindStepLengthBinomial <- function(X, y, lam, pmax=10000, nmin=1000) {
     ## add more sophistication later.
     Lmax <- norm(XX, "2")^2
     stplngth <- 8/Lmax
     return(stplngth)
}
