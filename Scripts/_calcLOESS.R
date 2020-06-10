calcLOESS <- function(Var1, Var2){
  # Var1 = x-column for LOESS (e.g. the metric)
  # Var2 = y-column for LOESS (e.g. age)
  
  # first check all SSE starting from a smoothing kernel of 5% till 100%
  param <- optim(par=c(0.05),x = Var1, y = Var2, calcSD, method="Brent", lower = 0.05, upper = 1)
  
  # run again with optimized parameters on the control data
  loessModOptim <- loess.sd(Var1 ~ Var2, span=param$par)
  var_mu <- loessModOptim$y
  var_sd <- loessModOptim$sd
  
  out <- as.data.frame(cbind(var_mu, var_sd))
  colnames(out) <- c("Mu", "SD")
  
  return(out)  
}