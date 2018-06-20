# brts;
pars <- c(0.3, 0, 1.5, 0.1);
cond <- 1; soc <- 2; tips_interval <- c(0, Inf);
alpha0 <- 10
methode <- 'expo'; abstol <- 1e-16; reltol <- 1e-10; minimum_multiple_births <- 0

Pc_alpha_calibration <- function(brts, pars, cond, soc, tips_interval, alpha0, methode,
                                 abstol, reltol, minimum_multiple_births){
alpha <- alpha0
deltaAlpha <- 1; count <- 0; same_result_count <- 0; Pc.notanumber <- 1;
while (Pc.notanumber)
{

  Pc1 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, cond = cond,
                                                 soc = soc, tips_interval = tips_interval,
                                                 alpha = alpha, methode = methode,
                                                 abstol = abstol, reltol = reltol,
                                                 minimum_multiple_births = minimum_multiple_births)$Pc

  Pc.notanumber <- is.nan(Pc1)
  alpha <- alpha - Pc.notanumber
}
while (deltaAlpha != 0 && count < 100 && same_result_count < 5)
{
  Pc2 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, cond = cond,
                                                 soc = soc, tips_interval = tips_interval,
                                                 alpha = alpha + deltaAlpha, methode = methode,
                                                 abstol = abstol, reltol = reltol,
                                                 minimum_multiple_births = minimum_multiple_births)$Pc
  if (is.nan(Pc2))
  {
    deltaAlpha <- deltaAlpha - 1
  }else if (Pc2 < Pc1)
  {
    deltaAlpha <- deltaAlpha + 1
    same_result_count <- same_result_count + 1
  }else
  {
    same_result_count <- 0
    deltaPc <- abs(Pc2 - Pc1)/Pc1;
    # deltaAlpha = floor(  10*(-1 + 2/( 1 + exp(-(1/2*deltaPc)) ))  )
    deltaAlpha <- floor(10 * deltaPc)
    alpha <- alpha + deltaAlpha
    Pc1 <- Pc2
  }

  count = count + 1
  # print(alpha)
}
if (max_k * alpha >= 2000)
{#check to see whether alpha is too big to be handled without memory issues
  alpha <- floor(1500/max_k);
  Pc1 <- MBD:::calculate_conditional_probability(brts = brts, pars = pars, tips_interval = tips_interval, cond = cond, soc = soc, alpha = alpha, methode = methode,
                                                 abstol = abstol, reltol = reltol, minimum_multiple_births = minimum_multiple_births)$Pc
}
Pc <- Pc1
if (count >= 100){alpha <- 10}
if (Pc <= 0 | Pc == Inf | Pc == -Inf){Pc <- 1; print("there's a problem with Pc")}
return(list(Pc = Pc, alpha = alpha))
}
