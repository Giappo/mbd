# rm(list = ls())
if (1) {#functions
  draw.triangle <- function(tm, brt, brt.bar){
    ti.minus <- brt - brt.bar
    ti.plus  <- brt + brt.bar;
    half_interval <- abs(ti.plus - ti.minus)/2  #brts.bars?
    ti <- (ti.minus + ti.plus)/2 #brts?
    m  <- 1/(half_interval)^2
    s <- sign(tm - ti)
    s1 <- sign(abs(tm - ti) - half_interval) #<0 it is in the interval; >0 is outside
    
    out = ((1 - s1)/2) *
      (
        ((1 - s)/2) * ( m * (tm - ti.minus)          ) +
          ((1 + s)/2) * (-m * (tm - ti      ) + sqrt(m))
      )
    
    if (all(out == 0))
    {
      x <- abs(tm - ti)
      out[which(x == min(x))] <- 1
    }
    return(out)
  }
  draw.triangles<- function(tm, brts, brts.bars){
    brts2 <- brts[-1]
    brts.bars2 <- brts.bars[-1]
    triangles <- vector("list", N <- length(brts2))
    for (i in 1:N)
    {
      triangles[[i]] <- draw.triangle(tm = tm, brt = brts2[i], brt.bar = brts.bars2[i])
    }
    return(triangles)
  }
  sample.brts <- function(tm, brts, brts.bars) {
    triangles <- draw.triangles(tm = tm, brts = brts, brts.bars = brts.bars)
    sampled.brts <- c(max(abs(brts)), rep(NA, length(brts) - 1))
    for (i in 2:length(brts))
    {
      sampled.brts[i] <- sample(x = tm, size = 1, prob = triangles[[i - 1]]) 
    }
    return(sampled.brts)
  }
}

N <- 5
age <- -10
tm <- seq(age, 0, 0.1)
brts <- c(age, sort(runif(n = N, min = -9.9, max = 0.1), decreasing = F))
brts.bars <- c(0, abs(rnorm(n = N, mean = 0, sd = 0.5)))

max_step <- 20; res <- rep(0, max_step); repetitions <- 1E3
for (i in 1:max_step)
{
  Nsteps <- (10 * i)
  tm2 <- seq(age, 0, abs(age - 0)/Nsteps)
  for (r in 1:repetitions)
  {
    test <- sample.brts(tm = tm2, brts = brts, brts.bars = brts.bars)
    births <- MBD:::brts2time_intervals_and_births(test)$births
    res[i] <- res[i] + sum(births > 1)
  }
  res[i] <- res[i]/(repetitions * length(births))
}
plot(res~(Nsteps <- (10 * 1:max_step)))

