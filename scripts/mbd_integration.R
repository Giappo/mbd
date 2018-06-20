install.packages("primes")
install.packages("gmp")
library(primes)
library(gmp)


if (1){
  draw.triangle <- function(tm, brts, brts.bars){
    ti.minus <- brts - brts.bars
    ti.plus  <- brts + brts.bars;
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
    return(out)
  }
  merge.triangles <- function(triangles){
    triangles.matrix <- matrix(unlist(triangles), nrow = length(triangles), byrow = T)
    out <- apply(X = triangles.matrix, MARGIN = 2, FUN = prod)
    return(out)
  }
  which.ones.intersect <- function(tm, brts, brts.bars){ -----
    t.minus <- brts-brts.bars; t.plus <- brts+brts.bars; t.minus <- t.minus[-1]; t.plus <- t.plus[-1];
    intervals <- rbind(t.minus, t.plus); N <- ncol(intervals)
    # for (i in 2:N){
    #  j <- which(t.minus[i] <= intervals[2, 1:(i-1)] & t.minus[i] >= intervals[1, 1:(i-1)])
    # }
    countm <- rep(0, length(tm))
    for (i in 1:N)
    {
      countm <- rbind(countm, i * ((tm <= intervals[2, i]) & (tm >= intervals[1, i])))
    }
    intersected <- apply(countm, MARGIN = 2, sort, decreasing = T)

    return(intersected)
  }
  # scenarios <- function(N){ ----
  #   C       <- t(partitions::restrictedparts(N, N))
  #   test    <- do.call(rbind, lapply(1:nrow(C),function(i) iterpc:::getall(iterpc::iterpc(table(C[i,]), order = T))));
  #   test2   <- split(x <- t(test), rep(1:ncol(x), each = nrow(x)))
  #   test3   <- lapply(test2, FUN = function(x){return(x[x != 0])})
  #   test4   <- unique(test3);test4
  #   ll      <- unlist(lapply(test4,length))
  #   test4.b <- lapply(test4, cumsum)
  #
  #   test5   <- rep(list(list()), length(test4))
  #   for (s in 1:length(test4))
  #   {
  #     test5[[s]] <- vector("list",length = ll[s])
  #     j<-1
  #     for (i in 1:ll[s])
  #     {
  #       test5[[s]][[i]] <- j:test4.b[[s]][i]
  #       j <- test4.b[[s]][i] + 1
  #       # test5[[s]][[i]] = (test4.b[[s]][i]):j
  #       # j <-test4.b[[s]][i]+1
  #     }
  #   }
  #   return(test5)
  # }
}

tm <- seq(-10, 0, 0.05)
brts <- c(-8, -7, -4, -1)
brts.bars <- c(1, 0.5, 2, 0.1)
vec <- 1:(Lb <- length(brts))
test <- list()
for (mm in 1:Lb){
  test[[mm]] <- combn(x = vec, m = mm)
}

test2 <- scenarios(10); test2[[3]]

tempo <- vector("list", N <- 10 )
for (n in 1:N)
{
 tempo[[n]] <- system.time(test2 <- scenarios(n))
};tempo

system.time(test2 <- scenarios(2))[3]

