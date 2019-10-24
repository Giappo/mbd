context("cond_probs")

test_that("nu differential is the same with every approach", {

q <- 0.1
lx <- 60
pvec <- matrix(0, nrow = lx, ncol = lx)
pvec[2:8, 2:3] <- 1; pvec <- pvec / sum(pvec)
pvec <- matrix(pvec, nrow = lx ^ 2, ncol = 1)
lx2 <- length(pvec)
lx <- sqrt(lx2)

mm <- 2:(lx + 1)
mm_plus_one <- mm + 1
mm_minus_one <- mm - 1

pp <- matrix(pvec, nrow = lx, ncol = lx)
# pp2 <- matrix(0, nrow = (lx + 2), ncol = (lx + 2))
# pp2[mm, mm] <- pp

nu_matrix <- matrix(0, nrow = lx, ncol = lx)
# alt: for (m1 in 0:(lx - 1)) {
# alt: for (n1 in 0:m1) {
for (n1 in 0:(lx - 1)) {
  for (m1 in n1:(lx - 1)) {
    aux <- lchoose(n1, max(0, m1 - n1))
    aux <- aux + (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
    nu_matrix[m1 + 1, n1 + 1] <- exp(aux)
  }
}
rownames(nu_matrix) <- paste0("m1=", 0:(lx - 1))
colnames(nu_matrix) <- paste0("n1=", 0:(lx - 1))

nu_matrix2 <- matrix(0, nrow = lx, ncol = lx)
for (m1 in 0:(lx - 1)) {
  for (a1 in 0:floor((m1 + 1) / 2)) {
    aux <- lchoose(m1 - a1, a1)
    aux <- aux + (a1) * log(q) + (m1 - 2 * a1) * log(1 - q)
    nu_matrix2[m1 + 1, m1 - a1 + 1] <- exp(aux)
  }
}
rownames(nu_matrix2) <- paste0("m1=", 0:(lx - 1))
colnames(nu_matrix2) <- paste0("n1=", 0:(lx - 1))

nu_matrix3 <- matrix(0, nrow = lx, ncol = lx)
for (m1 in 0:(lx - 1)) {
  for (n1 in 0:(m1 - 1)) {
    aux <- lchoose(n1, max(0, m1 - n1))
    aux <- aux + (m1 - n1) * log(q) + (2 * n1 - m1) * log(1 - q)
    nu_matrix3[m1 + 1, n1 + 1] <- exp(aux)
  }
}
rownames(nu_matrix3) <- paste0("m1=", 0:(lx - 1))
colnames(nu_matrix3) <- paste0("n1=", 0:(lx - 1))

nu_matrix4 <- matrix(0, nrow = lx, ncol = lx)
for (m1 in 0:(lx - 1)) {
  for (a1 in 1:floor((m1 + 1) / 2)) {
    aux <- lchoose(m1 - a1, a1)
    aux <- aux + (a1) * log(q) + (m1 - 2 * a1) * log(1 - q)
    nu_matrix4[m1 + 1, m1 - a1 + 1] <- exp(aux)
  }
}
rownames(nu_matrix4) <- paste0("m1=", 0:(lx - 1))
colnames(nu_matrix4) <- paste0("n1=", 0:(lx - 1))

# approaches with the zero terms produce same nu_matrix
testthat::expect_equal(
  nu_matrix,
  nu_matrix2
)
# approaches without zero terms produce same nu_matrix
testthat::expect_equal(
  nu_matrix3,
  nu_matrix4
)

dp1 <- nu_matrix %*% pp %*% t(nu_matrix) - pp
dp2 <- nu_matrix2 %*% pp %*% t(nu_matrix2) - pp
sum(dp1)
sum(dp2)

diag_matrix <- diag((1 - q) ^ (0:(lx - 1)))

loss_term <-
  diag_matrix %*% pp %*% t(nu_matrix) +
  nu_matrix %*% pp %*% t(diag_matrix) -
  diag_matrix %*% pp %*% t(diag_matrix) -
  pp
dp3 <- nu_matrix3 %*% pp %*% t(nu_matrix3) + loss_term
sum(dp3)

loss_term <-
  diag_matrix %*% pp %*% t(nu_matrix3) +
  nu_matrix3 %*% pp %*% t(diag_matrix) +
  diag_matrix %*% pp %*% t(diag_matrix) -
  pp
dp4 <- nu_matrix3 %*% pp %*% t(nu_matrix3) + loss_term
sum(dp4)

# approaches with the zero terms produce same differential
testthat::expect_equal(
  dp1,
  dp2
)
# all approaches produce same differential, as the zero term is compensated by
# the negative term
testthat::expect_equal(
  dp1,
  dp3
)
testthat::expect_equal(
  dp1,
  dp4
)
})
