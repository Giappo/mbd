% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbd_condprob.R
\name{condprob_dq_nu}
\alias{condprob_dq_nu}
\title{Calculates the nu component of dq}
\usage{
condprob_dq_nu(qq, nu_matrix)
}
\arguments{
\item{qq}{the matrix q_{m1, m2}}

\item{nu_matrix}{matrix for the nu component. Component {m,n} is equal to
choose(n, m - n) x q ^ (m - n) x (1 - q) ^ (2 * n - m)}
}
\description{
Calculates the nu component of dq
}
