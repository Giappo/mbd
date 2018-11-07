#' This function does nothing. It is intended to inherit is parameters'
#' documentation.
#' @param a_abstol something
#' @param a_reltol something
#' @param abstol absolute error tolerance for
#' the numerical integration using deSolve.
#' @param account_name something
#' @param after something
#' @param age the age of the tree.
#' @param alpha something
#' @param alpha0 something
#' @param b the number of simultaneous births on a given branching time.
#' @param brts A set of branching times of a phylogeny.
#' @param changeloglikifnoconv something
#' @param checked_functions something
#' @param colormap something
#' @param cond set 1 if you want to condition on stem or crown age
#'   and non-extinction of the phylogeny_ Set 0 otherwise.
#' @param function_name something
#' @param idparsfix something
#' @param idparsopt the ids of the parameters that must be optimized.
#'   The ids are defined as follows:
#'   \itemize{
#'     \item pars[1] is lambda, the sympatric speciation rate;
#'     \item pars[2] is mu, the extinction rate;
#'     \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'     \item pars[4] is q, the single-lineage speciation probability_
#'   }
#' @param initparsopt something
#' @param input_trees_path Path to \code{.trees} file.
#' @param input_xml something
#' @param interval_max something
#' @param interval_min something
#' @param iterations something
#' @param k the number of visible species in the phylogeny at a given time.
#' @param lambda the sympatric speciation rate.
#' @param logs something
#' @param lx it is the number of ODEs considered for the computation.
#' @param lx0 something
#' @param matrix something
#' @param matrix_builder function used to build the transition matrix.
#' Default option is \code{hyper_a_hanno}
#' @param max_iter Sets the maximum number of iterations in the optimization
#' @param max_k something
#' @param max_number_of_species something
#' @param max_repetitions something
#' @param max_sims something
#' @param maxiter something
#' @param mbd_lambda something
#' @param methode
#'   specifies how the integration must be performed:
#'   \itemize{
#'     \item \code{sexpm}: use \code{sexpm}
#'     \item \code{expo}: use \code{expoRkit}
#'     \item \code{lsoda}: use \code{lsoda} and \code{deSolve::ode}
#'   }
#' @param minimum_multiple_births minimum amount of multiple births
#' that have to be present in the simulated phylogeny.
#' @param missnumspec The number of species that are in the clade,
#'   but missing in the phylogeny.
#' @param mu the extinction rate.
#' @param mutation_rate something
#' @param n_0 the number of lineages at time equals zero.
#' @param n_species number of species.
#' @param init_n_lineages the number of lineages at time equals zero.
#' @param n_steps something
#' @param n_subs something
#' @param nu the multiple allopatric speciation trigger rate.
#' @param optimmethod optimization routine: choose between subplex and simplex
#' @param params transition matrix for the rhs of the ode system.
#' @param pars vector of parameters:
#' \itemize{
#'   \item pars[1] is lambda, the sympatric speciation rate;
#'   \item pars[2] is mu, the extinction rate;
#'   \item pars[3] is nu, the multiple allopatric speciation trigger rate;
#'   \item pars[4] is q, the single-lineage speciation probability;
#' }
#' @param pars_transform something
#' @param parsfix the values of the parameters that should not be optimized.
#' @param precision something
#' @param printit something
#' @param print_errors something
#' @param q_vector the \code{Q} vector from
#' 'Etienne et al. 2012 - Proc. R. Soc. B'.
#' @param q the single-lineage speciation probability at a triggered event.
#' @param quantiles_choice something
#' @param recursive something
#' @param reltol relative error tolerance for
#' the numerical integration using deSolve.
#' @param res something
#' @param results something
#' @param s something
#' @param safety_threshold adds a threshold for the evaluation of q. This is due
#' because you never want \code{q} to actually be equal to zero or one.
#' @param sample_interval something
#' @param sequence_length something
#' @param sim_pars vector of parameters:
#' \itemize{
#'   \item id == 1 corresponds to lambda (speciation rate)
#'   \item id == 2 corresponds to mu (extinction rate)
#'   \item id == 3 corresponds to nu (multiple speciation trigger rate)
#'   \item id == 4 corresponds to q (single-lineage speciation probability)
#' }
#' @param sim_phylo something
#' @param subsamp something
#' @param t something
#' @param t1 something
#' @param t2 something
#' @param time_interval something
#' @param tips_interval sets tips boundaries constrain on simulated dataset.
#'   You can also define the tips_interval as you can usually
#'   do with a standard usage of mbd_sim.
#' @param tol something
#' @param tr a phylogeny of class \code{phylo}
#' @param transition_matrix something
#' @param trparsfix something
#' @param trparsopt something
#' @param values something
#' @param verbose choose if you want to print the output or not
#' @param x something
#' @param x_name something
#' @param x_splits something
#' @param y something
#' @param y_name something
#' @param y_splits something
#' @param z something
#' @param z_name something
#' @param start_pars starting parameter values for the ml process.
#' @param true_pars true parameter values when running the ml process.
#' @param optim_ids ids of the parameters you want to optimize.
#' @author Documentation by Giovanni Laudanno,
#'   use of this function by Richel J.C. Bilderbeek
#' @note This is an internal function, so it should be marked with
#'   \code{@noRd}. This is not done, as this will disallow all
#'   functions to find the documentation parameters
default_params_doc <- function(
  a_abstol,
  a_reltol,
  abstol,
  account_name,
  after,
  age,
  alpha,
  alpha0,
  b,
  brts,
  changeloglikifnoconv,
  checked_functions,
  colormap,
  cond,
  function_name,
  idparsfix,
  idparsopt,
  initparsopt,
  input_trees_path,
  input_xml,
  interval_max,
  interval_min,
  iterations,
  k,
  lambda,
  logs,
  lx,
  lx0,
  matrix,
  matrix_builder,
  max_iter,
  max_k,
  max_number_of_species,
  max_repetitions,
  max_sims,
  maxiter,
  mbd_lambda,
  methode,
  minimum_multiple_births,
  missnumspec,
  mu,
  mutation_rate,
  n_0,
  n_species,
  init_n_lineages,
  n_steps,
  n_subs,
  nu,
  optimmethod,
  params,
  pars,
  pars_transform,
  parsfix,
  precision,
  print_errors,
  printit,
  q_vector,
  q,
  quantiles_choice,
  recursive,
  reltol,
  res,
  results,
  s,
  safety_threshold,
  sample_interval,
  sequence_length,
  sim_pars,
  sim_phylo,
  subsamp,
  t,
  t1,
  t2,
  tips_interval,
  time_interval,
  tol,
  tr,
  transition_matrix,
  trparsfix,
  trparsopt,
  values,
  verbose,
  x,
  x_name,
  x_splits,
  y,
  y_name,
  y_splits,
  z,
  z_name,
  start_pars,
  true_pars,
  optim_ids
) {
  # Nothing
}
