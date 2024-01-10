#' A Bootstrap Proportion Test for Brand Lift Testing (Liu et al., 2023)
#' @description This function generates binomial random samples for the control group
#' (with sample size \eqn{n_1} and success probability \eqn{p_1}) and the treatment
#' group (with sample size \eqn{n_2} and success probability \eqn{p_2}).
#' 
#' @import stats
#' @usage
#' gen.simu.data(n1, n2, p1, p2, summary=TRUE)
#' @param n1 sample size of the control group
#' @param n2 sample size of the treatment group
#' @param p1 success probability of the control group
#' @param p2 success probability of the treatment group
#' @param summary boolean variable. if TRUE it returns 2x2 contingency table; if FALSE it returns raw binomial random samples.
#' By default, summary=TRUE.
#' 
#' @return
#' A list of simulated data for the control group and the treatment group if \emph{summary=FALSE} or
#' a 2x2 contingency table if \emph{summary=TRUE}
#' @details 
#' The a 2x2 contingency table is of the following form
#' 
#' | col1 | col 2 |
#' |-----------------------|-----------------------------------|
#' | control sample size   | control positive response count   |
#' | treatment sample size | treatment positive response count |
#' @md
#' 
#' @export
#' @examples
#' n1 <- 100; n2 <- 100; p1 <- 0.1; p2 <- 0.2
#' set.seed(1)
#' sim.data <- gen.simu.data(n1, n2, p1, p2)
#' sim.data


gen.simu.data <- function (n1, n2, p1, p2, summary = TRUE){
  # n1: sample size in the control group
  # n2: sample size in the treatment group
  # p1: positive response rate in the control group
  # p2: positive response rate in the treatment group
  # summary: if True return 2x2 contingency table; if False return raw data
  # return: simulated survey response data
  
  control.data <- rbinom(n1, 1, p1)
  treatment.data <- rbinom(n2, 1, p2)
  simu.data <- list(control.data=control.data, treatment.data= treatment.data)
  
  if(summary){
    control.success <- sum(control.data)
    treatment.success <- sum(treatment.data)
    simu.data <- matrix(c(n1, n2, control.success, treatment.success), 2, 2)
  }
  
  return(simu.data)
}