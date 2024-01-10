#' A Bootstrap Proportion Test for Brand Lift Testing (Liu et al., 2023)
#' @description This function implements several proportion tests that can be applied 
#' to Brand Lift Testing, including
#' 
#' \enumerate{
#' \item \eqn{\mathbf{clt}}: Absolute lift based Z-test and relative lift based Z-test. The limiting 
#' distribution of Z-statistics are derived from the central limit theorem.
#' \item \eqn{\mathbf{bootstrap}}: Absolute lift based bootstrap test (BS-A) and relative lift based 
#' bootstrap test (BS-R), see Liu et al., (2023).
#' \item \eqn{\mathbf{bootstrapmean}}: Absolute lift based bootstrap mean test and relative lift 
#' based bootstrap mean test. (Efron and Tibshirani 1994).
#' \item \eqn{\mathbf{permutation}}: Absolute lift based permutation test and relative lift based 
#' permutation test. (Efron and Tibshirani 1994).
#' }
#' Learn more about the proportion tests in the section Details.
#' 
#' @import stats
#' @usage
#' proportion.test(data, method, B)
#' @param data A 2x2 matrix with first column being (control sample size, treatment sample size)
#' and the 2nd column being (control positive response count, treatment positive response count).
#' @param method The method should be one of ("clt", "bootstrap", 
#' "bootstrapmean", "permutation")
#' @param B Number of replications for bootstrap test or permutation test.
#' Only required for methods "bootstrap", "bootstrapmean", "permutation".
#' 
#' @return
#' A list of absolute lift, relative lift, standardized absolute lift and their corresponding
#' p-values. Standardized absolute lift equals absolute lift divided by its standard deviation.
#' Only absolute lift and relative lift are available for method clt.
#' @details 
#' \eqn{\mathbf{clt}}: the classic Z-test based on normal approximation. The absolute 
#' lift based Z-test is defined as
#' \deqn{
#' Z = \frac{\hat p_1 - \hat p_0}{\sqrt{{s_0^2}/{n_0} + {s_1^2}/{n_1}}},
#' }
#' and he relative lift based Z-test is defined as
#' \deqn{
#' Z_r = \frac{\hat p_1 / \hat p_0 - 1}{\sqrt{s_1^2/(n_1\hat p_0^2) + 
#' \hat p_1^2s_0^2/(n_0 \hat p_0^4)}},
#' }
#' where \eqn{s_0^2 = \hat p_0(1-\hat p_0)} and \eqn{s_1^2 = \hat p_1(1-\hat p_1)}.
#' 
#' \eqn{\mathbf{bootstrap}}: the bootstrap proportion tests proposed in Liu et al., (2023),
#' see Algorithm 1 in their paper. There are two bootstrap tests: the absolute lift 
#' based bootstrap test BS-A and the relative lift based bootstrap test BS-R. Note that
#' this type of bootstrap test is testing whether the distribution of the control group is
#' the same as the distribution of the treatment group. In the binomial distribution case,
#' it is equivalent to test whether the mean of the control group is
#' the same as the mean of the treatment group.
#' 
#' \eqn{\mathbf{bootstrapmean}} the bootstrap test to test whether the mean of the 
#' control group is the same as the mean of the treatment group. See Algorithm 16.2
#' of Efron and Tibshirani (1994).
#' 
#' \eqn{\mathbf{permutation}} the permutation test to test whether the distribution of the 
#' control group is the same as the distribution of the treatment group. 
#' See Algorithm 15.1 of Efron and Tibshirani (1994).
#' 
#' @export
#' @references
#' Wanjun Liu, Xiufan Yu, Jialiang Mao, Xiaoxu Wu, and Justin Dyer. 2023.
#' Quantifying the Effectiveness of Advertising: A Bootstrap Proportion Test
#' for Brand Lift Testing. \emph{In Proceedings of the 32nd ACM International Conference 
#' on Information and Knowledge Management (CIKM ’23)}
#' 
#' Efron, Bradley, and Robert J. Tibshirani. 
#' \emph{An introduction to the bootstrap}. CRC press, 1994.
#' 
#' @examples
#' n1 <- 100; n2 <- 100; p1 <- 0.1; p2 <- 0.2
#' set.seed(1)
#' sim.data <- gen.simu.data(n1, n2, p1, p2, summary = TRUE)
#' result <- proportion.test(sim.data, method = "bootstrap", B = 1000)
#' relative.lift <- result$lift$relative
#' relative.lift.pval <- result$pvalue$relative



## Perform two sample proportion test with different methods
proportion.test <- function(data, method, B){
  # data:2x2 contingency table
  # |-----------------------|-----------------------------------|
  # | control sample size   | control positive response count   |
  # |-----------------------|-----------------------------------|
  # | treatment sample size | treatment positive response count |
  # |-----------------------|-----------------------------------|
  # return: absolute lift and relative lift
  # method = ("clt", "bootstrap", "bootstrapmean", "permutation")
  #   -clt: standard proportion test with normal approximation
  #   -bootstrap: bootstrap test to test whether two distributions are the same
  #   -bootstrapmean: bootstrap test to test whether the means of two populations are the same
  #   -permutation: permutation test to test whether two distributions are the same
  # B: number of replications for bootstrap test and permutation test
  # return: absolute lift and relative lift together with p-values
  
  lift <- lift.calculator(data)
  absolute_lift <- lift$absolute_lift
  relative_lift <- lift$relative_lift
  std_absolute_lift <- lift$std_absolute_lift
  
  n1 <- data[1,1]
  n2 <- data[2,1]
  m1 <- data[1,2]
  m2 <- data[2,2]
  
  if(method == "clt"){
    
    p1 <- m1/n1
    p2 <- m2/n2
    v1 <- p1*(1-p1)
    v2 <- p2*(1-p2)
    
    if(m1==0 & m2==0){
      t_absolute <- 0
      t_relative <- 0
    } else if (m1==0){
      t_absolute <- (p2-p1)/sqrt(v1/n1 + v2/n2)
      t_relative <- sign(p2 - p1)*Inf
    } else {
      t_absolute <- (p2-p1)/sqrt(v1/n1 + v2/n2)
      vr <- v2/(n2*p1^2) + p2^2*v1/(n1*p1^4)
      t_relative <- relative_lift/sqrt(vr)
    }
    pvalue_absolute <- 2*(1-pnorm(abs(t_absolute)))
    pvalue_relative <- 2*(1-pnorm(abs(t_relative)))
    
    lift <- list(absolute=absolute_lift, 
                 relative=relative_lift)
    pvalue <- list(absolute=pvalue_absolute,
                   relative=pvalue_relative)
    
  }
  
  if(method == 'bootstrap'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    std_absolute_lift_B <- rep(NA, B)
    
    for (b in 1:B){
      n1b <- n1
      n2b <- n2
      m1b <- rbinom(1, n1b, (m1+m2)/(n1+n2))
      m2b <- rbinom(1, n2b, (m1+m2)/(n1+n2))
      datab <- matrix(c(n1b, n2b, m1b, m2b), 2, 2)
      liftb <- lift.calculator(datab)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      std_absolute_lift_B[b] <- liftb$std_absolute_lift
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B
    pvalue_std_absolute <- sum(abs(std_absolute_lift_B) > abs(std_absolute_lift))/B  
    
    lift <- list(absolute=absolute_lift, 
                 relative=relative_lift,
                 std_absolute=std_absolute_lift)
    pvalue <- list(absolute=pvalue_absolute,
                   relative=pvalue_relative, 
                   std_absolute=pvalue_std_absolute)
  }
  
  
  if(method == 'bootstrapmean'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    std_absolute_lift_B <- rep(NA, B)
    
    x1 <- c(rep(1, m1), rep(0, n1-m1))
    x2 <- c(rep(1, m2), rep(0, n2-m2))
    p1 <- mean(x1)
    p2 <- mean(x2)
    p12 <- mean(c(x1, x2))
    x1n <- x1 - p1 + p12
    x2n <- x2 - p2 + p12
    
    for (b in 1:B){
      
      x1b <- sample(x1n, n1, replace=T)
      x2b <- sample(x2n, n2, replace=T)
      liftb <- lift.calculator(x1b, x2b)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      std_absolute_lift_B[b] <- liftb$std_absolute_lift
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B  
    pvalue_std_absolute <- sum(abs(std_absolute_lift_B) > abs(std_absolute_lift))/B  
    
    lift <- list(absolute=absolute_lift, 
                 relative=relative_lift,
                 std_absolute=std_absolute_lift)
    pvalue <- list(absolute=pvalue_absolute,
                   relative=pvalue_relative, 
                   std_absolute=pvalue_std_absolute)
  }
  
  
  if(method == 'permutation'){
    
    absolute_lift_B <- rep(NA, B)
    relative_lift_B <- rep(NA, B)
    std_absolute_lift_B <- rep(NA, B)
    x_pooled <- c(rep(1, m1+m2), rep(0, n1+n2-m1-m2))
    
    for (b in 1:B){
      xb <- sample(x_pooled)
      n1b <- n1
      n2b <- n2
      m1b <- sum(xb[1:n1])
      m2b <- sum(xb[(n1+1):(n1+n2)])
      datab <- matrix(c(n1b, n2b, m1b, m2b), 2, 2)
      liftb <- lift.calculator(datab)
      absolute_lift_B[b] <- liftb$absolute_lift
      relative_lift_B[b] <- liftb$relative_lift
      std_absolute_lift_B[b] <- liftb$std_absolute_lift
    }
    pvalue_absolute <- sum(abs(absolute_lift_B) > abs(absolute_lift))/B
    pvalue_relative <- sum(abs(relative_lift_B) > abs(relative_lift))/B  
    pvalue_std_absolute <- sum(abs(std_absolute_lift_B) > abs(std_absolute_lift))/B
    
    lift <- list(absolute=absolute_lift, 
                 relative=relative_lift,
                 std_absolute=std_absolute_lift)
    pvalue <- list(absolute=pvalue_absolute,
                   relative=pvalue_relative, 
                   std_absolute=pvalue_std_absolute)
  }
  
  return(list(lift=lift, pvalue=pvalue))
}


## Calculate the absolute lift and relative lift
lift.calculator <- function(x, y = NULL){
  # data:2x2 contingency table
  # |-----------------------|-----------------------------------|
  # | control sample size   | control positive response count   |
  # |-----------------------|-----------------------------------|
  # | treatment sample size | treatment positive response count |
  # |-----------------------|-----------------------------------|
  # return: absolute lift and relative lift
  
  if (is.matrix(x)){
    if(ncol(x) != 2L || nrow(x) != 2L)
      stop("'x' must have 2 columns and 2 rows")
    n1 <- x[1,1]
    n2 <- x[2,1]
    m1 <- x[1,2]
    m2 <- x[2,2]
  }
  else if (is.vector(x) && is.vector(y)){
    n1 <- length(x)
    n2 <- length(y)
    m1 <- sum(x)
    m2 <- sum(y)
  } 
  else {
    stop("Input data needs to be a 2x2 matirx (x) or two vectors (x and y)")
  }
  
  prr1 <- m1/n1
  prr2 <- m2/n2
  absolute_lift <- prr2 - prr1
  
  if(m1==0 & m2==0){
    relative_lift <- 0
  } else if (m1==0) {
    relative_lift <- Inf
  } else {
    relative_lift <- (prr2 - prr1)/prr1
  }
  
  v1 <- prr1 * (1 - prr1)
  v2 <- prr2 * (1 - prr2)
  
  if(prr1==prr2 & v1==0 & v2==0){
    std_absolute_lift <- 0.0
  } else if (v1==0 & v2==0){
    std_absolute_lift <- sign(prr2 - prr1)*Inf
  } else {
    std_absolute_lift <- (prr2-prr1)/sqrt(v1/n1 + v2/n2)
  }
  
  return(list(absolute_lift=absolute_lift, 
              relative_lift=relative_lift,
              std_absolute_lift=std_absolute_lift))
  
}