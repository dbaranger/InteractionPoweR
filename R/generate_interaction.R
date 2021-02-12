
#' Generate interaction data set
#'
#' Simulate a single data set with an interaction (y ~ x1 + x2 + x1*x2). All values other than 'N' are population-level effects - the values within any single simulated data set will vary around the defined values.
#'
#' @param N Sample size. Must be a positive integer. Has no default value.
#' @param r.x1.y Pearson's correlation between x1 and y. Must be between -1 and 1.. Has no default value.
#' @param r.x2.y Pearson's correlation between x2 and y. Must be between -1 and 1.. Assumed to be the 'moderator' in some functions. Has no default value.
#' @param r.x1x2.y Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1.. Has no default value.
#' @param r.x1.x2 Pearson's correlation between x1 and x2. Must be between -1 and 1.. Has no default value.
#' @param sd.x1 Standard deviation of x1. Defaults to 1.
#' @param sd.x2 Standard deviation of x2. Defaults to 1.
#' @param sd.y Standard deviation of y. Defaults to 1.
#' @param mean.x1 Mean of x1. Defaults to 0.
#' @param mean.x2 Mean of x2. Defaults to 0.
#' @param mean.y Mean of y. Defaults to 0.
#'
#' @return A data frame containing variables 'x1', 'x2', 'y', and 'x1x2'. 'x1x2' is x1*x2. The correlations between these variables are drawn from the defined population-level values.
#' @export
#'
#' @examples
#' \dontrun{dataset <- generate_interaction(N = 250,r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' cor(dataset)
#' }
generate_interaction <- function(N,r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,sd.x1=1,sd.x2=1,sd.y=1,mean.x1=0,mean.x2=0,mean.y=0) {

  # compute sd of x1x2 (interaction term) using x1 and x2

  cor.mat = matrix(data = c(1,r.x1.x2,
                            r.x1.x2,1),nrow = 2,byrow = T)
  sd = c(sd.x1,sd.x2  )

  cov_x1x2<- diag(sd) %*% cor.mat %*% diag(sd) # cor 2 cov

  sd.x1x2     =  base::sqrt((cov_x1x2[1,1]*cov_x1x2[2,2]) + ((cov_x1x2[1,2])^2  ))
  r.x1.x1x2   =  0
  r.x2.x1x2   =  0



  ## Path tracing - compute values of coefficients in regression model y~ x1 + x2 + x1*x2

  #construct covariance matrix

  #x1        x2          x1x2        y
  cormat <- base::matrix(data = c(1,         r.x1.x2,    r.x1.x1x2,  r.x1.y,   #x1
                                  r.x1.x2,   1,          r.x2.x1x2,  r.x2.y,   #x2
                                  r.x1.x1x2, r.x2.x1x2,  1,          r.x1x2.y, #x1x2
                                  r.x1.y,    r.x2.y,     r.x1x2.y,   1),       #y
                         ncol = 4, byrow = TRUE)

  if(base::max(base::abs(cormat[lower.tri(cormat)]))>1){
    print("All correlations must be within [-1,1]")
    stop()
  }
  sd = c(sd.x1,sd.x2,sd.x1x2,sd.y   )
  #covmat<- cor2cov_1(cor.mat =cormat,)
  covmat<- diag(sd) %*% cormat %*% diag(sd) # cor 2 cov


  # actual path tracing
  b<-covmat[4,-4] # the y-row
  A<-covmat[-4,-4] # everything but the y-row
  betas<-base::solve(A,b)


  # coefficients
  b1<-betas[1]
  b2<-betas[2]
  b3<-betas[3]

  # simulate x1 and x2 as continuous normal variables

  data = MASS::mvrnorm(n=N, mu=c(mean.x1, mean.x2), Sigma=covmat[c(1:2),c(1:2)], empirical=TRUE)

  x1    = data[, 1]
  x2    = data[, 2]

  var.y<-covmat[4,4]
  var.x1<-covmat[1,1]
  var.x2<-covmat[2,2]

  # simulate y
  # sqrt because rnorm() takes SD, not variance

  ysd <- base::sqrt(var.y -
                      ((b1^2)*var.x1) -
                      ((b2^2)*var.x2) -
                      ((b3^2)*( (var.x1*var.x2) + (covmat[1,2])^2) ) -
                      (2*(b1*b2)* covmat[1,2] )
  )

  ymean <- (mean.y + b1*x1 + b2*x2 + b3*x1*x2)

  y <- stats::rnorm(N, ymean , ysd)

  x1x2<-x1*x2

  # output dataframe
  dat<-as.data.frame(cbind(x1,x2,y,x1x2))
  colnames(dat)<-c("x1","x2","y","x1x2")
  return(dat)
}
