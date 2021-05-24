
#' Generate interaction data set
#'
#' Simulate a single data set with an interaction (y ~ x1 + x2 + x1*x2). All values other than 'N' are population-level effects - the values within any single simulated data set will vary around the defined values.
#'
#' @param N Sample size. Must be a positive integer. Has no default value.
#' @param r.x1.y Pearson's correlation between x1 and y. Must be between -1 and 1. Has no default value.
#' @param r.x2.y Pearson's correlation between x2 and y. Must be between -1 and 1. Assumed to be the 'moderator' in some functions. Has no default value.
#' @param r.x1x2.y Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1. Has no default value.
#' @param r.x1.x2 Pearson's correlation between x1 and x2. Must be between -1 and 1. Has no default value.
#' @param rel.x1 Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.x2 Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.y Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param skew.x1 Skew of x1. Default is 0 (normally distributed).
#' @param skew.x2 Skew of x2. Default is 0 (normally distributed).
#' @param skew.y Skew of y. Default is 0 (normally distributed).
#' @param adjust.correlations If variables are skewed or binary, should correlations be adjusted so that output data has the specified correlation structure? Default is TRUE.
#' @param transform.x1 Transform x1? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma" if variable is skewed.
#' @param transform.x2 Transform x2? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma" if variable is skewed.
#' @param transform.y Transform y? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma"if variable is skewed.
#' @param r.x1.y.adjust Internal use only
#' @param r.x2.y.adjust Internal use only
#' @param r.x1x2.y.adjust Internal use only
#' @param r.x1.x2.adjust Internal use only
#' @return A data frame containing variables 'x1', 'x2', 'y', and 'x1x2'. 'x1x2' is x1*x2. The correlations between these variables are drawn from the defined population-level values.
#' @export
#'
#' @examples
#' \dontrun{dataset <- generate_interaction(N = 250,r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' cor(dataset)
#' }

generate_interaction <- function(N,
                                     r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,
                                     rel.x1=1,rel.x2=1,rel.y=1,
                                     skew.x1 = 0,
                                     skew.x2 = 0,
                                     skew.y = 0,
                                     transform.x1 = "default",
                                     transform.x2 = "default",
                                     transform.y = "default",
                                     adjust.correlations = T,
                                     r.x1.y.adjust=NULL,
                                     r.x2.y.adjust=NULL,
                                     r.x1x2.y.adjust=NULL,
                                     r.x1.x2.adjust=NULL
) {

  # order of operations
  #  1. Check for input errors
  #  2. Set defaults
  #  3. Adjust correlations for distribution transformations (if requested)
  #  4. Path tracing
  #  5. Simulate data as continuous normal variables
  #  6. Add in noise from reliability
  #  7. Transform variables (skew & number of levels)
  #  8. Output




  # check for input errors
  if(rel.x1 <= 0 | rel.x1 >1 |rel.x2 <= 0 | rel.x2 >1 | rel.y <= 0 | rel.y >1){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}



  if(max(abs(c( r.x1.y,r.x2.y,r.x1.x2,r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1]")}

  # defaults

  sd.x1 = 1
  sd.x2 = 1
  sd.y  = 1
  mean.x1 = 0
  mean.x2 = 0
  mean.y  = 0
  r.x1.x1x2 = 0
  r.x2.x1x2 = 0


  if(transform.x1 == "default"){
    if(skew.x1 != 0){transform.x1 = "gamma"}
    # if(skew.x1 < 0 & skew.x1 >= -1){transform.x1 ="sn"}
    #if(skew.x1 < -1){transform.x1 = "beta"}
  }
  if(transform.x2 == "default"){
    if(skew.x2 != 0){transform.x2 = "gamma"}
    # if(skew.x2 < 0 & skew.x2 >= -1){transform.x2 ="sn"}
   # if(skew.x2 < -1){transform.x2 = "beta"}
  }
  if(transform.y == "default"){
    if(skew.y != 0){transform.y = "gamma"}
    #  if(skew.y < 0 & skew.y >= -1){transform.y = "sn"}
   # if(skew.y < -1){transform.y = "beta"}
  }




  if(skew.x1 == 0 & skew.x2 == 0 & skew.y == 0 &
     transform.x1 == "default" & transform.x2 == "default" & transform.y == "default"){
    adjust.correlations<-F}

  if(adjust.correlations == T){
    adjustments<-compute_adjustment(r.x1.y = r.x1.y,r.x2.y = r.x2.y,
                                    r.x1x2.y = r.x1x2.y,r.x1.x2 = r.x1.x2,
                                    skew.x1 = skew.x1,skew.x2 = skew.x2,skew.y = skew.y,
                                    #levels.x1 = levels.x1, levels.x2 = levels.x2, levels.y = levels.y,
                                    transform.x1 = transform.x1,transform.x2 = transform.x2,transform.y = transform.y
    )

    r.x1.y.adjust = adjustments[2]
    r.x2.y.adjust = adjustments[4]
    r.x1.x2.adjust = adjustments[1]
    r.x1x2.y.adjust = adjustments[6]

  }

  # add adjustment to correlations for subsequent skew attenuation

  if(!is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust)){
    r.x1.y = r.x1.y + r.x1.y.adjust
    r.x2.y = r.x2.y + r.x2.y.adjust
    r.x1.x2 = r.x1.x2 + r.x1.x2.adjust
    r.x1x2.y = r.x1x2.y + r.x1x2.y.adjust

  }


  if(max(abs(c( r.x1.y,r.x2.y,r.x1.x2,r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1], this could be due to skew being too high.")
  }


  # compute sd of x1x2 (interaction term) using x1 and x2

  cor.mat = matrix(data = c(1,r.x1.x2,
                            r.x1.x2,1),nrow = 2,byrow = T)
  sd = c(sd.x1,sd.x2  )

  cov_x1x2<- diag(sd) %*% cor.mat %*% diag(sd) # cor 2 cov

  sd.x1x2     =  base::sqrt((cov_x1x2[1,1]*cov_x1x2[2,2]) + ((cov_x1x2[1,2])^2  ))




  ## Path tracing - compute values of coefficients in regression model y~ x1 + x2 + x1*x2

  #construct covariance matrix

  #x1        x2          x1x2        y
  cormat <- base::matrix(data = c(1,         r.x1.x2,    r.x1.x1x2,  r.x1.y,   #x1
                                  r.x1.x2,   1,          r.x2.x1x2,  r.x2.y,   #x2
                                  r.x1.x1x2, r.x2.x1x2,  1,          r.x1x2.y, #x1x2
                                  r.x1.y,    r.x2.y,     r.x1x2.y,   1),       #y
                         ncol = 4, byrow = TRUE)


  sd = c(sd.x1,sd.x2,sd.x1x2,sd.y   )
  covmat<- diag(sd) %*% cormat %*% diag(sd) # cor 2 cov


  #  path tracing
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

  yvar <-              var.y -
    ((b1^2)*var.x1) -
    ((b2^2)*var.x2) -
    ((b3^2)*( (var.x1*var.x2) + (covmat[1,2])^2) ) -
    (2*(b1*b2)* covmat[1,2] )

  if(yvar < 0){
    #print()
    stop("Settings produce a negative y-variance")
  }

  ysd <- base::sqrt(yvar)
  ymean <- (mean.y + b1*x1 + b2*x2 + b3*x1*x2)

  y <- stats::rnorm(N, ymean , ysd)



  # Add measurement error from reliability

  std_x1_err=sqrt((1-rel.x1)/rel.x1)
  std_x2_err=sqrt((1-rel.x2)/rel.x2)
  std_y_err=sqrt((1-rel.y)/rel.y)

  error_x1=stats::rnorm(N, 0, std_x1_err)
  error_x2=stats::rnorm(N, 0, std_x2_err)
  error_y=stats::rnorm(N, 0, std_y_err)

  y <- base::scale(x = (y+error_y),center = T,scale = T)
  x1<- base::scale(x = (x1+error_x1),center = T,scale = T)
  x2<- base::scale(x = (x2+error_x2),center = T,scale = T)

  # add levels


  #x1x2<-x1*x2

  ### skew ###

  if(transform.x1 == "binary"){x1 = norm2binary(x = x1,skew = skew.x1)}
  if(transform.x2 == "binary"){x2 = norm2binary(x = x2,skew = skew.x2)}
  if(transform.y == "binary"){y = norm2binary(x = y,skew = skew.y)}

  if(transform.x1 == "gamma"){x1 = norm2gamma(x = x1,skew = skew.x1)}
  if(transform.x2 == "gamma"){x2 = norm2gamma(x = x2,skew = skew.x2)}
  if(transform.y == "gamma"){y = norm2gamma(x = y,skew = skew.y)}

  # if(transform.x1 == "beta"){x1 = norm2gamma2beta(x = x1,skew = skew.x1)}
  # if(transform.x2 == "beta"){x2 = norm2gamma2beta(x = x2,skew = skew.x2)}
  # if(transform.y == "beta"){y = norm2gamma2beta(x = y,skew = skew.y)}



  x1x2<-x1*x2

  #################################


  # output dataframe


  dat<-as.data.frame(cbind(x1,x2,y,x1x2))
  colnames(dat)<-c("x1","x2","y","x1x2")

  return(dat)
}

