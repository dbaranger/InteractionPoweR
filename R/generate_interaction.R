
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
#' @param k.x1 Number of discrete values for x1. Can be used to make a variable binary or ordinal.
#' @param k.x2 Number of discrete values for x2. Can be used to make a variable binary or ordinal.
#' @param k.y Number of discrete values for y.Can be used to make a variable binary or ordinal.
#' @param adjust.correlations If variables are ordinal or binary, should correlations be adjusted so that output data has the specified correlation structure? Default is TRUE.
#' @param tol Correlation adjustment tolerance. When adjust.correlations = TRUE, correlations are adjusted so that the population correlation is within r='tol' of the target. Default = 0.005.
#' @param iter Max number of iterations to run the correlation adjustment for. Typically only a couple are needed. Default = 10.
#' @param N.adjustment Sample size to use when adjusting correlations. Default = 1000000.
#' @param r.x1.y.adjust Internal use only.
#' @param r.x2.y.adjust Internal use only.
#' @param r.x1x2.y.adjust Internal use only.
#' @param r.x1.x2.adjust Internal use only.
#' @param internal.adjust Internal use only.
#' @param skew.x1 No longer supported.
#' @param skew.x2 No longer supported.
#' @param skew.y No longer supported.

#' @return A data frame containing variables 'x1', 'x2', 'y', and 'x1x2'. 'x1x2' is x1*x2. The correlations between these variables are drawn from the defined population-level values. Output variables are all z-scored (mean=0, sd=1).
#' @export
#'
#' @examples
#' dataset <- generate_interaction(N = 10,r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)

generate_interaction <- function(N,
                                     r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,
                                     rel.x1=1,rel.x2=1,rel.y=1,
                                     k.x1 = 0,
                                     k.x2 = 0,
                                     k.y = 0,
                                     adjust.correlations = TRUE,
                                     tol=0.005,iter=10,N.adjustment=1000000,
                                     r.x1.y.adjust = NULL,r.x2.y.adjust=NULL,
                                     r.x1.x2.adjust=NULL,r.x1x2.y.adjust=NULL,
                                     internal.adjust=FALSE,skew.x1=NA,skew.x2=NA,skew.y=NA
) {

  # order of operations
  #  1. Check for input errors
  #  2. Set defaults
  #  3. Adjust correlations for distribution transformations (if requested)
  #  4. Path tracing
  #  5. Simulate data as continuous normal variables
  #  6. Add in noise from reliability
  #  7. Transform variables (number of levels)
  #  8. Output


  if (!is.na(skew.x1)|!is.na(skew.x2)|!is.na(skew.y)){
    stop("Skew is no longer supported")}

  if(adjust.correlations == TRUE){internal.adjust=TRUE}

  # check for input errors
  if(rel.x1 <= 0 | rel.x1 >1 |rel.x2 <= 0 | rel.x2 >1 | rel.y <= 0 | rel.y >1){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}

  if(max(abs(c( r.x1.y,r.x2.y,r.x1.x2,r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1].")}


  if(k.x1 != round(k.x1) | k.x1 == 1 | k.x1 < 0){stop("k.x1 must be a positive integer greater than 1, or 0.")}
  if(k.x2 != round(k.x2)| k.x2 == 1 | k.x2 < 0){stop("k.x2 must be a positive integer greater than 1, or 0.")}
  if(k.y != round(k.y)| k.y == 1 | k.y < 0){stop("k.y must be a positive integer greater than 1, or 0.")}

  # defaults

  sd.x1 = 1
  sd.x2 = 1
  sd.y  = 1
  mean.x1 = 0
  mean.x2 = 0
  mean.y  = 0
  r.x1.x1x2 = 0
  r.x2.x1x2 = 0


  ###### reliability ##########

  rel.x1x2 = ((rel.x1 *   rel.x2) + (  r.x1.x2^2))/(1 + (  r.x1.x2^2))

  obs.r.x1.y = r.x1.y*sqrt(rel.x1*rel.y)
  obs.r.x2.y = r.x2.y*sqrt(rel.x2*rel.y)
  obs.r.x1.x2 = r.x1.x2*sqrt(rel.x1*rel.x2)
  obs.r.x1x2.y = r.x1x2.y*sqrt(rel.x1x2*rel.y)

  r.x1.y   = obs.r.x1.y
  r.x2.y   = obs.r.x2.y
  r.x1.x2  = obs.r.x1.x2
  r.x1x2.y = obs.r.x1x2.y

  ################################

  if(k.x1>0 | k.x2>0 | k.y >0 ){needs.adjustment = TRUE}else{needs.adjustment = FALSE}

  if(adjust.correlations == TRUE & needs.adjustment == TRUE){
    adjustments<-compute_adjustment(r.x1.y = r.x1.y,r.x2.y = r.x2.y,tol = tol,iter = iter,N.adjustment=N.adjustment,
                                    r.x1x2.y = r.x1x2.y,r.x1.x2 = r.x1.x2,

                                    k.x1 = k.x1, k.x2 = k.x2, k.y = k.y
    )

    r.x1.y.adjust = adjustments[2]
    r.x2.y.adjust = adjustments[4]
    r.x1.x2.adjust = adjustments[1]
    r.x1x2.y.adjust = adjustments[6]

  }

  # add adjustment to correlations for subsequent attenuation

  if(!is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust) & !is.null(r.x1.y.adjust)){
    r.x1.y = r.x1.y + r.x1.y.adjust
    r.x2.y = r.x2.y + r.x2.y.adjust
    r.x1.x2 = r.x1.x2 + r.x1.x2.adjust
    r.x1x2.y = r.x1x2.y + r.x1x2.y.adjust

  }


  if(max(abs(c( r.x1.y,r.x2.y,r.x1.x2,r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1].")
  }


  # compute sd of x1x2 (interaction term) using x1 and x2

  cor.mat = matrix(data = c(1,r.x1.x2,
                            r.x1.x2,1),nrow = 2,byrow = TRUE)
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

  if(min(base::eigen(x = cormat,only.values = TRUE)$values) < 0){
    stop("Correlation matrix is not positive semi-definite. Try reducing one or more correlations.")
  }

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

  x1<- base::scale(x = (x1),center = TRUE,scale = TRUE)
  x2<- base::scale(x = (x2),center = TRUE,scale = TRUE)

  if(adjust.correlations == TRUE | internal.adjust == TRUE){

    if(k.x1 >= 2 & k.x2 == 0){
      x1 = norm2ordinal(x = x1,k = k.x1)
      x1<- base::scale(x = (x1),center = TRUE,scale = TRUE)

      xvar = 1 - (r.x1.x2^2)*1
      xsd <- base::sqrt(xvar)
      xmean <- (1 + r.x1.x2 * x1)

      x2 <- stats::rnorm(N, xmean , xsd)
      x2 <- base::scale(x = (x2),center = TRUE,scale = TRUE)


    }
    if(k.x2 >= 2 & k.x1 == 0){
      x2 = norm2ordinal(x = x2,k = k.x2)
      x2<- base::scale(x = (x2),center = TRUE,scale = TRUE)

      xvar = 1 - (r.x1.x2^2)*1
      xsd <- base::sqrt(xvar)
      xmean <- (1 + r.x1.x2 * x2)

      x1 <- stats::rnorm(N, xmean , xsd)
      x1 <- base::scale(x = (x1),center = TRUE,scale = TRUE)

    }

    if(k.x2 >= 2 & k.x1 >= 2){
      x2 = norm2ordinal(x = x2,k = k.x2)
      x2<- base::scale(x = (x2),center = TRUE,scale = TRUE)

      x1 = norm2ordinal(x = x1,k = k.x1)
      x1<- base::scale(x = (x1),center = TRUE,scale = TRUE)

    }

}

  var.y<-covmat[4,4]
  var.x1<-covmat[1,1]
  var.x2<-covmat[2,2]

  # simulate y

  yvar <- var.y -
    ((b1^2)*var.x1) -
    ((b2^2)*var.x2) -
    ((b3^2)*( (var.x1*var.x2) + (covmat[1,2])^2) ) -
    (2*(b1*b2)* covmat[1,2] )

  if(yvar < 0){
    #print()
    stop("Settings produce a negative y-variance. Try reducing one or more correlations.")
  }

  ysd <- base::sqrt(yvar)
  ymean <- (mean.y + b1*x1 + b2*x2 + b3*x1*x2)

  y <- stats::rnorm(N, ymean , ysd)

  y <- base::scale(x = (y),center = TRUE,scale = TRUE)

  if(adjust.correlations == FALSE & internal.adjust == FALSE){

    if(k.x1 >= 2){
      x1 = norm2ordinal(x = x1,k = k.x1)
      x1<- base::scale(x = (x1),center = TRUE,scale = TRUE)
    }
    if(k.x2 >= 2){
      x2 = norm2ordinal(x = x2,k = k.x2)
      x2<- base::scale(x = (x2),center = TRUE,scale = TRUE)
      }

  }

  if(k.y >= 2){
  y = norm2ordinal(x = y,k = k.y)
  y <- base::scale(x = (y),center = TRUE,scale = TRUE)
  }


  x1x2<-x1*x2

  #################################

  # output data frame

  dat<-data.frame(x1=x1,x2=x2,y=y,x1x2=x1x2)
  return(dat)
}

