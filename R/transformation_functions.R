#' norm2binary
#'
#' Transforms a vector with a normal distribution to a binomial distribution with two values.
#'
#' @param x Input vector
#' @param skew Desired output skew
#'
#' @return A binary variable
#' @export
#'
#' @examples
#' \dontrun{
#' norm2binary(x = rnorm(n = 100,mean = 0,sd = 1), skew = 1)
#'}
norm2binary = function(x,skew){

  fnToFindRoot = function(s,n,p) {return(((1-2*p)/sqrt(n*p*(1-p)))-s)}
  binomial_p <- stats::uniroot(f = fnToFindRoot, interval = c(0, 1), tol = 0.0001, s=skew,n=1)$root
  p <- stats::pnorm(x, 0, 1)
  p2<-stats::qbinom(p, 1, binomial_p )
  x.new<-base::scale(p2,center = T,scale = T)
  return(x.new)

}





#' norm2likert
#'
#' Transforms a vector with a normal distribution to a binomial distribution with two values.
#'
#' @param x Input vector
#' @param skew Desired output skew
#' @param k Number of discrete values (e.g., 2=binary, 5=likert scale)

#'
#' @return A likert or binary variable
#' @export
#'
#' @examples
#' \dontrun{
#' norm2likert(x = rnorm(n = 100,mean = 0,sd = 1), skew = 1,k=2)
#'}
norm2likert = function(x,skew,k){

  k=k-1

  fnToFindRoot = function(s,n,p) {return(((1-2*p)/sqrt(n*p*(1-p)))-s)}
  binomial_p <- stats::uniroot(f = fnToFindRoot, interval = c(0, 1), tol = 0.0001, s=skew,n=k)$root
  p <- stats::pnorm(x, 0, 1)
  p2<-stats::qbinom(p, k, binomial_p )
  x.new<-base::scale(p2,center = T,scale = T)
  return(x.new)

}





#' binary.p2skew
#'
#' Converts the probability parameter of a binomial distribution to the skew, assuming n=1.
#'
#' @param p The binomial probability
#'
#' @return Skew
#' @export
#'
#' @examples
#'  \dontrun{
#' binary.p2skew(p=.5)
#'}
binary.p2skew = function(p){
  if(min(p) <= 0 |max(p) >= 1){stop("p must be greater than 0 and less than 1")}
  q=1-p
  n=1
  skew = (q-p)/base::sqrt(n*p*q)
  return(skew)
}




#' norm2gamma
#'
#' Transforms a vector with a normal distribution to a gamma distribution.
#'
#' @param x Input vector
#' @param skew Desired skew
#'
#' @return A vector with a (skewed) gamma distribution
#' @export
#'
#' @examples
#' \dontrun{
#' norm2gamma(x = rnorm(n = 100,mean = 0,sd = 1), skew = 1)
#'}
norm2gamma = function(x,skew){
  if(skew ==0){skew<-0.00000001}
  if(skew<0){x=x*-1}
  shape = (2/skew)^2
  p <- stats::pnorm(x, 0, 1)
  p2<- stats::qgamma(p, shape, scale = 1)
  x.new<-base::scale(p2,center = T,scale = T)
  if(skew<0){x.new=x.new*-1}
  return(x.new)
}



#' compute_adjustment
#'
#' Computes how much variable correlations need to be adjusted so that they have the desired correlation structure after transformation. Intended for internal use only.
#'
#' @param r.x1.y Internal use only
#' @param r.x2.y Internal use only
#' @param r.x1x2.y Internal use only
#' @param r.x1.x2 Internal use only
#' @param N.adjustment Internal use only
#' @param tol Internal use only
#' @param iter Internal use only
#' @param skew.x1 Internal use only
#' @param skew.x2 Internal use only
#' @param skew.y Internal use only
#' @param transform.x1 Internal use only
#' @param transform.x2 Internal use only
#' @param transform.y Internal use only
#' @param k.x1 Internal use only
#' @param k.x2 Internal use only
#' @param k.y Internal use only
#'
#' @return Correlation adjustments.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_adjustment(r.x1.y = .2,r.x2.y = .2,r.x1x2.y = .1,r.x1.x2 = .2,transform.y = "binary")
#'}
compute_adjustment<-function(r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,N.adjustment=1000000,tol=0.005,iter=10,
                             skew.x1,
                             skew.x2,
                             skew.y,
                             # levels.x1,
                             # levels.x2,
                             # levels.y,
                             transform.x1,
                             transform.x2,
                             transform.y,
                             k.x1,k.x2,k.y
){

  target_matrix<-c(r.x1.x2,
                   r.x1.y,
                   0,
                   r.x2.y,
                   0,
                   r.x1x2.y)

  adjustments<-c(0,0,0,0,0,0)
  sim_error<-1
  i=1

  while(sim_error>tol & i <iter){


    a2<-generate_interaction(N = N.adjustment,
                                 adjust.correlations = F,
                                 r.x1.y = (r.x1.y+adjustments[2]), # wrong
                                 r.x2.y = (r.x2.y+adjustments[4]), # right
                                 r.x1x2.y = (r.x1x2.y+adjustments[6]), #right
                                 r.x1.x2 = (r.x1.x2+adjustments[1]), # right
                                 transform.x1 = transform.x1,
                                 skew.x1= skew.x1,
                                 transform.x2 = transform.x2,
                                 skew.x2= skew.x2,
                                 transform.y = transform.y,
                                 skew.y= skew.y,
                             k.x1 = k.x1,k.x2=k.x2,k.y=k.y)
    c1<-stats::cor(a2)

    sim_error<-max(abs((target_matrix - c1[lower.tri(c1)]) [-c(3,5) ]))
    m2b<-c1[lower.tri(c1)]
    m1b<-target_matrix
    b<-matrix(t(c(m1b,m2b)),ncol = 2)


    adjust<-apply(b,MARGIN = 1,FUN = function(x){

      d<-base::data.frame(y =c(0,x[1]),x= c(0,x[2]))
      f1<-stats::lm( y~x,data = d )
      p0<-polynom::polynomial(stats::coefficients(f1))
      return(stats::predict(p0,x[1]))
    })

    adjustments<-adjustments + (adjust - b[,1])
    i=i+1
  }
  return(adjustments)
}


