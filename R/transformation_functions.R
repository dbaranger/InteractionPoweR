#' norm2ordinal
#'
#' Transforms a vector with a normal distribution to a binomial distribution with two values.
#'
#' @param x Input vector
#' @param k Number of discrete values (e.g., 2=binary, 5=likert scale)

#'
#' @return A ordinal or binary variable
#' @export
#'
#' @examples
#' norm2ordinal(x = rnorm(n = 100,mean = 0,sd = 1),k=2)
norm2ordinal = function(x,k){
  skew=0
  k=k-1

  fnToFindRoot = function(s,n,p) {return(((1-2*p)/sqrt(n*p*(1-p)))-s)}
  binomial_p <- stats::uniroot(f = fnToFindRoot, interval = c(0, 1), tol = 0.0001, s=skew,n=k)$root
  p <- stats::pnorm(x, 0, 1)
  p2<-stats::qbinom(p, k, binomial_p )
  x.new<-base::scale(p2,center = T,scale = T)
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
#' @param k.x1 Internal use only
#' @param k.x2 Internal use only
#' @param k.y Internal use only
#'
#' @return Correlation adjustments.
#' @export
#'
#' @examples
#' \donttest{
#' compute_adjustment(r.x1.y = .2,r.x2.y = .2,r.x1x2.y = .1,r.x1.x2 = .2,
#'k.x1 = 0,k.x2=0,k.y=2)
#'}
compute_adjustment<-function(r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,
                             N.adjustment=1000000,tol=0.005,iter=10,
                             k.x1,
                             k.x2,
                             k.y
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
                                 adjust.correlations = FALSE,
                                 r.x1.y = (r.x1.y+adjustments[2]),
                                 r.x2.y = (r.x2.y+adjustments[4]),
                                 r.x1x2.y = (r.x1x2.y+adjustments[6]),
                                 r.x1.x2 = (r.x1.x2+adjustments[1]),

                             internal.adjust=TRUE,
                             k.x1 = k.x1,k.x2=k.x2,k.y=k.y)
    c1<-stats::cor(a2)

    sim_error<-max(abs((target_matrix - c1[lower.tri(c1)]) [-c(3,5) ]))


    m1b<-target_matrix
    m2b<-c1[lower.tri(c1)]
    m3b = c((r.x1.x2+adjustments[1]),
            (r.x1.y+adjustments[2]),
            0,
            (r.x2.y+adjustments[4]),
            0,
            (r.x1x2.y+adjustments[6]))

    b<-matrix(t(c(m1b,m2b,m3b)),ncol = 3) # target, observed, settings


    adjust<-apply(b,MARGIN = 1,FUN = function(x){

      d<-base::data.frame(y =c(0,x[3]),x= c(0,x[2]))
      f1<-stats::lm( y~x,data = d )
      p0<-polynom::polynomial(stats::coefficients(f1))
      return(stats::predict(p0,x[1]))
    })



    cor.mat = base::matrix(data = 0,nrow = 4,ncol = 4) %>% base::as.data.frame()
    base::diag(cor.mat) = 1
    ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )
    cor.mat[ind]=adjust
    ind2=ind
    ind2[,1]=ind[,2]
    ind2[,2]=ind[,1]
    cor.mat[ind2]=adjust

    while(min(eigen(cor.mat,only.values = T)$values)<0){
      cor.pd = Matrix::nearPD(x = as.matrix(cor.mat),keepDiag = T,corr = T,do2eigen = T)$mat
      cor.pd = cor.pd[base::lower.tri(cor.pd,diag=F)]

      max.diff = which(abs(cor.pd - adjust)[-c(3,5)] == max(abs(cor.pd - adjust)[-c(3,5)]))

      adjust[-c(3,5)] [max.diff] = cor.pd[-c(3,5)] [max.diff]

      cor.mat = base::matrix(data = 0,nrow = 4,ncol = 4) %>% base::as.data.frame()
      base::diag(cor.mat) = 1
      ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )
      cor.mat[ind]=adjust
      ind2=ind
      ind2[,1]=ind[,2]
      ind2[,2]=ind[,1]
      cor.mat[ind2]=adjust
      }

      adjustments<-(adjust - b[,1])


    if(max(adjustments+target_matrix)>1){

          stop("Data cannot be transformed to desired distribution")
    }

    i=i+1

    }

  if(sim_error>tol ){
    stop(paste("Data cannot be transformed to desired distribution. Error (",round(sim_error,4),") is greater than tol (",tol,") Try reducing correlations, increasing tolerance, or increasing the number of iterations.",sep=""))

  }

  return(adjustments)
}


