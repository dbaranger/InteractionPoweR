
#' Test interaction
#'
#' Test the interaction from a single simulated data set.
#'
#' @param data Simulated data set. Output of 'generate_interaction()'.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param q Simple slopes. How many quantiles should x2 be split into for simple slope testing? Default is 2.
#'          Simple slope testing returns the effect-size (slope) of y~x1 for the two most extreme quantiles of x2.
#'          If q=3 then the two slopes are y~x1 for the bottom 33% of x2, and the top 33% of x2.
#' @param simple For internal use. Default is FALSE.
#' @return Either a named list or a data frame containing the results of the regression y~x1+x2+x1*x2, the pearson's correlation between y, x1,x2, and x1x2, and the slopes of the simple slopes.
#' @export
#'
#' @examples
#' \dontrun{
#' dataset <- generate_interaction(N = 250,r.x1.y = 0,r.x2.y = .1,r.x1x2.y = -.2,r.x1.x2 = .3)
#' test_interaction(data = dataset, alpha=0.05, q=2)
#' }
test_interaction<-function(data,alpha=0.05,q=2,simple=F){

  #out = list()

  if(length(table(data$y))>2){
    mod<-stats::lm(y ~ x1 + x2 + x1x2, data = data)
  }else{mod<-stats::glm(as.factor(y) ~ x1 + x2 + x1x2, data = data,family = "binomial")
}


  results_out<-stats::coefficients(base::summary(mod))[-1,]
  results<-stats::coefficients(base::summary(mod))[-1,]

  #rownames(results)[3]<-c("x1_x2")
  if(length(table(data$y))>2){colnames(results)<-c("est","se","t","p")}else{colnames(results)<-c("est","se","z","p")}
  results2<-base::as.data.frame(t(c(results[1,],results[2,],results[3,])))
  colnames(results2)<-c(paste(rownames(results)[1], colnames(results),sep="_"),
                        paste(rownames(results)[2], colnames(results),sep="_"),
                        paste(rownames(results)[3], colnames(results),sep="_"))
  rownames(results2)<-NULL
  results2$sig_int <- (results2$x1x2_p < alpha)*1
  # data.correlation<-stats::cor(data)[lower.tri(stats::cor(data))]


  dd = data.frame(cor = as.vector(stats::cor(data)),
                  v1=colnames(data),
                  v2=rep(colnames(data), each=ncol(data)))
  dd = dd[as.vector(upper.tri(stats::cor(data))),]

  dd_out = dd

  cols_dd<-paste("r_",dd$v1,"_",dd$v2,sep="")
  dd<-t(dd[,-c(2,3)])
  dd<-as.data.frame(dd)
  colnames(dd)<-cols_dd

  ## simple slopes
  if(length(table(data$x2)) == 2){

    data$groups<-as.numeric(as.factor(data$x2))
  }else{

  data$groups<-as.numeric(cut(data$x2,
                              include.lowest = T,
                              breaks = stats::quantile(data$x2,probs = seq(0,1,by = 1/q))))
}
  g1<-dplyr::filter(data,data$groups == min(data$groups))
  g2<-dplyr::filter(data,data$groups == max(data$groups))
  if(length(table(data$y))>2){
  mod1<-stats::lm(y~x1,data = g1)
  mod2<-stats::lm(y~x1,data = g2)
  }else{
    mod1<-stats::glm(as.factor(y)~x1,data = g1,family = "binomial")
    mod2<-stats::glm(as.factor(y)~x1,data = g2,family = "binomial")

  }

  slopes = data.frame(lower.slope = stats::coefficients(mod1)[2],
                      upper.slope = stats::coefficients(mod2)[2]
                      )

  results2$est_min<-stats::coefficients(mod1)[2]
  results2$est_max<-stats::coefficients(mod2)[2]


  ####
  results2<-cbind(results2,dd)
  results2<-as.data.frame(results2)
  #results4<-cbind(results2,results3)

  out=list(linear.model = results_out, correlation = dd_out, simple.slopes =slopes )


  if(simple == T){return(results2)}else{return(out)}

}
