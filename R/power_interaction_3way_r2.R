#' Analytic power analysis for 3-way interactions
#'
#' Power analysis for 3-way interaction models, computed via change in R2. Valid for interactions with continuous, normally distributed, variables. Either `b.x1x2x3` or `f2` can be used to specify the magnitude of the interaction effect size.
#'
#' @param N Sample size. Must be a positive integer. Has no default value. Can be a single value or a vector of values.
#' @param b.x1x2x3 Regression coefficient of the 3-way interaction term x1x2x3. Should not be specified if `f2` is specified. Must be between -1 and 1. Default is NULL. Can be a single value or a vector of values.
#' @param r.x1.y Pearson's correlation between x1 and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x2.y Pearson's correlation between x2 and y. Must be between -1 and 1. Assumed to be the 'moderator' in some functions. Has no default value. Can be a single value or a vector of values.
#' @param r.x3.y Pearson's correlation between x3 and y. Must be between -1 and 1. Assumed to be the 'moderator' in some functions. Has no default value. Can be a single value or a vector of values.
#' @param r.x1x2.y Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x1x3.y Pearson's correlation between the interaction term x1x2 (x1 * x3) and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x2x3.y Pearson's correlation between the interaction term x1x2 (x2 * x3) and y. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x1.x2 Pearson's correlation between x1 and x2. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x1.x3 Pearson's correlation between x1 and x3. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param r.x2.x3 Pearson's correlation between x2 and x3. Must be between -1 and 1. Has no default value. Can be a single value or a vector of values.
#' @param rel.x1 Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.x2 Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.x3 Reliability of x3 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.y Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param detailed_results Default is FALSE. Should detailed results be reported? Returns regression slopes, f2, r2, and the full correlation matrix.
#' @param cl Number of clusters to use for running simulations in parallel. Default is NULL (i.e. not in parallel). Useful when running several thousand analyses at once.
#'
#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%do%"
#' @importFrom stats "lm"
#' @importFrom rlang ".data"
#' @importFrom rlang ".env"
#'
#' @return A data frame containing the power for each unique setting combination.
#' @export
#'
#' @examples
#' power_interaction_3way_r2(N=1000,r.x1.y = .1,r.x2.y = .2,r.x3.y = .3,
#' r.x1x2.y =  .05,r.x1x3.y =  .07,r.x2x3.y =  .09,b.x1x2x3 =0.01,
#' r.x1.x2 = .2,r.x1.x3 = .4,r.x2.x3 = .3)
#'
power_interaction_3way_r2<-function(N,
                                    b.x1x2x3,
                                    r.x1.y,
                                    r.x2.y,
                                    r.x3.y,
                                    r.x1x2.y,
                                    r.x1x3.y,
                                    r.x2x3.y,
                                    r.x1.x2,
                                    r.x1.x3,
                                    r.x2.x3,
                                    rel.x1=1,rel.x2=1,rel.x3=1,rel.y=1,
                                    alpha=0.05,detailed_results = FALSE,cl=NULL){


    settings<-expand.grid(list( N=N,
                                b.x1x2x3 = b.x1x2x3,
                                r.y.x1 = r.x1.y,
                                r.y.x2 = r.x2.y,
                                r.y.x3 = r.x3.y,
                                r.y.x1x2 =  r.x1x2.y,
                                r.y.x1x3 =  r.x1x3.y,
                                r.y.x2x3 =  r.x2x3.y,
                                r.x1.x2 = r.x1.x2,
                                r.x1.x3 = r.x1.x3,
                                r.x2.x3 = r.x2.x3,
                                rel.x1=rel.x1,
                                rel.x2=rel.x2,
                                rel.x3=rel.x3,
                                rel.y=rel.y,
                                alpha = alpha))
    if(max(abs(settings$b.x1x2x3))> 1 ){
      stop("b.x1x2x3 must be within [-1,1]")}


  settings$N <- round(settings$N )
  settings = unique(settings)


  if(min(c(settings$rel.x1,settings$rel.x2,settings$rel.x3,settings$rel.y)) <= 0 |
     max(c(settings$rel.x1,settings$rel.x2,settings$rel.x3,settings$rel.y)) > 1 ){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}

  if(max(abs(settings[,grep(x = colnames(settings),pattern = "r[.]")]))> 1 ){
    stop("All correlations must be within [-1,1]")}

  if(min(c(settings$alpha)) < 0 |
     max(c(settings$alpha)) > 1 ){
    stop("alpha must be between 0 and 1, inclusive")}


  message(base::paste("Performing",(base::dim(settings)[1]) ,"analyses",sep=" "))


  adjust.cor = function(cor1,rel1,rel2){
    newcor = cor1*sqrt(rel1*rel2)
    return(newcor)
  }


  two.way.rel = function(cor1,rel1,rel2){
    newrel=((rel1*rel2)+cor1^2)/(1+cor1^2)
    return(newrel)
  }


  # https://stats.stackexchange.com/questions/605858/reliability-of-3-way-interaction-term-between-correlated-variables
  three.way.rel <- function(r.X.M, r.X.Z, r.M.Z, rel.X, rel.M, rel.Z)
  {
    S = matrix(data = c(1,r.X.M,r.X.Z,
                        r.X.M,1,r.M.Z,
                        r.X.Z,r.M.Z,1),nrow = 3,byrow = TRUE)

    (S[1,1]*(S[2,2]*S[3,3] + 2*S[2,3]^2) +
        2*S[1,2]*(S[1,2]*S[3,3] + 2*S[1,3]*S[2,3]) +
        2*S[1,3]*(2*S[1,2]*S[2,3] + S[1,3]*S[2,2])) /
      ((S[1,1]/rel.X)*((S[2,2]/rel.M)*(S[3,3]/rel.Z) + 2*S[2,3]^2) +
         2*S[1,2]*(S[1,2]*(S[3,3]/rel.Z) + 2*S[1,3]*S[2,3]) +
         2*S[1,3]*(2*S[1,2]*S[2,3] + S[1,3]*(S[2,2]/rel.M)))
  }




  if(!is.null(cl)){
    clus <- parallel::makeCluster(cl)
    doParallel::registerDoParallel(clus)
    chunks<-  base::sample(base::rep(c(1:cl),length=base::dim(settings)[1]),replace = FALSE)
    chunks=chunks[base::order(chunks,decreasing = F)]

  }else{
    chunks<-base::rep(1,dim(settings)[1])
  }


  d=NULL

  power_test<-foreach::foreach(d = 1:base::length(unique(chunks)),
                               .combine = 'rbind',
                               .packages = c('dplyr')) %dopar% {


                                 invalid.settings=list()
                                 power.out=list()

                                 full.dat2 = settings[chunks == d,]

                                 for(i in 1: dim(full.dat2)[1]){

                                   settings2 = full.dat2[i,]

                                   cor.mat = base::matrix(data = 0,nrow = 8,ncol = 8) %>% base::as.data.frame()

                                   base::diag(cor.mat) = 1
                                   var.names = c("y","x1","x2","x3","x1x2","x1x3","x2x3","x1x2x3")

                                   base::colnames(cor.mat) = var.names
                                   base::rownames(cor.mat) = var.names

                                   ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )

                                   cor.mat2 = base::data.frame( col = dimnames(cor.mat)[[2]][ind[,2]] ,
                                                                row = dimnames(cor.mat)[[1]][ind[,1]] ,
                                                                val = cor.mat[ ind ] ) %>% as.data.frame()


                                   cor.mat2$V1 = base::paste("r",cor.mat2$col,cor.mat2$row,sep = ".")

                                   rel.mat = data.frame(V1 = var.names,rel=1,sd=1,sd1=1,V2 = paste("rel.",var.names,sep=""))

                                   # bring settings into correlation matrix
                                   cor.mat2$val[!is.na(match(cor.mat2$V1,colnames(settings2)))] = settings2[1,stats::na.omit(match(cor.mat2$V1,colnames(settings2)))]%>% unlist()
                                   rel.mat$rel[!is.na(match(rel.mat$V2,colnames(settings2)))] = settings2[1,stats::na.omit(match(rel.mat$V2,colnames(settings2)))]%>% unlist()

                                   ################# find corresponding r of desired beta

                                   ##### compute standard deviations

                                   rel.mat$sd[which(rel.mat$V1 == "x1x2")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]^2)
                                   rel.mat$sd[which(rel.mat$V1 == "x1x3")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")]^2)
                                   rel.mat$sd[which(rel.mat$V1 == "x2x3")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")]^2)

                                   ## https://doi.org/10.5539/ijsp.v5n6p73

                                   rel.mat$sd[which(rel.mat$V1 == "x1x2x3")] =  sqrt( (1 + 2*cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")]^2) +
                                                                                        (2*cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] * (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] +
                                                                                                                                              (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * 2))) +
                                                                                        (2*cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] +
                                                                                                                                              (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] * 2))) )


                                   #### compute final covariances that are determined
                                   ## https://doi.org/10.5539/ijsp.v5n6p73

                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x2.x1x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")])
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x2.x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")])
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x3.x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")])

                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]) * 2)
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x2.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]) * 2)
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x3.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")]) * 2)


                                   ## convert covariance to correlation

                                   simple.vars = c("r.x1x2.x1x3","r.x1x2.x2x3","r.x1x3.x2x3","r.x1.x1x2x3","r.x2.x1x2x3","r.x3.x1x2x3")

                                   for(a in 1:length(simple.vars)){

                                     var.temp = simple.vars[a]
                                     var.temp2 = base::strsplit(x = var.temp,split = "[.]") %>% unlist()

                                     cov1a = cor.mat2$val[which(cor.mat2$V1 == var.temp)]
                                     sd1a = rel.mat$sd[which(rel.mat$V1 == var.temp2[2])]
                                     sd2a = rel.mat$sd[which(rel.mat$V1 == var.temp2[3])]

                                     cov.mat = base::diag(c(sd1a,sd2a)^2)
                                     cov.mat[cov.mat == 0] =  cov1a


                                     cor.mat2$val[which(cor.mat2$V1 == var.temp)] = stats::cov2cor(cov.mat)[2,1]

                                   }


                                   cor.mat[ind]=cor.mat2$val
                                   ind2=ind
                                   ind2[,1]=ind[,2]
                                   ind2[,2]=ind[,1]
                                   cor.mat[ind2]=cor.mat2$val


                                   int.stats = function(cor.mat,rel.mat,new.cor){
                                     cor.mat[1,8]=new.cor
                                     cor.mat[8,1]=new.cor

                                     ##  path tracing
                                     covmat_all<- base::diag(rel.mat$sd) %*% base::as.matrix(cor.mat) %*% base::diag(rel.mat$sd) # cor 2 cov
                                     b<-covmat_all[1,-1] # the y-row
                                     A<-covmat_all[-1,-1] # everything but the y-row
                                     betas_all<-base::solve(A,b) # get the standardized regression effect sizes


                                     return(data.frame(new.cor = new.cor, beta = betas_all[7]))

                                   }


                                   cor.patterns =  apply(X = matrix(seq(-.01,.01,.005)),1,
                                                         FUN = function(X){int.stats(cor.mat = cor.mat,rel.mat=rel.mat,new.cor = X)},simplify = T) %>%
                                     unlist() %>% matrix(ncol = 2,byrow = T) %>% as.data.frame()

                                   colnames(cor.patterns) = c("new.cor","beta")

                                   breaks=2
                                   while(dim(cor.patterns)[1]<2){
                                     cor.patterns =  apply(X = matrix(seq(-1,1,.1/breaks)),1,
                                                           FUN = function(X){int.stats(cor.mat = cor.mat,rel.mat=rel.mat,new.cor = X)},simplify = T) %>%
                                       unlist() %>% matrix(ncol = 2,byrow = T) %>% as.data.frame()

                                     colnames(cor.patterns) = c("new.cor","beta")

                                     breaks=breaks+1

                                   }


                                     target = settings2$b.x1x2x3

                                     min.row =  which(abs((cor.patterns$beta - target)) == (min(abs(cor.patterns$beta - target))))
                                     if(min.row == 1){second.row = 2}else{
                                       if(min.row == dim(cor.patterns)[1]){second.row = dim(cor.patterns)[1]-1}else{
                                         second.row = c(min.row-1,min.row+1)[which(abs((cor.patterns$beta[c(min.row-1,min.row+1)] - target)) == (min(abs(cor.patterns$beta[c(min.row-1,min.row+1)] - target))))]
                                       }}

                                     out.mat=cor.patterns[c(second.row,min.row),]

                                     f1<-stats::lm( new.cor~beta,data = out.mat )
                                     out.mat=rbind(cor.patterns[min.row,],int.stats(cor.mat = cor.mat,rel.mat=rel.mat,new.cor =  stats::predict(f1,data.frame(beta=target))))

                                     cor.dif = abs(target-out.mat$beta[2])

                                     tol=0.0000001

                                     while(cor.dif >tol){

                                       f1<-stats::lm( new.cor~beta,data = out.mat )
                                       out.mat[1,] = out.mat[2,]
                                       out.mat[2,] = int.stats(cor.mat = cor.mat,rel.mat=rel.mat,new.cor =  stats::predict(f1,data.frame(beta=target)))

                                       cor.dif = abs(target-out.mat$beta[2])
                                     }





                                   cor.mat2$val[which(cor.mat2$V1 == "r.y.x1x2x3")] = out.mat$new.cor[2]

                                   #################
                                   #now do actual power analysis



                                   #### compute 3-way reliability

                                   rel.mat$rel[which(rel.mat$V1 == "x1x2x3")] =  three.way.rel(r.X.M = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")],
                                                                                               r.X.Z = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")],
                                                                                               r.M.Z = cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")],
                                                                                               rel.X = rel.mat$rel[which(rel.mat$V1 == "x1")],
                                                                                               rel.M = rel.mat$rel[which(rel.mat$V1 == "x2")],
                                                                                               rel.Z = rel.mat$rel[which(rel.mat$V1 == "x3")])


                                   #### adjust correlations between y, x1, x2, x3 for reliability
                                   simple.vars = c("r.y.x1","r.y.x2","r.y.x3","r.x1.x2","r.x1.x3","r.x2.x3")
                                   for(a in 1:length(simple.vars)){

                                     var.temp = simple.vars[a]
                                     var.temp2 = base::strsplit(x = var.temp,split = "[.]") %>% unlist()

                                     cor1a = cor.mat2$val[which(cor.mat2$V1 == var.temp)]
                                     rel1a = rel.mat$rel[which(rel.mat$V1 == var.temp2[2])]
                                     rel2a = rel.mat$rel[which(rel.mat$V1 == var.temp2[3])]

                                     cor.mat2$val[which(cor.mat2$V1 == var.temp)] = adjust.cor(cor1 = cor1a,rel1 = rel1a,rel2 = rel2a)

                                   }


                                   #### compute 2-way reliabilities

                                   two.way.vars = c("x1x2","x1x3","x2x3")
                                   for(a in 1:length(two.way.vars)){
                                     var.temp = two.way.vars[a]

                                     var.simp = base::strsplit(x = var.temp,split = "") %>% unlist()
                                     v1 = paste(var.simp[1],var.simp[2],sep = "")
                                     v2 = paste(var.simp[3],var.simp[4],sep = "")

                                     cor1a = cor.mat2$val[which(cor.mat2$V1 == paste("r",v1,v2,sep = "."))]

                                     rel1a = rel.mat$rel[which(rel.mat$V1 == v1)]
                                     rel2a = rel.mat$rel[which(rel.mat$V1 == v2)]

                                     rel.mat$rel[which(rel.mat$V1 == var.temp)] =  two.way.rel(cor1 = cor1a,rel1 = rel1a,rel2 = rel2a)

                                   }


                                   #### adjust correlations between y, x1, x2, x3 and interactions
                                   simple.vars = c("r.y.x1x2","r.y.x1x3","r.y.x2x3","r.y.x1x2x3")
                                   for(a in 1:length(simple.vars)){

                                     var.temp = simple.vars[a]
                                     var.temp2 = base::strsplit(x = var.temp,split = "[.]") %>% unlist()

                                     cor1a = cor.mat2$val[which(cor.mat2$V1 == var.temp)]
                                     rel1a = rel.mat$rel[which(rel.mat$V1 == var.temp2[2])]
                                     rel2a = rel.mat$rel[which(rel.mat$V1 == var.temp2[3])]

                                     cor.mat2$val[which(cor.mat2$V1 == var.temp)] = adjust.cor(cor1 = cor1a,rel1 = rel1a,rel2 = rel2a)

                                   }


                                   ##### compute standard deviations

                                   rel.mat$sd[which(rel.mat$V1 == "x1x2")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]^2)
                                   rel.mat$sd[which(rel.mat$V1 == "x1x3")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")]^2)
                                   rel.mat$sd[which(rel.mat$V1 == "x2x3")] =  sqrt(1 + cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")]^2)

                                   ## https://doi.org/10.5539/ijsp.v5n6p73

                                   rel.mat$sd[which(rel.mat$V1 == "x1x2x3")] =  sqrt( (1 + 2*cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")]^2) +
                                                                                        (2*cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] * (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] +
                                                                                                                                              (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * 2))) +
                                                                                        (2*cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] +
                                                                                                                                              (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] * 2))) )


                                   #### compute final covariances that are determined
                                   ## https://doi.org/10.5539/ijsp.v5n6p73

                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x2.x1x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")])
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x2.x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")])
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1x3.x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] + (cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")])

                                   cor.mat2$val[which(cor.mat2$V1 == "r.x1.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]) * 2)
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x2.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")]) * 2)
                                   cor.mat2$val[which(cor.mat2$V1 == "r.x3.x1x2x3")] = cor.mat2$val[which(cor.mat2$V1 == "r.x1.x2")] + ((cor.mat2$val[which(cor.mat2$V1 == "r.x2.x3")] * cor.mat2$val[which(cor.mat2$V1 == "r.x1.x3")]) * 2)


                                   ## convert covariance to correlation

                                   simple.vars = c("r.x1x2.x1x3","r.x1x2.x2x3","r.x1x3.x2x3","r.x1.x1x2x3","r.x2.x1x2x3","r.x3.x1x2x3")

                                   for(a in 1:length(simple.vars)){

                                     var.temp = simple.vars[a]
                                     var.temp2 = base::strsplit(x = var.temp,split = "[.]") %>% unlist()

                                     cov1a = cor.mat2$val[which(cor.mat2$V1 == var.temp)]
                                     sd1a = rel.mat$sd[which(rel.mat$V1 == var.temp2[2])]
                                     sd2a = rel.mat$sd[which(rel.mat$V1 == var.temp2[3])]

                                     cov.mat = base::diag(c(sd1a,sd2a)^2)
                                     cov.mat[cov.mat == 0] =  cov1a


                                     cor.mat2$val[which(cor.mat2$V1 == var.temp)] = stats::cov2cor(cov.mat)[2,1]

                                   }



                                   ###################

                                   cor.mat[ind]=cor.mat2$val
                                   ind2=ind
                                   ind2[,1]=ind[,2]
                                   ind2[,2]=ind[,1]
                                   cor.mat[ind2]=cor.mat2$val


                                   if(base::min(eigen(cor.mat,only.values = T)$values)<0){
                                     invalid.settings[[length(invalid.settings)+1]]=settings2
                                   }else{


                                     covmat_all<- base::diag(rel.mat$sd1) %*% base::as.matrix(cor.mat) %*% base::diag(rel.mat$sd1) # cor 2 cov

                                     cor.mat.no.int = base::as.matrix(cor.mat)[-base::which(rownames(cor.mat)=="x1x2x3"),-base::which(colnames(cor.mat)=="x1x2x3")]
                                     covmat.no.int <- base::diag(rel.mat$sd1[-base::which(rel.mat$V1 == "x1x2x3")]) %*% cor.mat.no.int %*% base::diag(rel.mat$sd1[-base::which(rel.mat$V1 == "x1x2x3")])# cor 2 cov


                                     ###  path tracing
                                     b<-covmat_all[1,-1] # the y-row
                                     A<-covmat_all[-1,-1] # everything but the y-row

                                     betas_all<-base::solve(A,b) # get the standardized regression effect sizes
                                     totalr2 = base::sum(betas_all*base::as.matrix(cor.mat)[1,-1])
                                     ##
                                     b<-covmat.no.int[1,-1] # the y-row
                                     A<-covmat.no.int[-1,-1] # everything but the y-row

                                     betas_noint<-base::solve(A,b) # get the standardized regression effect sizes
                                     nullr2 = base::sum(betas_noint*base::as.matrix(cor.mat.no.int)[1,-1])

                                     ### power

                                     df_denom <- (settings2$N - (base::dim(rel.mat)[1]-2) )-1
                                     f2 <- (totalr2 -nullr2 )/(1-nullr2)
                                     lambda = f2*df_denom
                                     minusalpha<-1-settings2$alpha
                                     Ft<-stats::qf(minusalpha, 1, df_denom)
                                     Power<-1-stats::pf(Ft, 1, df_denom,lambda)

                                     settings2=cbind(Power,settings2)
                                     colnames(settings2)[1]="pwr"

                                     if(detailed_results == TRUE){

                                       ind <- base::which( base::lower.tri(cor.mat,diag=F) , arr.ind = TRUE )

                                       cor.mat3 = base::data.frame( col = dimnames(cor.mat)[[2]][ind[,2]] ,
                                                                    row = dimnames(cor.mat)[[1]][ind[,1]] ,
                                                                    val = cor.mat[ ind ] ) %>% as.data.frame()

                                       cor.mat3$V1 = base::paste("r",cor.mat3$col,cor.mat3$row,sep = ".")
                                       cor.mat4 = cor.mat3$val %>% t() %>% as.data.frame()
                                       colnames(cor.mat4) = paste("obs.",cor.mat3$V1,sep="")



                                       covmat_all<- base::diag(rel.mat$sd) %*% base::as.matrix(cor.mat) %*% base::diag(rel.mat$sd) # cor 2 cov
                                       ###  path tracing
                                       b<-covmat_all[1,-1] # the y-row
                                       A<-covmat_all[-1,-1] # everything but the y-row
                                       betas_all<-base::solve(A,b) # get the standardized regression effect sizes

                                       details = c(f2,totalr2,nullr2,df_denom,betas_all) %>% t() %>% base::as.data.frame()
                                       colnames(details)=c("f2","totalr2","nullr2","df",base::paste("obs.b",rel.mat$V1[-1],sep=".") )
                                       settings2=base::cbind(settings2,details,cor.mat4)

                                     }


                                     power.out[[base::length(power.out)+1]]=base::as.data.frame(settings2)


                                   }

                                 }

                                 out.power=list()
                                 out.power$power = power.out
                                 out.power$error = invalid.settings
                                 return(out.power)


                               } # end dopar



  if(!is.null(cl)){
    parallel::stopCluster(clus)
    foreach::registerDoSEQ()

    power.out = power_test[,1]
    invalid.settings = power_test[,2]

    ##########

    power.out2=base::list()
    for(X in 1:base::length(power.out)){
      power.out2[[X]] =   base::as.data.frame(base::do.call(base::rbind, power.out[[X]]))}

    power.final = base::as.data.frame(base::do.call(base::rbind, power.out2))



    #########

    remove=0
    for(k in 1:base::length(invalid.settings)){
      if(base::length(invalid.settings[[k]]) == 0){remove=c(remove,k)}
    }
    if(base::length(remove) >1){
      invalid.settings = invalid.settings[-remove]
    }

  }else{
    power.out = power_test$power
    invalid.settings = power_test$error

    power.final = base::as.data.frame(base::do.call(base::rbind, power.out))

  }



  if(length(invalid.settings) > 0){

    invalid2=base::list()
    for(X in 1:base::length(invalid.settings)){
      invalid2[[X]] =   base::as.data.frame(base::do.call(base::rbind, invalid.settings[[X]]))}

    invalid.settings = base::as.data.frame(base::do.call(base::rbind, invalid2))

    warning(base::paste("N=",base::length(invalid.settings)),"power analyses could not be completed due to a correlation matrix that is not postive semi-definite.")
    warning(invalid2)
  }




  return(power.final)
}

