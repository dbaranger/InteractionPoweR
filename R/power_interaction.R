#' Power analysis for interactions
#'
#' Power analysis for interaction models, by simulation. A set of n.iter simulations is run for each unique combination of model settings.
#'
#' @param n.iter Number of iterations. The number of simulations to run for each unique setting combination. Must be a positive integer.
#' @param N Sample size. Must be a positive integer. Has no default value. Can be a single value or a vector of values.
#' @param r.x1.y Pearson's correlation between x1 and y. Must be between -1 and 1.. Has no default value. Can be a single value or a vector of values.
#' @param r.x2.y Pearson's correlation between x2 and y. Must be between -1 and 1.. Assumed to be the 'moderator' in some functions. Has no default value. Can be a single value or a vector of values.
#' @param r.x1x2.y Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1.. Has no default value. Can be a single value or a vector of values.
#' @param r.x1.x2 Pearson's correlation between x1 and x2. Must be between -1 and 1.. Has no default value. Can be a single value or a vector of values.
#' @param rel_x1 Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel_x2 Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel_y Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param q Simple slopes. How many quantiles should x2 be split into for simple slope testing? Default is 2. Simple slope testing returns the effect-size (slope) of y~x1 for the two most extreme quantiles of x2. If q=3 then the two slopes are y~x1 for the bottom 33% of x2, and the top 33% of x2.
#' @param ss.IQR Simple slope IQR. Multiplier when estimating the distribution of simple slopes within each simulation setting. Default is 1.5.
#' @param cl Number of clusters to use for running simulations in parallel (recommended). Default is 1 (i.e. not in parallel).
#' @param detailed_results Default is FALSE. Should detailed results be reported (mean of effect size and correlations for each setting)?
#' @param full_simulation Default is FALSE. If TRUE, will return a list that includes the full per-simulation results.
#'
#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
#' @importFrom stats "lm"
#' @importFrom rlang ".data"
#' @importFrom rlang ".env"
#'
#' @return A data frame containing the power (% significant results) for each unique setting combination. If full_simulation = TRUE will return a list, with one data frame that includes power, and a second that includes raw simulation results.
#' @export
#'
#' @examples
#' \dontrun{
#' power_interaction(n.iter=1000, N=seq(100,300,by=10),r.x1.y=0.2, r.x2.y=.2,r.x1x2.y=0.5,r.x1.x2=.2)
#'}
#'
power_interaction<-function(n.iter,N,r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,rel_x1=1,rel_x2=1,rel_y=1,alpha=0.05,q=2,cl=NULL,ss.IQR=1.5,detailed_results=FALSE,full_simulation=FALSE){

  settings<-expand.grid(list( N=N,
                              r.x1.y = r.x1.y,
                              r.x2.y = r.x2.y,
                              r.x1x2.y = r.x1x2.y,
                              r.x1.x2 = r.x1.x2,
                              rel_x1=rel_x1,
                              rel_x2=rel_x2,
                              rel_y=rel_y,
                              # sd.x1 = sd.x1,
                              # sd.x2 = sd.x2,
                              # sd.y = sd.y,
                              # mean.x1 = mean.x1,
                              # mean.x2 = mean.x2,
                              # mean.y = mean.y,
                              alpha = alpha,
                              q = q))


  omit<-NA
  i<-NULL
  for(i in 1: dim(settings)[1]){
    out<-base::tryCatch(expr = {
      data<-generate_interaction(N = settings$N[i],
                                 r.x1.y = settings$r.x1.y[i],
                                 r.x2.y = settings$r.x2.y[i],
                                 r.x1x2.y= settings$r.x1x2.y [i],
                                 r.x1.x2  = settings$r.x1.x2[i],
                                 rel_x1=settings$rel_x1[i],
                                 rel_x2=settings$rel_x2[i],
                                 rel_y=settings$rel_y[i]
                                 # mean.x1 = settings$mean.x1[i],
                                 # sd.x1 = settings$sd.x1[i],
                                 # mean.x2 = settings$mean.x2[i],
                                 # sd.x2 = settings$sd.x2[i],
                                 # mean.y = settings$mean.y[i],
                                 # sd.y = settings$sd.y[i]
                                 #sd.x1x2 = settings$sd.x1x2[i],
                                 #r.x1.x1x2=settings$r.x1.x1x2[i],
                                 #r.x2.x1x2=settings$r.x2.x1x2[i],
                                 #binary.x1 = settings$ binary.x1[i],
                                 #binary.x2 = settings$ binary.x2[i]

      )

    },
    warning = function(cond){i},
    error = function(cond){i}
    )
    omit<-c(omit,out)
  }


    l<-base::lapply(omit, base::is.integer)
    bad_ones<-stats::na.omit(base::unlist(omit[l==TRUE]))
if(length(bad_ones) >0){
    print(paste(round(length(bad_ones)/dim(settings)[1]*100,2)," % of requested simulations are invalid, n=",length(bad_ones),
                ". Removing from list.",sep=""))

    settings<-settings[-bad_ones,]


    if(dim(settings)[1] == 0){
      #print()
      stop("No valid settings")
    }
}

  print(paste("Performing",(dim(settings)[1]*n.iter) ,"tests",sep=" "))
  #
  # defaults<-data.frame(alpha=0.05,q=2,#binary.x1=0,binary.x2=0,
  #                      r.x1.x1x2=0,r.x2.x1x2=0,
  #                      sd.x1=1,sd.x2=1,sd.y=1,sd.x1x2=1,
  #                      mean.x1=0,mean.x2=0,mean.y=0)

  # drops<-intersect(colnames(settings),colnames(defaults))
  # if(length(drops)>0){defaults<-defaults[,-which(colnames(defaults) == drops)]}
  # defaults<-as.list(defaults)
  # settings<-expand.grid(c(params,defaults))
  settings$N <- round(settings$N )


  if(!is.null(cl)){
  clus <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(clus)
  }
i<-NULL
    if(dim(settings)[1] > 1){

    power_test<-foreach::foreach(i = 1: dim(settings)[1],
                        .combine = 'rbind',
                        .packages = c('dplyr','MASS'),
                        .export=c("test_interaction","generate_interaction")) %dopar% {
                          out.mat<-sapply(X = c(1:n.iter),FUN = function(X){test_interaction(data =
                                                                                               generate_interaction(N = settings$N[i],
                                                                                                                    r.x1.y = settings$r.x1.y[i],
                                                                                                                    r.x2.y = settings$r.x2.y[i],
                                                                                                                    r.x1x2.y= settings$r.x1x2.y [i],
                                                                                                                    r.x1.x2  = settings$r.x1.x2[i],
                                                                                                                    rel_x1=settings$rel_x1[i],
                                                                                                                    rel_x2=settings$rel_x2[i],
                                                                                                                    rel_y=settings$rel_y[i]
                                                                                                                    # mean.x1 = settings$mean.x1[i],
                                                                                                                    # sd.x1 = settings$sd.x1[i],
                                                                                                                    # mean.x2 = settings$mean.x2[i],
                                                                                                                    # sd.x2 = settings$sd.x2[i],
                                                                                                                    # mean.y = settings$mean.y[i],
                                                                                                                    # sd.y = settings$sd.y[i]
                                                                                                                    #sd.x1x2 = settings$sd.x1x2[i],
                                                                                                                    #r.x1.x1x2=settings$r.x1.x1x2[i],
                                                                                                                    #r.x2.x1x2=settings$r.x2.x1x2[i],
                                                                                                                    #binary.x1 = settings$ binary.x1[i],
                                                                                                                    #binary.x2 = settings$ binary.x2[i]
                                                                                               ),
                                                                                             alpha = settings$alpha[i],
                                                                                             q = settings$q[i])})

                          out.mat<-t(out.mat)
                          out.mat<-cbind(out.mat,settings[i,],row.names=NULL)
                          out.cols<-colnames(out.mat)
                          out.mat<-as.data.frame(matrix(unlist(out.mat),nrow = n.iter))
                          colnames(out.mat)<-out.cols





                          return(out.mat)
                        }

  }


  if(dim(settings)[1] == 1){

    out.mat<-foreach::foreach(i = 1:n.iter,.combine = 'rbind',
                     .packages = c('dplyr','MASS'),
                     .export=c("test_interaction","generate_interaction")) %dopar% {
                       simulation <- test_interaction(data =
                                                        generate_interaction(N = settings$N[1],
                                                                             r.x1.y = settings$r.x1.y[1],
                                                                             r.x2.y = settings$r.x2.y[1],
                                                                             r.x1x2.y= settings$r.x1x2.y [1],
                                                                             r.x1.x2  = settings$r.x1.x2[1],
                                                                             rel_x1=settings$rel_x1[1],
                                                                             rel_x2=settings$rel_x2[1],
                                                                             rel_y=settings$rel_y[1]
                                                                             # mean.x1 = settings$mean.x1[1],
                                                                             # sd.x1 = settings$sd.x1[1],
                                                                             # mean.x2 = settings$mean.x2[1],
                                                                             # sd.x2 = settings$sd.x2[1],
                                                                             # mean.y = settings$mean.y[1],
                                                                             # sd.y = settings$sd.y[1]
                                                                             #sd.x1x2 = settings$sd.x1x2[1],
                                                                             #r.x1.x1x2=settings$r.x1.x1x2[1],
                                                                             #r.x2.x1x2=settings$r.x2.x1x2[1],
                                                                             #binary.x1 = settings$ binary.x1[1],
                                                                             #binary.x2 = settings$ binary.x2[1]
                                                                             ),
                                                      alpha = settings$alpha[1],
                                                      q = settings$q[1])

                       return(simulation)
                     }

    #out.mat<-t(out.mat)
    out.mat<-cbind(out.mat,settings[1,],row.names=NULL)
    out.cols<-colnames(out.mat)
    out.mat<-as.data.frame(matrix(unlist(out.mat),nrow = n.iter))
    colnames(out.mat)<-out.cols
    power_test <- out.mat

  }





  if(!is.null(cl)){
  parallel::stopCluster(clus)
  foreach::registerDoSEQ()
  }


  settings<-power_test[,c(22:29)]
  dimnum<- sapply(X=c(1:dim(settings)[2]), FUN=function(x){length(table(settings[,x]))})
  grouping_variables<-colnames(settings)[dimnum>1]
  #grouping_variables<-colnames(results(power_test)[grep("+test", colnames(results(power_test)))])

  power_results<-power_test %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups  = "drop_last",
                     pwr = mean(.data$sig_int ))

  power_results2<-power_test %>%
    dplyr::filter(.data$sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "drop_last",
                     min.lwr = unname(stats::quantile(.data$est_min)[3]- (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     min.upr = unname(stats::quantile(.data$est_min)[3]+ (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     max.lwr = unname(stats::quantile(.data$est_max)[3]- (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                     max.upr = unname(stats::quantile(.data$est_max)[3]+ (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  )

  power_results3<-power_test %>%
    # dplyr::filter(sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "drop_last",
                     x1_est_mean   =mean(.data$x1_est),
                     x2_est_mean = mean(.data$x2_est),
                     x1x2_est_mean = mean(.data$x1x2_est),
                     r_x1_y_mean = mean(.data$r_x1_y),
                     r_x2_y_mean = mean(.data$r_x2_y),
                     r_x1_x2_mean = mean(.data$r_x1_x2),
                     r_y_x1x2_mean = mean(.data$r_y_x1x2),
                     r_x1_x1x2_mean = mean(.data$r_x1_x1x2),
                     r_x2_x1x2_mean = mean(.data$r_x2_x1x2)


    )


  results<-merge(power_results,power_results2,all = T)

  if(detailed_results == TRUE){results<-merge(results,power_results3,all=T)}

  if(full_simulation == TRUE){results<-list(results = results, simulation=power_test)  }


  return(results)

}
