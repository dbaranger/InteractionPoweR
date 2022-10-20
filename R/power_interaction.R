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
#' @param rel.x1 Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.x2 Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.y Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param k.x1 Number of discrete values for x1. Can be used to make a variable binary or ordinal.
#' @param k.x2 Number of discrete values for x2. Can be used to make a variable binary or ordinal.
#' @param k.y Number of discrete values for y. Can be used to make a variable binary or ordinal.
#' @param adjust.correlations If variables are ordinal or binary, should correlations be adjusted so that output data has the specified correlation structure? Default is TRUE.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param q Simple slopes. How many quantiles should x2 be split into for simple slope testing? Default is 2. Simple slope testing returns the effect-size (slope) of y~x1 for the two most extreme quantiles of x2. If q=3 then the two slopes are y~x1 for the bottom 33% of x2, and the top 33% of x2.
#' @param ss.IQR Simple slope IQR. Multiplier when estimating the distribution of simple slopes within each simulation setting. Default is 1.5.
#' @param cl Number of clusters to use for running simulations in parallel (recommended). Default is 1 (i.e. not in parallel).
#' @param detailed_results Default is FALSE. Should detailed results be reported?
#' @param full_simulation Default is FALSE. If TRUE, will return a list that includes the full per-simulation results.
#' @param tol Correlation adjustment tolerance. When adjust.correlations = TRUE, correlations are adjusted so that the population correlation is within r='tol' of the target. Default = 0.005.
#' @param N.adjustment Sample size for simulations where correlation matrix is corrected to allow for binary/ordinal variables. Default is 1000000
#' @param iter Max number of iterations to run the correlation adjustment for. Typically only a couple are needed. Default = 10.
#'
#' @importFrom dplyr "%>%"
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%do%"
#' @importFrom stats "lm"
#' @importFrom rlang ".data"
#' @importFrom rlang ".env"
#'
#' @return A data frame containing the power (% significant results) for each unique setting combination. If full_simulation = TRUE will return a list, with one data frame that includes power, and a second that includes raw simulation results.
#' @export
#'
#' @examples
#' power_interaction(n.iter=10, N=10,r.x1.y=0.2, r.x2.y=.2,r.x1x2.y=0.5,r.x1.x2=.2)
power_interaction<-function(n.iter,N,r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,
                                 rel.x1=1,rel.x2=1,rel.y=1,
                                 k.x1 = 0,
                                 k.x2 = 0,
                                 k.y = 0,
                                 adjust.correlations = TRUE,
                                 alpha=0.05,q=2,cl=NULL,ss.IQR=1.5,N.adjustment =1000000,
                                 detailed_results=FALSE,full_simulation=FALSE,tol=0.005,iter=10){



  settings<-expand.grid(list( N=N,
                              r.x1.y = r.x1.y,
                              r.x2.y = r.x2.y,
                              r.x1x2.y = r.x1x2.y,
                              r.x1.x2 = r.x1.x2,
                              rel.x1=rel.x1,
                              rel.x2=rel.x2,
                              rel.y=rel.y,
                              k.x1 = k.x1,
                              k.x2 = k.x2,
                              k.y =  k.y,
                              alpha = alpha,
                              q = q))


  settings$N <- round(settings$N )


  if(min(c(settings$rel.x1,settings$rel.x2,settings$rel.y)) <= 0 |
     max(c(settings$rel.x1,settings$rel.x2,settings$rel.y)) > 1 ){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}

  if(max(abs(c( settings$r.x1.y,settings$r.x2.y,settings$r.x1.x2,settings$r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1]")}


  if(any((c(settings$k.x1,settings$k.x2,settings$k.y) - floor(c(settings$k.x1,settings$k.x2,settings$k.y))) != 0 )){
    stop("k.x1,k.x2,and k.y must all be positive integers >=2, or 0")
  }
  if(any(c(settings$k.x1,settings$k.x2,settings$k.y) < 0 )){
    stop("k.x1,k.x2,and k.y must all be positive integers >=2, or 0")
  }


  set2<-unique(settings[,c(2:5)])

  omit<-NA
  i<-NULL

  for(i in 1: dim(set2)[1]){
    out<-base::tryCatch(expr = {
      data<-generate_interaction(N = 10,
                                     adjust.correlations = FALSE,
                                     r.x1.y = set2$r.x1.y[i],
                                     r.x2.y = set2$r.x2.y[i],
                                     r.x1x2.y= set2$r.x1x2.y [i],
                                     r.x1.x2  = set2$r.x1.x2[i]
      )
      data<-0
    },
    warning = function(cond){i},
    error = function(cond){i}
    )
    omit<-c(omit,out)
  }

  omit<-stats::na.omit(omit)
  bad_ones<-omit[omit>0]

  if(length(bad_ones) >0){

    to_be_removed<-set2[bad_ones,]
    to_be_removed$Removed = 1
    settings2<-merge(settings,to_be_removed,all.x=TRUE)
    removed<-settings2[!is.na(settings2$Removed),]

    warning(paste(round(dim(removed)[1]/dim(settings)[1]*100,2)," % of requested simulations are impossible, n=",dim(removed)[1],
                ". Removing from list.",sep=""))
    warning(to_be_removed)

    settings<-settings2[is.na(settings2$Removed),-14]


    if(dim(settings)[1] == 0){
      #print()
      stop("No valid settings")
    }
  }

  if(adjust.correlations == TRUE){internal.adjust=TRUE}

  if(adjust.correlations == TRUE){
    if(max(settings$k.x1) == 0 &  max(settings$k.x2) == 0 & max(settings$k.y) == 0){
      adjust.correlations <- FALSE
    }
  }


  if(adjust.correlations == FALSE){
    settings$r.x1.y.adjust = NULL
    settings$r.x2.y.adjust = NULL
    settings$r.x1.x2.adjust = NULL
    settings$r.x1x2.y.adjust = NULL
  }


  if(!is.null(cl)){
    clus <- parallel::makeCluster(cl)
    doParallel::registerDoParallel(clus)
  }


  if(adjust.correlations == TRUE) {


    settingsa<-base::unique(dplyr::select(settings,dplyr::all_of(c("r.x1.y", "r.x2.y", "r.x1x2.y", "r.x1.x2",
                                                                   "k.x1", "k.x2", "k.y"))))



i=NULL
    new_settings<-foreach::foreach(i = 1: dim(settingsa)[1],.inorder = FALSE,
                                   .combine = 'rbind',
                                   .packages = c('dplyr','MASS'),
                                   .export=c("test_interaction","generate_interaction",
                                             "norm2ordinal",
                                             "compute_adjustment"  )) %dopar% {


                                                 settingsb<-settingsa[i,]


                                                 if(settingsb$k.x1 != 0 | settingsb$k.x2 != 0 | settingsb$k.y != 0){

                                                   adjustments<-base::tryCatch(expr = {
                                                   adjustments<-compute_adjustment(tol = tol,
                                                                                   iter = iter,
                                                                                   N.adjustment = N.adjustment,
                                                                                   r.x1.y = settingsb$r.x1.y,
                                                                                   r.x2.y = settingsb$r.x2.y,
                                                                                   r.x1x2.y = settingsb$r.x1x2.y,
                                                                                   r.x1.x2 = settingsb$r.x1.x2,
                                                                                   k.x1 = settingsb$k.x1,
                                                                                   k.x2 = settingsb$k.x2,
                                                                                   k.y = settingsb$k.y)
                                                   },
                                                   error = function(cond){adjustments<-rep(NA,6) })

                                                 }else{adjustments<-rep(0,6)}

                                                 settingsb$r.x1.y.adjust = adjustments[2]
                                                 settingsb$r.x2.y.adjust = adjustments[4]
                                                 settingsb$r.x1.x2.adjust = adjustments[1]
                                                 settingsb$r.x1x2.y.adjust = adjustments[6]


                                                 return(settingsb)

                                               }



    settings<-base::merge(settings,new_settings,all.x=TRUE)

  }

if( sum(base::is.na(settings$r.x1x2.y.adjust)) > 0){
  error_out = base::paste(sum(base::is.na(settings$r.x1x2.y.adjust))," correlations cannot be adjusted as adjusting results in |r|>1, ",
                    sum(base::is.na(settings$r.x1x2.y.adjust))/dim(settings)[1],"% of settings. These will be removed from simulations" ,sep="")

warning(error_out)
warning(base::unique(settings[base::is.na(settings$r.x1x2.y.adjust),c(1:7)]))
settings = stats::na.omit(settings)

if(dim(settings)[1] == 0){
  #print()
  stop("No valid settings")
}

}




  i<-NULL
  #d<-NULL




  new.order = function(input_dat,col_names){
    col_index = base::match(base::colnames(input_dat),col_names)[!base::is.na(base::match(base::colnames(input_dat),col_names))]

    for(i in 1: base::length(col_index)){
      order1 = base::order(base::as.matrix(input_dat[,col_index[i]]))
      input_dat = input_dat[order1,]
    }

    return(input_dat)
  }


  message(paste("Performing",(dim(settings)[1]*n.iter) ,"simulations",sep=" "))


  if(!is.null(cl)){
    settings$chunk<-  base::sample(base::rep(c(1:cl),length=dim(settings)[1]),replace = FALSE)

  }else{settings$chunk<-1}

  settings_chunks<-base::split(x = settings,f = settings$chunk)


  dimnum<- sapply(X=c(1:dim(settings)[2]), FUN=function(x){length(table(settings[,x]))})
  grouping_variables<-colnames(settings)[dimnum>1]

  adjust_col = base::grep(pattern = "adjust",x = grouping_variables)
  if(base::length(adjust_col) > 0){grouping_variables<-grouping_variables[-adjust_col]}

  adjust_col = base::grep(pattern = "chunk",x = grouping_variables)
  if(base::length(adjust_col) > 0){grouping_variables<-grouping_variables[-adjust_col]}

  quants = c(.025,.5,.975) #quantiles

  d=NULL

  if(dim(settings)[1] == 1){ grouping_variables = "N"}

  power_test<-foreach::foreach(d = 1: length(settings_chunks),
                              .combine = 'rbind',
                              .packages = c('dplyr','MASS'),
                              .export=c("test_interaction","generate_interaction","norm2ordinal" )) %dopar% {


                                          settingsd<-settings_chunks[[d]]

                                          for (i in 1:dim(settingsd)[1] ){
                                            if(adjust.correlations == TRUE){
                                              test_data =  generate_interaction(
                                              adjust.correlations = FALSE,
                                              internal.adjust=TRUE,
                                              N = settingsd$N[i] * n.iter,
                                              r.x1.y = settingsd$r.x1.y[i],
                                              r.x2.y = settingsd$r.x2.y[i],
                                              r.x1x2.y= settingsd$r.x1x2.y [i],
                                              r.x1.x2  = settingsd$r.x1.x2[i],
                                              rel.x1=settingsd$rel.x1[i],
                                              rel.x2=settingsd$rel.x2[i],
                                              rel.y=settingsd$rel.y[i],
                                              k.x1 = settingsd$k.x1[i],
                                              k.x2 = settingsd$k.x2[i],
                                              k.y = settingsd$k.y[i],
                                              r.x1.y.adjust = settingsd$r.x1.y.adjust[i],
                                              r.x2.y.adjust = settingsd$r.x2.y.adjust[i],
                                              r.x1x2.y.adjust = settingsd$r.x1x2.y.adjust[i],
                                              r.x1.x2.adjust = settingsd$r.x1.x2.adjust[i] )
                                            }
                                            if(adjust.correlations == FALSE){
                                              test_data =  generate_interaction(
                                                adjust.correlations = FALSE,
                                                internal.adjust=FALSE,
                                                N = settingsd$N[i] * n.iter,
                                                r.x1.y = settingsd$r.x1.y[i],
                                                r.x2.y = settingsd$r.x2.y[i],
                                                r.x1x2.y= settingsd$r.x1x2.y [i],
                                                r.x1.x2  = settingsd$r.x1.x2[i],
                                                rel.x1=settingsd$rel.x1[i],
                                                rel.x2=settingsd$rel.x2[i],
                                                rel.y=settingsd$rel.y[i],
                                                k.x1 = settingsd$k.x1[i],
                                                k.x2 = settingsd$k.x2[i],
                                                k.y = settingsd$k.y[i],
                                                r.x1.y.adjust = settingsd$r.x1.y.adjust[i],
                                                r.x2.y.adjust = settingsd$r.x2.y.adjust[i],
                                                r.x1x2.y.adjust = settingsd$r.x1x2.y.adjust[i],
                                                r.x1.x2.adjust = settingsd$r.x1.x2.adjust[i] )
                                            }


                                            test_data <- aperm(array(t(test_data),
                                                                     list(4,settingsd$N[i],n.iter)), perm = c(2,1,3))

                                             for(d in 1:n.iter){
                                              a1=test_data[,,d] %>% as.data.frame()
                                              colnames(a1) = c("x1","x2","y","x1x2")

                                              temp = test_interaction( alpha = settingsd$alpha[i],
                                                                       simple = TRUE,
                                                                       detailed_results = detailed_results,
                                                                       data = a1,
                                                                       q = settingsd$q[i])

                                              if(d ==1){
                                                out.f = matrix(data = NA,nrow = n.iter,ncol = dim(temp)[2]) %>% as.data.frame()
                                                colnames(out.f) = colnames(temp)
                                              }

                                              out.f[d,] = temp
                                            }

                                            out.f<-cbind(out.f,settingsd[i,],row.names=NULL)
                                            if(full_simulation == TRUE) {out.f2 = out.f
                                            }else{
                                              settings_keep = settingsd[i,match(x = grouping_variables,colnames(settingsd))] %>% as.data.frame()
                                              colnames(settings_keep) = grouping_variables

                                              power_results<-
                                                out.f %>%
                                                dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                dplyr::summarise(.groups  = "keep",
                                                                 pwr = mean(.data$sig_int ) )

                                              if(detailed_results == TRUE){
                                                power_results3<-
                                                  out.f %>%
                                                  dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                  dplyr::summarise(.groups  = "keep",
                                                                   x1_pwr = mean( .data$x1_p > alpha),
                                                                   x2_pwr = mean( .data$x2_p > alpha)
                                                  )
                                                power_results2<-out.f %>% # effect size
                                                  dplyr::filter(.data$sig_int == 1) %>%
                                                  dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                  dplyr::summarise(.groups = "keep",
                                                                   x1x2_est_mean = mean(.data$x1x2_est),
                                                                   x1x2_r2_mean= mean(.data$x1x2_r2),
                                                                   crossover_mean = mean(.data$crossover),
                                                                   shape_mean = mean(.data$shape),
                                                                   shape_q_2.5 = unname(stats::quantile(.data$shape,quants)[1]),
                                                                   shape_q_97.5 = unname(stats::quantile(.data$shape,quants)[3]),
                                                                   crossover_q_2.5 = unname(stats::quantile(.data$crossover,quants)[1]),
                                                                   crossover_q_97.5 = unname(stats::quantile(.data$crossover,quants)[3]),
                                                                   min.lwr = unname(stats::quantile(.data$est_min)[3]- (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                                                                   min.upr = unname(stats::quantile(.data$est_min)[3]+ (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                                                                   max.lwr = unname(stats::quantile(.data$est_max)[3]- (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                                                                   max.upr = unname(stats::quantile(.data$est_max)[3]+ (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                                                                   x1x2_95_CI_2.5_mean = mean(.data$x1x2_95confint_25 ),
                                                                   x1x2_95_CI_97.5_mean = mean(.data$x1x2_95confint_975 ),
                                                                   x1x2_95_CI_width_mean = mean(.data$x1x2_95confint_975 - .data$x1x2_95confint_25),

                                                                   r_y_x1x2_q_2.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[1]),
                                                                   r_y_x1x2_q_50.0 = unname(stats::quantile(.data$r_y_x1x2,quants)[2]),
                                                                   r_y_x1x2_q_97.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[3]),
                                                                   x1_est_mean   = mean(.data$x1_est),
                                                                   x2_est_mean = mean(.data$x2_est),
                                                                   r_x1_y_mean = mean(.data$r_x1_y),
                                                                   r_x2_y_mean = mean(.data$r_x2_y),
                                                                   r_x1_x2_mean = mean(.data$r_x1_x2),
                                                                   r_y_x1x2_mean = mean(.data$r_y_x1x2),
                                                                   r_x1_x1x2_mean = mean(.data$r_x1_x1x2),
                                                                   r_x2_x1x2_mean = mean(.data$r_x2_x1x2))


                                                power_results4 = merge(power_results3,power_results2)
                                                power_results = merge(power_results,power_results4)
                                              }
                                              out.f2 = merge(settings_keep,power_results)
                                            }

                                            if(i == 1){out.mat = out.f2}else{out.mat = rbind(out.mat,out.f2)}

                                          } # close settingsd loop
                                        return(out.mat)
                                        } # close dopar



  if(full_simulation == TRUE) {

  power_results<-
    power_test %>%
      dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
      dplyr::summarise(.groups  = "keep",
                     pwr = mean(.data$sig_int )
                     )
  if(detailed_results == TRUE){
   power_results3<-
    power_test %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups  = "keep",
                     x1_pwr = mean( .data$x1_p < alpha),
                     x2_pwr = mean( .data$x2_p < alpha)
    )
  power_results2<-power_test %>% # effect size
    dplyr::filter(.data$sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "keep",
                     x1x2_est_mean = mean(.data$x1x2_est),
                     x1x2_r2_mean= mean(.data$x1x2_r2),
                     crossover_mean = mean(.data$crossover),
                     shape_mean = mean(.data$shape),
                     shape_q_2.5 = unname(stats::quantile(.data$shape,quants)[1]),
                     shape_q_97.5 = unname(stats::quantile(.data$shape,quants)[3]),
                     crossover_q_2.5 = unname(stats::quantile(.data$crossover,quants)[1]),
                     crossover_q_97.5 = unname(stats::quantile(.data$crossover,quants)[3]),
                     min.lwr = unname(stats::quantile(.data$est_min)[3]- (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     min.upr = unname(stats::quantile(.data$est_min)[3]+ (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     max.lwr = unname(stats::quantile(.data$est_max)[3]- (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                     max.upr = unname(stats::quantile(.data$est_max)[3]+ (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                     x1x2_95_CI_2.5_mean = mean(.data$x1x2_95confint_25 ),
                     x1x2_95_CI_97.5_mean = mean(.data$x1x2_95confint_975 ),
                     x1x2_95_CI_width_mean = mean(.data$x1x2_95confint_975 - .data$x1x2_95confint_25),

                     r_y_x1x2_q_2.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[1]),
                     r_y_x1x2_q_50.0 = unname(stats::quantile(.data$r_y_x1x2,quants)[2]),
                     r_y_x1x2_q_97.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[3]),
                     x1_est_mean   = mean(.data$x1_est),
                     x2_est_mean = mean(.data$x2_est),
                     r_x1_y_mean = mean(.data$r_x1_y),
                     r_x2_y_mean = mean(.data$r_x2_y),
                     r_x1_x2_mean = mean(.data$r_x1_x2),
                     r_y_x1x2_mean = mean(.data$r_y_x1x2),
                     r_x1_x1x2_mean = mean(.data$r_x1_x1x2),
                     r_x2_x1x2_mean = mean(.data$r_x2_x1x2))

  power_results4 = merge(power_results3,power_results2)

  power_results<-base::merge(power_results,power_results4)
}
  results = list()
  results$results = power_results
  #power_test = power_test %>% dplyr::arrange(.vars = dplyr::vars(dplyr::all_of(grouping_variables)))
  power_test = new.order(input_dat = power_test,col_names = grouping_variables)

  results$simulation  = power_test

  }else{
    #power_test = power_test %>% dplyr::arrange(.vars = dplyr::vars(dplyr::all_of(grouping_variables)))
    power_test = new.order(input_dat = power_test,col_names = grouping_variables)

    results = power_test

    }


  if(!is.null(cl)){
    parallel::stopCluster(clus)
    foreach::registerDoSEQ()
  }



  return(results)
}
