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
#' @param skew.x1 Skew of x1. Default is 0 (normally distributed).
#' @param skew.x2 Skew of x2. Default is 0 (normally distributed).
#' @param skew.y Skew of y. Default is 0 (normally distributed).
#' @param adjust.correlations If variables are skewed or binary, should correlations be adjusted so that output data has the specified correlation structure? Default is TRUE.
#' @param transform.x1 Transform x1? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma" if variables are skewed.
#' @param transform.x2 Transform x2? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma"if variables are skewed.
#' @param transform.y Transform y? Options are "default", "binary", or "gamma". "binary" will cause variable to be binarized  - 2 unique values. Default ("default") will pick "gamma" if variables are skewed.
#' @param alpha The alpha. At what p-value is the interaction deemed significant? Default is 0.05.
#' @param q Simple slopes. How many quantiles should x2 be split into for simple slope testing? Default is 2. Simple slope testing returns the effect-size (slope) of y~x1 for the two most extreme quantiles of x2. If q=3 then the two slopes are y~x1 for the bottom 33% of x2, and the top 33% of x2.
#' @param ss.IQR Simple slope IQR. Multiplier when estimating the distribution of simple slopes within each simulation setting. Default is 1.5.
#' @param cl Number of clusters to use for running simulations in parallel (recommended). Default is 1 (i.e. not in parallel).
#' @param detailed_results Default is FALSE. Should detailed results be reported?
#' @param full_simulation Default is FALSE. If TRUE, will return a list that includes the full per-simulation results.
#' @param seed Simulation seed. Default is NULL, in which case a seed will be chosen at random and echoed to the user. This seed can then be used to repeat the simulation with identical results.
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
#' \dontrun{
#' power_interaction(n.iter=1000, N=seq(100,300,by=10),r.x1.y=0.2, r.x2.y=.2,r.x1x2.y=0.5,r.x1.x2=.2)
#'}
#'
power_interaction<-function(n.iter,N,r.x1.y,r.x2.y,r.x1x2.y,r.x1.x2,
                                 rel.x1=1,rel.x2=1,rel.y=1,
                                 skew.x1 = 0,
                                 skew.x2 = 0,
                                 skew.y = 0,
                                 transform.x1 = "default",
                                 transform.x2 = "default",
                                 transform.y = "default",
                                 adjust.correlations = T,
                                 alpha=0.05,q=2,cl=NULL,ss.IQR=1.5,
                                 detailed_results=FALSE,full_simulation=FALSE,
                                 seed=NULL){

if(is.null(seed)){
  seed =  base::sample(c(1:1000000),1)
  print(paste("Seed is",seed))
}
  base::set.seed(seed = seed)

  settings<-expand.grid(list( N=N,
                              r.x1.y = r.x1.y,
                              r.x2.y = r.x2.y,
                              r.x1x2.y = r.x1x2.y,
                              r.x1.x2 = r.x1.x2,
                              rel.x1=rel.x1,
                              rel.x2=rel.x2,
                              rel.y=rel.y,
                              skew.x1 = skew.x1,
                              skew.x2 = skew.x2,
                              skew.y = skew.y,
                              alpha = alpha,
                              q = q))


  settings$N <- round(settings$N )

  base::print("Checking for errors in inputs...")

  if(min(c(settings$rel.x1,settings$rel.x2,settings$rel.y)) <= 0 |
     max(c(settings$rel.x1,settings$rel.x2,settings$rel.y)) > 1 ){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}

  if(max(abs(c( settings$r.x1.y,settings$r.x2.y,settings$r.x1.x2,settings$r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1]")}


  set2<-unique(settings[,c(2:5)])

  omit<-NA
  i<-NULL

  for(i in 1: dim(set2)[1]){
    out<-base::tryCatch(expr = {
      data<-generate_interaction(N = 10,
                                     adjust.correlations = F,
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
    settings2<-merge(settings,to_be_removed,all.x=T)
    removed<-settings2[!is.na(settings2$Removed),]

    print(paste(round(dim(removed)[1]/dim(settings)[1]*100,2)," % of requested simulations are impossible, n=",dim(removed)[1],
                ". Removing from list.",sep=""))
    print(to_be_removed)

    settings<-settings2[is.na(settings2$Removed),-14]


    if(dim(settings)[1] == 0){
      #print()
      stop("No valid settings")
    }
  }






  if(transform.x1 == "default"){
    if(base::max(settings$skew.x1) != 0){transform.x1 = "gamma"}
    # if(skew.x1 < 0 & skew.x1 >= -1){transform.x1 ="sn"}
   # if(min(settings$skew.x1) < 0){transform.x1 = "beta"}
  }
  if(transform.x2 == "default"){
    if(base::max(settings$skew.x2) != 0){transform.x2 = "gamma"}
    # if(skew.x2 < 0 & skew.x2 >= -1){transform.x2 ="sn"}
  #  if(min(settings$skew.x2) < 0){transform.x2 = "beta"}
  }
  if(transform.y == "default"){
    if(base::max(settings$skew.y) != 0){transform.y = "gamma"}
    #  if(skew.y < 0 & skew.y >= -1){transform.y = "sn"}
   # if(min(settings$skew.y) < 0){transform.y = "beta"}
  }



  if(adjust.correlations == T){
    if(min(settings$skew.x1) == 0 & max(settings$skew.x1) == 0 &
       min(settings$skew.x2) == 0 & max(settings$skew.x2) == 0 &
       min(settings$skew.y) == 0 & max(settings$skew.y)   == 0 &
       transform.x1 == "default" & transform.x2 == "default" & transform.y == "default"){
      adjust.correlations <- F
    }
  }


  if(adjust.correlations == F){
    settings$r.x1.y.adjust = NULL
    settings$r.x2.y.adjust = NULL
    settings$r.x1.x2.adjust = NULL
    settings$r.x1x2.y.adjust = NULL
  }


  if(!is.null(cl)){
    clus <- parallel::makeCluster(cl)
    doParallel::registerDoParallel(clus)
  }


  if(adjust.correlations == T) {

    base::print("Adjusting correlations for variable transformations...")

    settingsa<-base::unique(settings[,c(2:5,9:11)])

    seed.list1 = base::sample(c(1:1000000),dim(settingsa)[1],replace = F)

    # if(!is.null(cl)){settingsa$chunk<-base::rep(c(1:cl),length=dim(settingsa)[1])}else{settingsa$chunk<-1}
    # chunks<-base::split(x = settingsa,f = settingsa$chunk)
i=NULL
    new_settings<-foreach::foreach(i = 1: dim(settingsa)[1],.inorder = F,
                                   .combine = 'rbind',
                                   .packages = c('dplyr','MASS'),
                                   .export=c("test_interaction","generate_interaction",
                                             "norm2binary","norm2gamma",
                                             "compute_adjustment"  )) %dopar% {


                                                base::set.seed(seed.list1[i])
                                                 settingsb<-settingsa[i,]


                                                 if(settingsb$skew.x1 != 0 | settingsb$skew.x2 != 0 | settingsb$skew.y != 0 |
                                                    transform.x1 != "default" | transform.x2 != "default" | transform.y != "default"){

                                                   adjustments<-base::tryCatch(expr = {
                                                   adjustments<-compute_adjustment(r.x1.y = settingsb$r.x1.y,
                                                                                   r.x2.y = settingsb$r.x2.y,
                                                                                   r.x1x2.y = settingsb$r.x1x2.y,
                                                                                   r.x1.x2 = settingsb$r.x1.x2,
                                                                                   skew.x1 = settingsb$skew.x1,
                                                                                   skew.x2 = settingsb$skew.x2,
                                                                                   skew.y = settingsb$skew.y,
                                                                                   transform.x1 = transform.x1,
                                                                                   transform.x2 = transform.x2,
                                                                                   transform.y =  transform.y)
                                                   },
                                                   error = function(cond){adjustments<-rep(NA,6) })

                                                 }else{adjustments<-rep(0,6)}

                                                 settingsb$r.x1.y.adjust = adjustments[2]
                                                 settingsb$r.x2.y.adjust = adjustments[4]
                                                 settingsb$r.x1.x2.adjust = adjustments[1]
                                                 settingsb$r.x1x2.y.adjust = adjustments[6]


                                                 return(settingsb)

                                               }
                                                #)

                                             #   temp_names<-rownames(adjustments_chunk)
                                             #   adjustments_chunk<-as.data.frame(matrix(unlist(adjustments_chunk),
                                             #                                           byrow = T,nrow = dim(settingsc)[1]))
                                             #   colnames(adjustments_chunk)<-temp_names
                                             #
                                             #   return(adjustments_chunk)
                                             #
                                             # }



    settings<-base::merge(settings,new_settings,all.x=T)

  }

if( sum(base::is.na(settings$r.x1x2.y.adjust)) > 0){
  error_out = base::paste(sum(base::is.na(settings$r.x1x2.y.adjust))," correlations cannot be adjusted as adjusting results in |r|>1, ",
                    sum(base::is.na(settings$r.x1x2.y.adjust))/dim(settings)[1],"% of settings. These will be removed from simulations" ,sep="")

print(error_out)
print(base::unique(settings[base::is.na(settings$r.x1x2.y.adjust),c(1:7)]))
settings = stats::na.omit(settings)

if(dim(settings)[1] == 0){
  #print()
  stop("No valid settings")
}

}




  i<-NULL
  #d<-NULL

  print(paste("Performing",(dim(settings)[1]*n.iter) ,"simulations",sep=" "))




  if(dim(settings)[1] == 1 | full_simulation == T){


    seed.list2 = base::sample(c(1:1000000),dim(settings)[1],replace = F)

    power_test<-foreach::foreach(i = 1: dim(settings)[1],
                                 .combine = 'rbind',
                                 .packages = c('dplyr','MASS'),
                                 .export=c("test_interaction","generate_interaction",
                                           "norm2binary","norm2gamma" )) %dopar% {

                                            # settingsd<-settings[d,]
                                             #settingsd<-settingsd[,-18]
                                           #  a = c(1:dim(settingsd)[1])

                                            # out.mat.outer<-sapply(a,FUN = function(a){

                                             #  i = a
                                               base::set.seed(seed.list2[i])

                                             out.mat<-sapply(X = c(1:n.iter),FUN = function(X){
                                                     # i = X
                                                 test_interaction(simple = T,data =
                                                                    generate_interaction(
                                                                      adjust.correlations = F,
                                                                      N = settings$N[i],
                                                                      r.x1.y = settings$r.x1.y[i],
                                                                      r.x2.y = settings$r.x2.y[i],
                                                                      r.x1x2.y= settings$r.x1x2.y [i],
                                                                      r.x1.x2  = settings$r.x1.x2[i],
                                                                      rel.x1=settings$rel.x1[i],
                                                                      rel.x2=settings$rel.x2[i],
                                                                      rel.y=settings$rel.y[i],
                                                                      skew.x1 = settings$skew.x1[i],
                                                                      skew.x2 = settings$skew.x2[i],
                                                                      skew.y = settings$skew.y[i],
                                                                      transform.x1 = transform.x1,
                                                                      transform.x2 = transform.x2,
                                                                      transform.y =  transform.y,
                                                                      r.x1.y.adjust = settings$r.x1.y.adjust[i],
                                                                      r.x2.y.adjust = settings$r.x2.y.adjust[i],
                                                                      r.x1x2.y.adjust = settings$r.x1x2.y.adjust[i],
                                                                      r.x1.x2.adjust = settings$r.x1.x2.adjust[i]
                                                                    ),
                                                                  alpha = settings$alpha[i],
                                                                  q = settings$q[i])



                                               }) # end of inner sapply

                                               out.mat<-t(out.mat)
                                               out.mat<-cbind(out.mat,settings[i,],row.names=NULL)
                                               out.cols<-colnames(out.mat)
                                               out.mat<-as.data.frame(matrix(unlist(out.mat),nrow = n.iter))
                                               colnames(out.mat)<-out.cols

                                               return(out.mat)

                                            # }) # end of outer sapply

                                            # out.mat.outer<-t(out.mat.outer)
                                             #out.mat.outer<-cbind(out.mat.outer,settings[i,],row.names=NULL)
                                           #  out.cols<-colnames(out.mat.outer)
                                            # out.mat.outer<-as.data.frame(matrix(unlist(out.mat.outer),ncol = length(out.cols)))
                                            # colnames(out.mat.outer)<-out.cols

                                             #return(out.mat)

                                           } # end of dopar




#
#   if(dim(settings)[1] == 1){
#
#     # if(!is.null(cl)){
#     #   chunk_settings = data.frame(iter = c(1:n.iter),
#     #                               chunks = rep(1:cl,length=n.iter))
#     # }else{
#     #   chunk_settings = data.frame(iter = c(1:n.iter),
#     #                               chunks = rep(1:1,length=n.iter))
#
#    # }
#
#     out.mat<-foreach::foreach(i = 1:n.iter,.combine = 'rbind',
#                               .packages = c('dplyr','MASS'),
#                               .export=c("test_interaction","generate_interaction",
#                                         "norm2binary","norm2gamma",
#                                         "compute_adjustment")) %do% {
#
#                                          # n_sims = base::length(chunk_settings$iter[chunk_settings$chunks == i])
#
#
#                                         #  simulation <- base::sapply(X = c(1:n_sims),FUN = function(X){
#
#                                             sim<-test_interaction(simple = T,data =
#                                                                     generate_interaction(N = settings$N[1],
#                                                                                              r.x1.y = settings$r.x1.y[1],
#                                                                                              r.x2.y = settings$r.x2.y[1],
#                                                                                              r.x1x2.y= settings$r.x1x2.y [1],
#                                                                                              r.x1.x2  = settings$r.x1.x2[1],
#                                                                                          rel.x1=settings$rel.x1[i],
#                                                                                          rel.x2=settings$rel.x2[i],
#                                                                                          rel.y=settings$rel.y[i],
#                                                                                          skew.x1 = settings$skew.x1[i],
#                                                                                          skew.x2 = settings$skew.x2[i],
#                                                                                          skew.y = settings$skew.y[i],
#                                                                                          transform.x1 = transform.x1,
#                                                                                          transform.x2 = transform.x2,
#                                                                                          transform.y =  transform.y,
#                                                                                          r.x1.y.adjust = settings$r.x1.y.adjust[i],
#                                                                                          r.x2.y.adjust = settings$r.x2.y.adjust[i],
#                                                                                          r.x1x2.y.adjust = settings$r.x1x2.y.adjust[i],
#                                                                                          r.x1.x2.adjust = settings$r.x1.x2.adjust[i]
#
#                                                                     ),
#                                                                   alpha = settings$alpha[1],
#                                                                   q = settings$q[1])
#
#                                            # return(sim)
#                                           #}) # end sapply
#
#
#                                          # simulation<-t(simulation)
#                                           #out.mat.outer<-cbind(out.mat.outer,settings[i,],row.names=NULL)
#                                           #out.cols<-colnames(simulation)
#                                           #simulation<-as.data.frame(matrix(unlist(simulation),ncol = length(out.cols)))
#                                          # colnames(simulation)<-out.cols
#
#                                           return(sim)
#
#                                         } # end do par
#
#     #out.mat<-t(out.mat)
#     out.mat<-cbind(out.mat,settings[1,],row.names=NULL)
#     out.cols<-colnames(out.mat)
#     out.mat<-as.data.frame(matrix(unlist(out.mat),nrow = n.iter))
#     colnames(out.mat)<-out.cols
#     power_test <- out.mat
#
#   } # end if only 1 setting row




 # power_test<-power_test[,is.na(match(colnames(power_test),"chunk"))]

  new_settings<-power_test[,c(which(colnames(power_test) == "N"):dim(power_test)[2])]

  dimnum<- sapply(X=c(1:dim(new_settings)[2]), FUN=function(x){length(table(new_settings[,x]))})
  grouping_variables<-colnames(new_settings)[dimnum>1]

  adjust_col = base::grep(pattern = "adjust",x = grouping_variables)
  if(base::length(adjust_col) > 0){grouping_variables<-grouping_variables[-adjust_col]}

  #grouping_variables<-colnames(results(power_test)[grep("+test", colnames(results(power_test)))])

  power_results<-power_test %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups  = "drop_last",
                     pwr = mean(.data$sig_int ))

  power_results2<-power_test %>% # effect size
    dplyr::filter(.data$sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "drop_last",
                     x1x2_est_mean = mean(.data$x1x2_est),

                     x1x2_r2_mean= mean(.data$x1x2_r2),
                     crossover_mean = mean(.data$crossover),
                     shape_mean = mean(.data$shape),
                     min.lwr = unname(stats::quantile(.data$est_min)[3]- (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     min.upr = unname(stats::quantile(.data$est_min)[3]+ (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                     max.lwr = unname(stats::quantile(.data$est_max)[3]- (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                     max.upr = unname(stats::quantile(.data$est_max)[3]+ (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  )


  power_results3<-power_test %>% # mean effects
    dplyr::filter(.data$sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "drop_last",
                     x1_pwr = mean( .data$x1_p > alpha),
                     x2_pwr = mean( .data$x2_p > alpha)
                     # ,
                     #
                     # r_x1_y_mean = mean(.data$r_x1_y),
                     # r_x2_y_mean = mean(.data$r_x2_y),
                     # r_x1_x2_mean = mean(.data$r_x1_x2),
                     # r_y_x1x2_mean = mean(.data$r_y_x1x2),
                     # r_x1_x1x2_mean = mean(.data$r_x1_x1x2),
                     # r_x2_x1x2_mean = mean(.data$r_x2_x1x2)

    )

  quants = c(.025,.5,.975) #quantiles

  power_results4<-power_test %>% # precision
    dplyr::filter(.data$sig_int == 1) %>%
    dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
    dplyr::summarise(.groups = "drop_last",

                     x1x2_95_CI_2.5_mean = mean(.data$x1x2_95confint_25 ),
                     x1x2_95_CI_97.5_mean = mean(.data$x1x2_95confint_975 ),
                     x1x2_95_CI_width_mean = mean(.data$x1x2_95confint_975 - .data$x1x2_95confint_25),

                     r_y_x1x2_q_2.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[1]),
                     r_y_x1x2_q_50.0 = unname(stats::quantile(.data$r_y_x1x2,quants)[2]),
                     r_y_x1x2_q_97.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[3]),
                     x1_est_mean   = mean(.data$x1_est),
                     x2_est_mean = mean(.data$x2_est),
    )


  results<-power_results

  if(detailed_results == TRUE){results<-merge(results,power_results2,all=T)
                                              results<-merge(results,power_results4,all=T)
                                                 results<-merge(results,power_results3,all=T)}

  pwr_col = base::which(base::colnames(results) == "pwr") - 1
  results = results %>% dplyr::arrange(results[,c(1:pwr_col)])


  if(full_simulation == TRUE){results<-list(results = results, simulation=power_test)  }


  }


  if(dim(settings)[1] > 1 & full_simulation == F){ # start chunked loop


    if(!is.null(cl)){settings$chunk<-base::rep(c(1:cl),length=dim(settings)[1])}else{settings$chunk<-1}
    settings_chunks<-base::split(x = settings,f = settings$chunk)

    seed.list3 = base::sample(c(1:1000000),length(settings_chunks),replace = F)


    d=NULL

    out_final<-foreach::foreach(d = 1: length(settings_chunks),
                                .combine = 'rbind',
                                .packages = c('dplyr','MASS'),
                                .export=c("test_interaction","generate_interaction",
                                          "norm2binary","norm2gamma" )) %dopar% {

                                            base::set.seed(seed.list3[d])

                                            settingsd<-settings_chunks[[d]]
                                            a = c(1:dim(settingsd)[1])


                                            #settingsd<-settingsd[,-18]

                                            out.mat.outer<-sapply(a,FUN = function(a1){

                                              i = a1
                                              out.mat<-sapply(X = c(1:n.iter),FUN = function(X){
                                                # i = X
                                                t1 = test_interaction(simple = T,data =
                                                                        generate_interaction(
                                                                          adjust.correlations = F,
                                                                          N = settingsd$N[i],
                                                                          r.x1.y = settingsd$r.x1.y[i],
                                                                          r.x2.y = settingsd$r.x2.y[i],
                                                                          r.x1x2.y= settingsd$r.x1x2.y [i],
                                                                          r.x1.x2  = settingsd$r.x1.x2[i],
                                                                          rel.x1=settingsd$rel.x1[i],
                                                                          rel.x2=settingsd$rel.x2[i],
                                                                          rel.y=settingsd$rel.y[i],
                                                                          skew.x1 = settingsd$skew.x1[i],
                                                                          skew.x2 = settingsd$skew.x2[i],
                                                                          skew.y = settingsd$skew.y[i],
                                                                          transform.x1 = transform.x1,
                                                                          transform.x2 = transform.x2,
                                                                          transform.y =  transform.y,
                                                                          r.x1.y.adjust = settingsd$r.x1.y.adjust[i],
                                                                          r.x2.y.adjust = settingsd$r.x2.y.adjust[i],
                                                                          r.x1x2.y.adjust = settingsd$r.x1x2.y.adjust[i],
                                                                          r.x1.x2.adjust = settingsd$r.x1.x2.adjust[i]
                                                                        ),
                                                                      alpha = settingsd$alpha[i],
                                                                      q = settingsd$q[i])

                                                return(t1)

                                              }) # end of inner sapply

                                              out.mat<-t(out.mat)
                                              #    out.mat<-cbind(out.mat,settings[i,],row.names=NULL)
                                              out.cols<-colnames(out.mat)
                                              out.mat<-as.data.frame(matrix(unlist(out.mat),nrow = n.iter))
                                              colnames(out.mat)<-out.cols

                                              # return(out.mat)


                                              power_test=out.mat


                                              settings_e<-settingsd[i,]
                                              adjust_col = base::grep(pattern = c("adjust"),x = colnames(settings_e))
                                              chunk_col = base::grep(pattern = c("chunk"),x = colnames(settings_e))
                                              settings_e = settings_e[,-c(adjust_col,chunk_col)]


                                              power_results<-power_test %>%
                                               # dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                dplyr::summarise(.groups  = "drop_last",
                                                                 pwr = mean(.data$sig_int ))

                                              power_results2<-power_test %>%
                                                dplyr::filter(.data$sig_int == 1) %>%
                                               # dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                dplyr::summarise(.groups = "drop_last",
                                                                 x1x2_est_mean = mean(.data$x1x2_est),
                                                                 x1x2_r2_mean= mean(.data$x1x2_r2),
                                                                 crossover_mean = mean(.data$crossover),
                                                                 shape_mean = mean(.data$shape),

                                                                  min.lwr = unname(stats::quantile(.data$est_min)[3]- (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                                                                 min.upr = unname(stats::quantile(.data$est_min)[3]+ (diff(stats::quantile(.data$est_min)[c(2,4)])*ss.IQR))  ,
                                                                 max.lwr = unname(stats::quantile(.data$est_max)[3]- (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  ,
                                                                 max.upr = unname(stats::quantile(.data$est_max)[3]+ (diff(stats::quantile(.data$est_max)[c(2,4)])*ss.IQR))  )


                                              power_results3<-power_test %>%
                                               # dplyr::filter(.data$sig_int == 1) %>%
                                                #dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                dplyr::summarise(.groups = "drop_last",

                                                                 x1_pwr = mean( .data$x1_p < alpha),
                                                                 x2_pwr = mean( .data$x2_p < alpha)


                                                                 # r_x1_y_mean = mean(.data$r_x1_y),
                                                                 # r_x2_y_mean = mean(.data$r_x2_y),
                                                                 # r_x1_x2_mean = mean(.data$r_x1_x2),
                                                                 # r_y_x1x2_mean = mean(.data$r_y_x1x2),
                                                                 # r_x1_x1x2_mean = mean(.data$r_x1_x1x2),
                                                                 # r_x2_x1x2_mean = mean(.data$r_x2_x1x2),

                                                )

                                              quants = c(.025,.5,.975) #quantiles

                                              power_results4<-power_test %>%
                                                dplyr::filter(.data$sig_int == 1) %>%
                                               # dplyr::group_by_at(.vars = dplyr::vars(dplyr::all_of(grouping_variables))) %>%
                                                dplyr::summarise(.groups = "drop_last",

                                                                 x1x2_95_CI_2.5_mean = mean(.data$x1x2_95confint_25 ),
                                                                 x1x2_95_CI_97.5_mean = mean(.data$x1x2_95confint_975 ),
                                                                 x1x2_95_CI_width_mean = mean(.data$x1x2_95confint_975 - .data$x1x2_95confint_25),

                                                                 r_y_x1x2_q_2.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[1]),
                                                                 r_y_x1x2_q_50.0 = unname(stats::quantile(.data$r_y_x1x2,quants)[2]),
                                                                 r_y_x1x2_q_97.5 = unname(stats::quantile(.data$r_y_x1x2,quants)[3]),
                                                                 x1_est_mean   =mean(.data$x1_est),
                                                                 x2_est_mean = mean(.data$x2_est)
                                                )


                                              results<-power_results
                                              results<-merge(results,settings_e,all = T)

                                              if(detailed_results == TRUE){results<-merge(results,power_results2,all=T)
                                                                           results<-merge(results,power_results4,all=T)
                                                                           results<-merge(results,power_results3,all=T)
                                                                         }

                                              if(full_simulation == TRUE){results<-list(results = results, simulation=power_test)  }


                                              return(results)

                                            }) # end of outer sapply

                                            out.mat.outer<-t(out.mat.outer)
                                            #out.mat.outer<-cbind(out.mat.outer,settings[i,],row.names=NULL)
                                            out.cols<-colnames(out.mat.outer)
                                            out.mat.outer<-as.data.frame(matrix(unlist(out.mat.outer),ncol = length(out.cols)))
                                            colnames(out.mat.outer)<-out.cols

                                            return(out.mat.outer)

                                          } # end of dopar


    settings_f = out_final[,c(2:14)]  # ugh hard coded
    group_cols = base::apply(settings_f,2,function(X) base::length(base::table(X)))>1
    settings_e = base::as.data.frame(settings_f[,group_cols])
    colnames(settings_e) = colnames(settings_f)[group_cols]

    #results = base::cbind(settings_e,out_final$pwr,out_final[,c(2:5)])
    results = base::cbind(out_final$pwr,settings_e)
    colnames(results)[1] = "pwr"

    if(detailed_results == TRUE){
      results<-base::cbind(results,
                           out_final[,c(base::which(colnames(out_final)=="x1x2_est_mean"): base::which(colnames(out_final)=="x2_pwr")   )]  )}

    order_cols = stats::na.omit(match(colnames(settings_f),base::colnames(results)))

    # results = base::eval(
    #   base::parse(
    #     text = paste("results[base::order(",paste("results[,",order_cols,"]",collapse = ","), "),]")))

    num_cols = c(1:dim(results)[2])
    results_rest = results[,-c(1,match(order_cols,num_cols))]
    results = cbind(results[,c(order_cols,1)],results_rest)

   pwr_col = base::which(base::colnames(results) == "pwr") - 1

   results = results %>% dplyr::arrange(results[,c(1:pwr_col)])

    }# end chunked loop



  if(!is.null(cl)){
    parallel::stopCluster(clus)
    foreach::registerDoSEQ()
  }

  return(results)






}
