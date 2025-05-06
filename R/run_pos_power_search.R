#' Calculate the point of stability for a two-way interaction estimate
#’
#’ This function conducts a simulation-based binary search to find the minimum sample size (POS) at which a proportion of estimates for the interaction coefficient in (y ~ x1 + x2 + x1*x2) falls within a specified interval (COS) around the expected value.
#’

#' @param r.x1.y Pearson's correlation between x1 and y. Must be between -1 and 1. Has no default value.
#' @param r.x2.y Pearson's correlation between x2 and y. Must be between -1 and 1. Assumed to be the 'moderator' in some functions. Has no default value.
#' @param r.x1x2.y Pearson's correlation between the interaction term x1x2 (x1 * x2) and y. Must be between -1 and 1. If NULL or unspecified, it is automatically assigned as the average of the main effects divided by 2, considering reliability. The default value is NULL.
#' @param r.x1.x2 Pearson's correlation between x1 and x2. Must be between -1 and 1. Has no default value.
#' @param rel.x1 Reliability of x1 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.x2 Reliability of x2 (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param rel.y Reliability of xy (e.g. test-retest reliability, ICC, Cronbach's alpha). Default is 1 (perfect reliability). Must be greater than 0 and less than or equal to 1.
#' @param single.N Numeric or NULL. A specific sample size at which the stability statistics (POS, COS) are computed. A numeric value overrides the search for a POS. Must be a whole number greater than 0 or NULL. Default is NULL.
#' @param step Numeric or NULL. Determines the rounding (degree of precision) in the estimate of the POS. When set to NULL, dynamic adjustments are made with 1/100th the of the interval between the upper.bound and lower.bound on each iteration of the search algorithm. If it is non-null, it must be a whole natural number. Default is NULL.
#' @param start.power Numeric or NULL. Determines the starting sample size of the search by locating the sample size with power equal to start.power. If it is non-null, it must be greater than 0 and less than 1. Note that a non-null start.power overrides user inputs for the parameters lower.bound and upper.bound so both are NULL. The default is 0.8.
#' @param n.datasets Numeric. Number of simulations conducted at each sample size identified by the search process. Larger values result in greater computational demands. Default is 10000. Must be greater than 1.
#' @param lower.bound Numeric or NULL. The lower bound for the search. It is automatically assigned when set to NULL The user can manually input one. There is an automatic override to NULL if start.power is non-null. The default is NULL.
#' @param upper.bound Numeric or NULL. The upper bound for the search. If it is null, an upper bound is determined automatically; there is a ceiling of 1e7. Default is NULL.
#' @param cos.width Numeric. The width of the corridor of stability COS within which we require some number of estimates to fall to identify the POS. Must be greater than 0. Values above 1 are unlikely to be informative. Default is 0.5.
#' @param pos.percent Numeric. The minimum proportion of simulated estimates we require to fall inside the COS to consider the sample size to produce stable estimates. Default is 0.8.
#' @param n.cores Numeric. Number of cores to be used in parallel processing. Default is 1.
#'
#' @return A data frame containing the variables 'rel.x1', 'rel.x2', 'rel.x1x2', 'rel.y', 'r.x1x2.y', 'r.x1.y', 'r.x2.y', 'obs.r.x1x2.y', 'COS_pos.percents', 'COS_interval', 'POS_samplesize', and 'POS_power'
#'         OR a data frame containing the variables 'single.N', 'rel.x1', 'rel.x2', 'rel.x1x2', 'rel.y', 'r.x1x2.y', 'r.x1.y', 'r.x2.y', 'obs.r.x1x2.y', 'COS_pos.percents', 'within_COS_interval', 'POS_determined_COS', 'within_POS_determined_COS', 'POS_determined_COS.percent', sign_error_rate', and 'power'
#’
#'         If the user inputs single.N as NULL, the first dataframe is returned with the search results for the POS.
#'         If the user inputs single.N as a numeric whole number greater than 3, the second dataframe is returned with stability information at that specific sample size.
#'
#'         COS_pos.percents: A string representing the COS width and target POS percentage.
#’         COS_interval: The COS interval around the expected interaction effect.
#’         POS_samplesize: The sample size where interaction estimates are regarded as stable.
#’         POS_power: The statistical power, calculated analytically, at the POS.
#'         within_COS_interval: The proportion of estimates within the inputted COS at the sample size single.N.
#'         POS_determined_COS: The COS necessary to achieve the inputted pos.percent at the sample size single.N. There is a ceiling of 100 times the magnitude of the estimate.
#'         within_POS_determined_COS: The proportion of estimates within the POS_determined_COS at the sample size single.N.
#'         POS_determined_COS.percent: The percent deviation from the expected value defining POS_determined_COS.
#'         sign_error_rate: The number of estimates signed incorrectly at the sample size single.N.
#'         power: The statistical power, calculated analytically, at the sample size single.N.
#'         r.x1x2.y_to_be_stable_at_single.N: The minimum interaction effect required for stability at the sample size single.N, determined via binary search.
#'         power_to_be_stable: The statistical power associated with the stable_effect_size at sample size single.N.
#’
#' @export
#'
#' @examples
#' run_pos_power_search(r.x1.y = 0.2, r.x2.y = 0.2, r.x1x2.y = 0.15, r.x1.x2 = 0.1,
#' rel.x1 = 0.8, rel.x2 = 0.8,start.power = 0.8, step = NULL, n.datasets = 1000,
#' lower.bound = NULL, upper.bound = 500, cos.width = 0.5, pos.percent = 0.8)

run_pos_power_search <- function(r.x1.y,
                                 r.x2.y,
                                 r.x1x2.y = NULL,
                                 r.x1.x2,
                                 rel.x1,
                                 rel.x2,
                                 rel.y = 1,
                                 single.N = NULL,
                                 start.power = 0.8,
                                 step = NULL,
                                 n.datasets = 1000,
                                 lower.bound = NULL,
                                 upper.bound = NULL,
                                 cos.width = 0.5,
                                 pos.percent = 0.8,
                                 n.cores = 1) {


  if(length(r.x1.y) != 1) stop("r.x1.y must be a single value")
  if(length(r.x2.y) != 1) stop("r.x2.y must be a single value")
  if(length(r.x1.x2) != 1) stop("r.x1.x2 must be a single value")
  if(length(rel.x1) != 1) stop("rel.x1 must be a single value")
  if(length(rel.x2) != 1) stop("rel.x2 must be a single value")
  if(length(rel.y) != 1) stop("rel.y must be a single value")
  if(length(n.datasets) != 1) stop("n.datasets must be a single value")
  if(length(cos.width) != 1) stop("cos.width must be a single value")
  if(length(pos.percent) != 1) stop("pos.percent must be a single value")

  # check for input errors
  if(rel.x1 <= 0 | rel.x1 >1 |rel.x2 <= 0 | rel.x2 >1 | rel.y <= 0 | rel.y >1){
    stop("All reliabilities must be greater than 0 and less than or equal to 1")}

  if(max(abs(c( r.x1.y,r.x2.y,r.x1.x2,r.x1x2.y)))> 1 ){
    stop("All correlations must be within [-1,1].")}

  if(cos.width <= 0 | pos.percent >1 | pos.percent <= 0){
    stop("COS width must be greater than 0.")}

  if(pos.percent >1 | pos.percent <= 0){
    stop("POS percent must be within (0,1]")}

  if (!is.null(step)) {
    if(step <= 0 | step != round(step)){
      stop("Step must be a whole number greater than 0.")}
  }

  if (!is.null(single.N)) {
    if(single.N <= 0 | single.N != round(single.N)){
      stop("single.N must be a whole number greater than 0 or NULL.")}
  }

  if(n.datasets <= 0 | n.datasets != round(n.datasets)){
    stop("n.datasets must be a whole number greater than 0.")}

  if (!is.null(start.power)) {
    if(start.power <= 0 | start.power >= 1){
      stop("start.power must be in (0, 1) or NULL.")}
  }

  if (!is.null(upper.bound) & !is.null(lower.bound)) {
    if (lower.bound >= upper.bound | lower.bound < 0 | upper.bound < 0) {
      stop("lower.bound must be a whole number greater than 0 and less than upper.bound, or NULL")
    }
    if (lower.bound != round(lower.bound) | upper.bound != round(upper.bound)) {
      stop("lower.bound and upper.bound must be whole numbers greater than 0, or NULL.")
    }
    if (lower.bound < 0 | lower.bound != round(lower.bound)) {
      stop("lower.bound must be a whole number greater than 0.")
    }
  }


  #If
  if (is.null(r.x1x2.y)) {
    message("r.x1x2.y was not specified, and has been automatically set as the average of the main effects divided by 2.")
    r.x1x2.y <- mean(c(r.x1.y * rel.x1, r.x2.y * rel.x2)) / 2
  } else if (length(r.x1x2.y) != 1) {
    stop("r.x1x2.y must be a single value or NULL.")
  }

  #Reliability of x1x2
  rel.x1x2 = ((rel.x1*rel.x2) + (r.x1.x2^2))/(1 + (r.x1.x2^2))

  #If the user inputs a single.N for testing, we will conduct the tests here and skip the rest of the logic by exiting early.
  if(!is.null(single.N)){
    single_condition <- function(n, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width = 0.5, pos.percent = 0.8) {

      `%dopar%` <- foreach::`%dopar%`
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl))
      doParallel::registerDoParallel(cl)

      estimates <- foreach::foreach(i = 1:n.datasets, .packages = c("InteractionPoweR")) %dopar% {
        dataset <- generate_interaction(
          N = n,
          r.x1.y = r.x1.y,
          r.x2.y = r.x2.y,
          r.x1x2.y = r.x1x2.y,
          r.x1.x2 = r.x1.x2,
          rel.x1 = rel.x1,
          rel.x2 = rel.x2,
          rel.y = rel.y
        )
        cor.mat <- stats::cor(dataset)
        base::solve(cor.mat[-3, -3], cor.mat[3, -3])[3]
      }

      power_output <- InteractionPoweR::power_interaction_r2(
        N = n,
        r.x1.y = r.x1.y,
        r.x2.y = r.x2.y,
        r.x1x2.y = r.x1x2.y,
        r.x1.x2 = r.x1.x2,
        rel.x1 = rel.x1,
        rel.x2 = rel.x2,
        rel.y = rel.y,
        detailed_results = TRUE
      )
      expected_value <- power_output$obs.r.x1x2.y
      power <- power_output$pwr

      #This solves for the POS percentage given a COS width
      percent_within_fixed <- base::mean(estimates > expected_value - abs(expected_value) * cos.width &
                                           estimates <= expected_value + abs(expected_value) * cos.width)

      #Sign errors
      sign_error_rate <-  mean(sign(unlist(estimates)) != sign(expected_value))

      #This solves for the COS width given a POS percentage. Robust to negatives, positives, zero.
      f <- function(x) {
        mean(estimates > expected_value - abs(expected_value) * x &
               estimates <= expected_value + abs(expected_value) * x) - pos.percent
      }

      #There is a ceiling of 100x the COS width.
      if(expected_value < 0){
        lower <- 0
        upper <- 1
        max_upper <- 100

      }else{
        lower <- 0
        upper <- 1
        max_upper <- 100

      }
      success <- FALSE

      while (upper <= max_upper) {
        f_lower <- f(lower)
        f_upper <- f(upper)

        if (f_lower * f_upper < 0) {
          sol <- stats::uniroot(f, lower = lower, upper = upper)
          target_cos.width <- sol$root
          success <- TRUE
          break
        } else {
          upper <- upper * 2  # expand search window
        }
      }

      if (!success) {
        warning(paste0("COS width required for ", pos.percent*100, "% of the estimates to be in the COS bounds is greater than 100x the estimate." ))
        target_cos.width <- NA
      }

      percent_within_variable <- base::mean(estimates > expected_value - abs(expected_value) * target_cos.width &
                                              estimates <= expected_value + abs(expected_value) * target_cos.width)

      #This locates the necessary effect size to be stable at the single.N with all other conditions as-is using a search.
      lower_bound <- 0.001
      upper_bound <- 1
      tolerance <- 0.001
      stable_effect_size <- NA
      stable_power <- NA
      max_iter <- 20
      iter <- 0

      while ((upper_bound - lower_bound) > tolerance && iter < max_iter) {
        iter <- iter + 1
        mid_effect <- (lower_bound + upper_bound) / 2

        est_trial <- foreach::foreach(i = 1:n.datasets, .packages = c("InteractionPoweR")) %dopar% {
          dataset <- generate_interaction(
            N = n,
            r.x1.y = r.x1.y,
            r.x2.y = r.x2.y,
            r.x1x2.y = mid_effect,
            r.x1.x2 = r.x1.x2,
            rel.x1 = rel.x1,
            rel.x2 = rel.x2,
            rel.y = rel.y
          )
          cor.mat <- stats::cor(dataset)
          base::solve(cor.mat[-3, -3], cor.mat[3, -3])[3]
        }

        est_trial <- unlist(est_trial)
        power_trial <- InteractionPoweR::power_interaction_r2(
          N = n,
          r.x1.y = r.x1.y,
          r.x2.y = r.x2.y,
          r.x1x2.y = mid_effect,
          r.x1.x2 = r.x1.x2,
          rel.x1 = rel.x1,
          rel.x2 = rel.x2,
          rel.y = rel.y,
          detailed_results = TRUE
        )
        expected_trial <- power_trial$obs.r.x1x2.y
        prop_within <- mean(est_trial > expected_trial - abs(expected_trial) * cos.width &
                              est_trial <= expected_trial + abs(expected_trial) * cos.width)

        if (prop_within >= pos.percent) {
          stable_effect_size <- mid_effect
          stable_power <- power_trial$pwr
          upper_bound <- mid_effect
        } else {
          lower_bound <- mid_effect
        }
      }

      if (is.na(stable_effect_size)) {
        warning("Could not identify stable effect size within specified bounds.")
      }



      results_df <- data.frame(
        single.N = single.N,
        rel.x1 = rel.x1,
        rel.x2 = rel.x2,
        rel.x1x2 = rel.x1x2,
        rel.y = rel.y,
        r.x1x2.y = r.x1x2.y,
        r.x1.y = r.x1.y,
        r.x2.y = r.x2.y,
        r.x1.x2 = r.x1.x2,
        obs.r.x1x2.y = expected_value,
        COS_pos.percents = base::paste0("(", cos.width * 100, ", ", pos.percent * 100, ")"),
        COS_interval = base::paste0("[", base::round(expected_value - abs(expected_value) * cos.width, 4), ", ", base::round(expected_value + abs(expected_value) * cos.width, 4), "]"),
        within_COS_interval = percent_within_fixed,
        POS_determined_COS = base::paste0("[", base::round(expected_value - abs(expected_value) * target_cos.width, 4), ", ", base::round(expected_value + abs(expected_value) * target_cos.width, 4), "]"),        within_POS_determined_COS = percent_within_variable,
        POS_determined_COS.percent = base::round(target_cos.width, 4) * 100,
        sign_error_rate = sign_error_rate,
        power = power,
        r.x1x2.y_to_be_stable_at_single.N = stable_effect_size,
        power_to_be_stable = stable_power,
        stringsAsFactors = FALSE
      )


      return(results_df)
    }


    return(single_condition(single.N, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent))
  }


  #if (!is.null(upper.bound) | !is.null(lower.bound) ) {
  #  stop("Both lower.bound and upper.bound must be of the same type (numeric or NULL) ")
  #}

  if (is.null(lower.bound) & is.null(start.power)) {
    stop("Please select a target power or lower.bound ")
  }

  #This rounding is optional but makes interpretation cleaner.
  round_to_nearest_n <- function(x) {
    base::round(x / step) * step
  }

  #This stores the originally inputted lower.bound for an edge case; see below.
  if (!is.null(lower.bound)) {
    original_lower <- lower.bound
  }else{
    original_lower <- NULL
  }

  #This stores the step value to determine how steps are made later (dynamically or fixed)
  if(is.null(step)){
    dynamicStep <- TRUE
  } else {
    dynamicStep <- FALSE
    fixedStep <- step
  }


  #If a target power value is entered, it will override user selections lower.bound and upper.bound.
  if(!is.null(start.power)){
  #Calculate the necessary sample size for X% power
    find_sample_size_for_power <- function(start.power, r.x1.y, r.x2.y, r.x1x2.y, r.x1.x2,
                                           initial_N = 5, max_N = 1e7) {
      #find an upper bound by doubling the sample size until power meets or exceeds start.power
      current_N <- initial_N
      result <- InteractionPoweR::power_interaction_r2(N = current_N,
                                                       r.x1.y = r.x1.y,
                                                       r.x2.y = r.x2.y,
                                                       r.x1x2.y = r.x1x2.y,
                                                       r.x1.x2 = r.x1.x2,
                                                       rel.x1 = rel.x1,
                                                       rel.x2 = rel.x2)

      while (result$pwr < start.power && current_N < max_N) {
        current_N <- current_N * 2
        result <- InteractionPoweR::power_interaction_r2(N = current_N,
                                                         r.x1.y = r.x1.y,
                                                         r.x2.y = r.x2.y,
                                                         r.x1x2.y = r.x1x2.y,
                                                         r.x1.x2 = r.x1.x2,
                                                         rel.x1 = rel.x1,
                                                         rel.x2 = rel.x2)
      }
      if (current_N >= max_N)
        stop("Could not achieve target power; max sample size exceeded.")
      power.ss <- current_N

      #here is a floor to ensure we ID the true minimum.
      power.floor <- round(power.ss/2)

      #binary search for the s.s. for power of X%.
      while (power.floor < power.ss) {
        mid <- floor((power.floor + power.ss) / 2)
        result_mid <- InteractionPoweR::power_interaction_r2(N = mid,
                                                             r.x1.y = r.x1.y,
                                                             r.x2.y = r.x2.y,
                                                             r.x1x2.y = r.x1x2.y,
                                                             r.x1.x2 = r.x1.x2,
                                                             rel.x1 = rel.x1,
                                                             rel.x2 = rel.x2)
        if (result_mid$pwr >= start.power) {
          power.ss <- mid
        } else {
          power.floor <- mid + 1
        }
      }
      return(power.floor)
    }



    lower.bound <- find_sample_size_for_power(start.power, r.x1.y, r.x2.y, r.x1x2.y, r.x1.x2,
                                           initial_N = 5, max_N = 1e7)
    upper.bound = NULL

    power.at.lb <- InteractionPoweR::power_interaction_r2(N = lower.bound,
                                                          r.x1.y = r.x1.y,
                                                          r.x2.y = r.x2.y,
                                                          r.x1x2.y = r.x1x2.y,
                                                          r.x1.x2 = r.x1.x2,
                                                          rel.x1 = rel.x1,
                                                          rel.x2 = rel.x2)$pwr
  }






  #Replicate the simulation n.datasets times using InteractionPoweR's generate_interaction function. We also calculate power.
  check_condition <- function(n, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent) {
    estimates <- foreach::foreach(i = 1:n.datasets, .packages = c("InteractionPoweR")) %dopar% {
      dataset <- generate_interaction(N = n,
                                      r.x1.y = r.x1.y,
                                      r.x2.y = r.x2.y,
                                      r.x1x2.y = r.x1x2.y,
                                      r.x1.x2 = r.x1.x2,
                                      rel.x1 = rel.x1,
                                      rel.x2 = rel.x2,
                                      rel.y = rel.y)
      cor.mat <- stats::cor(dataset)
      base::solve(cor.mat[-3, -3], cor.mat[3, -3])[3]
    }

    power_output <- InteractionPoweR::power_interaction_r2(
      N = n,
      r.x1.y = r.x1.y,
      r.x2.y = r.x2.y,
      r.x1x2.y = r.x1x2.y,
      r.x1.x2 = r.x1.x2,
      rel.x1 = rel.x1,
      rel.x2 = rel.x2,
      rel.y = rel.y,
      detailed_results = TRUE
    )
    expected_value <- power_output$obs.r.x1x2.y
    percent_within <- base::mean(estimates > expected_value - abs(expected_value) * cos.width&
                                   estimates <= expected_value + abs(expected_value) * cos.width)
    return(list(percent_within = percent_within, power_output = power_output))
  }

  #Searching for the POS with median splits and increments. Caps at 10 million.
  binary_search_sample_size <- function(rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent, lower.bound, upper.bound) {
    if(is.null(upper.bound)) {
      initial_result <- check_condition(lower.bound, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent)
      if (initial_result$percent_within >= pos.percent) {
        current <- lower.bound
        repeat {
          new_n <- max(round(current / 2), 1)
          new_result <- check_condition(new_n, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent)
          if (new_result$percent_within < pos.percent) break
          current <- new_n
        }
        upper.bound <- lower.bound #save the original lower.bound as the upper limit
        lower.bound <- current    #the new lower bound where condition fails
      }
      current_n <- lower.bound
      repeat {
        result <- check_condition(current_n, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent)
        if(result$percent_within >= pos.percent) {
          upper.bound <- current_n
          break
        } else {
          current_n <- current_n * 2
          if(current_n > 1e7) stop("Failed to find an upper bound by 1e7; consider revising your parameters")
        }
      }
    }

    while (lower.bound < upper.bound) {
      if(dynamicStep){
        currentStep <- max(round((upper.bound - lower.bound)/100), 1)
        n <- round((lower.bound + upper.bound) / 2)
      } else {
        currentStep <- fixedStep
        n <- round_to_nearest_n((lower.bound + upper.bound) / 2)
      }
      result <- check_condition(n, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent)
      if (result$percent_within >= pos.percent) {
        upper.bound <- n
      } else {
        lower.bound <- n + currentStep}
    }

    #this calculates power at the updated lower.bound; if it is unchanged, it trips the logic to shrink lower.bound and search further.
    power.at.lb_rolling <- InteractionPoweR::power_interaction_r2(N = lower.bound,
                                                                  r.x1.y = r.x1.y,
                                                                  r.x2.y = r.x2.y,
                                                                  r.x1x2.y = r.x1x2.y,
                                                                  r.x1.x2 = r.x1.x2,
                                                                  rel.x1 = rel.x1,
                                                                  rel.x2 = rel.x2)$pwr

    if(power.at.lb == power.at.lb_rolling) {
      message(paste0('The sample size at the start.power of ', base::round(power.at.lb,2), ' is ',lower.bound ,'. Re-running the search with lower.bound set to ', round(lower.bound/2), ' and upper bound set to ', lower.bound, '. '))
      return(binary_search_sample_size(rel.x1,
                                       rel.x2,
                                       rel.y,
                                       r.x1x2.y,
                                       r.x1.y,
                                       r.x2.y,
                                       r.x1.x2,
                                       cos.width,
                                       pos.percent,
                                       lower.bound = round(lower.bound/2),
                                       upper.bound = lower.bound))
    }


    #this will restart the search if the user specifies a lower.bound and the POS is found to be equal to it.
    if(!is.null(original_lower)){
      if(lower.bound == original_lower && original_lower > 5 && is.null(start.power)) {
        message(base::paste0('The calculated POS equals the input lower.bound. Re-running the search with lower.bound set to ', round(lower.bound/2), ' and upper bound set to ', lower.bound, '. '))
        return(binary_search_sample_size(rel.x1,
                                         rel.x2,
                                         rel.y,
                                         r.x1x2.y,
                                         r.x1.y,
                                         r.x2.y,
                                         r.x1.x2,
                                         cos.width,
                                         pos.percent,
                                         lower.bound = lower.bound/2,
                                         upper.bound = lower.bound))
      }
    }

    #Consolidating results, calculating the COS bounds and f2
    final_result <- check_condition(lower.bound, rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent)
    expected_value <- final_result$power_output$obs.r.x1x2.y
    cos_val <- c(expected_value - abs(expected_value) * cos.width, expected_value + abs(expected_value) * cos.width)

    totalr2 <- base::sum(c(final_result$power_output$b1, final_result$power_output$b2, final_result$power_output$b3) *
                           c(final_result$power_output$obs.r.x1.y, final_result$power_output$obs.r.x2.y, final_result$power_output$obs.r.x1x2.y))
    nullr2 <- base::sum(c(final_result$power_output$b1, final_result$power_output$b2) *
                          c(final_result$power_output$obs.r.x1.y, final_result$power_output$obs.r.x2.y))
    f2 <- (totalr2 - nullr2) / (1 - nullr2)
    return(list(sample_size = lower.bound,
                cos = cos_val,
                f2 = f2,
                power = final_result$power_output$pwr,
                expected_val = expected_value))
  }
    #parallelizing
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(n.cores)
   # on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)

    search_results <- binary_search_sample_size(rel.x1, rel.x2, rel.y, r.x1x2.y, r.x1.y, r.x2.y, r.x1.x2, cos.width, pos.percent, lower.bound, upper.bound)

    #Placing everything into a df
    result_df <- base::data.frame(
      rel.x1 = rel.x1,
      rel.x2 = rel.x2,
      rel.x1x2 = rel.x1x2,
      rel.y = rel.y,
      r.x1x2.y = r.x1x2.y,
      r.x1.y = r.x1.y,
      r.x2.y = r.x2.y,
      r.x1.x2 = r.x1.x2,
      obs.r.x1x2.y = search_results$expected_val,
      COS_pos.percents = base::paste0("(", cos.width * 100, ", ", pos.percent * 100, ")"),
      COS_interval = base::paste0("[", base::round(search_results$cos[1], 4), ", ", base::round(search_results$cos[2], 4), "]"),
      POS_samplesize = search_results$sample_size,
      POS_power = search_results$power,
      stringsAsFactors = FALSE
    )

   # if(cl){
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    #}

    return(result_df)
}
