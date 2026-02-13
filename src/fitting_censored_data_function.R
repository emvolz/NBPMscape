# Description:
# Adjusting time of detection (TD and 3TD) for censoring due to maximum simulation 
# time (type 1 right-censoring) 
# Function fits distributions to time of infection (tinf) and period between
# tinf and time of detection (treport)
# Author: Kieran Drake
# Date: November 2025
# @param  input_folder    Folder containing the TD data
# @param  input_file      File containing the TD data
# @param  results_type    String describing the simulation/sampling method used to 
#                         generate the TD results. This is used to label output files.
# @param  sample_type     String describing the sampling method, e.g. ICU, used to 
#                         generate the TD results. This is used to label output files.
# @param  td_type         String describing whether the results represent the time to
#                         first detection (TD) or time to 3rd detection (3TD). This 
#                         is used to label output files.
# @param  output_folder   Folder to save output plots and csv files to
# @param  censoring_point Numeric value. Type 1 right censoring point, e.g. for
#                         outbreak simulation this is the maximum time for the
#                         simulation of new infections (maxtime).

# @return Primary outputs from the function are:
#                         - three plots as jpg files:
#                                                     (1) tinf
#                                                     (2) treport-tinf, and
#                                                     (3) treport(TD)
#                         - csv file containing the distribution fit results
#                         - csv file containing samples from the new fitted 
#                           treport distribution
#         No values are expected to be returned, however, the last variable in the
#         function is the distribution fitting results df.
# @examples
# right_censored_data_fit( input_folder = "~/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621_nrep1200_analysis/sample_sensitivity/"
#                        , input_file = "icu_tds_50_samples.csv"
#                        , results_type = "icu_actual_12sites"
#                        , sample_type = "ICU"
#                        , td_type = "TD" #"3TD"
#                        , output_folder = "~/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621_nrep1200_analysis/TD_censored_fit_results/"
#                        , censoring_point = 90
#                        )


library(fitdistrplus)

right_censored_data_fit <- function(  input_folder #= "~/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621_nrep1200_analysis/sample_sensitivity/"
                                    , input_file #= "icu_tds_50_samples.csv"
                                    , results_type #= "icu_actual_12sites"
                                    , sample_type #= "ICU"
                                    #, td_type = "TD" #"3TD"
                                    , output_folder #= "~/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621_nrep1200_analysis/TD_censored_fit_results/"
                                    , censoring_point #= 90
                                    ){

  
  ##### Read in data 
  
  # Read in data from file
  tinf_td_data <- read.csv( paste0(input_folder, input_file) )
  
  simreps <- nrow(tinf_td_data)
  
  ### Run for TD and 3TD
  for (td_type in c("TD","3TD")){
  
    ## Compile censored data in required format
    
    # Different columns for TD and 3TD
    if (td_type == "TD"){
      tinf_col_var = paste0("tinf_relating_to_", sample_type, "_", td_type)
    }
    if (td_type == "3TD"){
      tinf_col_var = paste0("last_tinf_relating_to_", sample_type, "_", td_type)
    }
  
    # Define tinf (time of infection) data
    tinf <- tinf_td_data[ , tinf_col_var] 
    # Define TD (time of detection) data
    td_col_var = paste0(sample_type, "_", td_type)
    td <- tinf_td_data[, td_col_var]
    
    # Create a dataframe for censored data:
    # Left column is NA for censored observations (Type 1 right censoring),
    # Right column is the censoring threshold for censored values,
    # otherwise observed values for uncensored
    censdata <- data.frame(
      left = ifelse(tinf <= censoring_point, tinf, NA),
      right = ifelse(tinf <= censoring_point, tinf, censoring_point )
    )
    
    # Compute median from simulated data (fitted distribution parameters added later)
    tinf_median           <- median(tinf)
    tinf_median_finite    <- median(tinf[is.finite(tinf)])
    
    
    #### Fit distribution to time of infection data
    
    # Create df to store results of fitting
    tinf_fit_results_df = data.frame( data_type = c("tinf","tinf")
                                      ,source = c("tinf_all", "tinf_ex_Inf")
                                      ,median = c(tinf_median,  tinf_median_finite)
                                      ,loglik = c(NA, NA)
                                      ,AIC =    c(NA, NA)
                                      ,BIC =    c(NA, NA)
                                      ,param_1 =c(NA, NA)
                                      ,param_2 =c(NA, NA)
                                      #, best_fit =c("","")
    )
    
    
    # Fit various distributions to censored data
    # Normal distribution fit
    tinf_norm_results_row <- tryCatch({
      # Compute fit to data
      tinf_fit_norm <- fitdistcens(censdata, "norm")
      # Mean
      tinf_fit_norm_median <- as.numeric(tinf_fit_norm$estimate["mean"]) # qnorm(0.5, mean=fit_norm$estimate["mean"], sd = fit_norm$estimate["sd"]) # Normal median equals the mean
      # loglik
      tinf_fit_norm_loglik <- as.numeric(logLik(tinf_fit_norm) )
      # AIC
      tinf_fit_norm_aic <- as.numeric(AIC(tinf_fit_norm))
      #BIC#
      tinf_fit_norm_bic <- as.numeric(BIC(tinf_fit_norm))
      # param 1
      tinf_fit_norm_p1 <- as.numeric(tinf_fit_norm$estimate["mean"][[1]])
      # param 2
      tinf_fit_norm_p2 <- as.numeric(tinf_fit_norm$estimate["sd"][[1]])
      # Add results to a row
      c("tinf","Normal", tinf_fit_norm_median, tinf_fit_norm_loglik, tinf_fit_norm_aic
        , tinf_fit_norm_bic, tinf_fit_norm_p1, tinf_fit_norm_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      tinf_norm_results_row <- c("tinf","Normal", NA, NA, NA, NA, NA, NA)
    })
      
    # Add to tinf results df
    tinf_fit_results_df <- rbind(tinf_fit_results_df, tinf_norm_results_row)
    
    # Gamma distribution fit
    tinf_gamma_results_row <- tryCatch({
      tinf_fit_gamma <- fitdistcens(censdata, "gamma")
      tinf_fit_gamma_median <- qgamma(0.5, shape = tinf_fit_gamma$estimate["shape"]
                                         , rate = tinf_fit_gamma$estimate["rate"])
      tinf_fit_gamma_loglik <- as.numeric(logLik(tinf_fit_gamma) )
      tinf_fit_gamma_aic <- as.numeric(AIC(tinf_fit_gamma)) 
      tinf_fit_gamma_bic <- as.numeric(BIC(tinf_fit_gamma))
      tinf_fit_gamma_p1 <- as.numeric(tinf_fit_gamma$estimate["shape"][[1]])
      tinf_fit_gamma_p2 <- as.numeric(tinf_fit_gamma$estimate["rate"][[1]])
      # Add results to a row
      c("tinf","Gamma", tinf_fit_gamma_median, tinf_fit_gamma_loglik, tinf_fit_gamma_aic
        , tinf_fit_gamma_bic, tinf_fit_gamma_p1, tinf_fit_gamma_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      tinf_gamma_results_row <- c("tinf","Gamma", NA, NA, NA, NA, NA, NA)
    })
    
    # Add to tinf results df
    tinf_fit_results_df <- rbind(tinf_fit_results_df, tinf_gamma_results_row)
    
    #tinf_fit_results_df<-tinf_fit_results_df[1:2,]
    
    # Weibull distribution fit
    tinf_weibull_results_row <- tryCatch({
      # Compute fit to data
      tinf_fit_weibull <- fitdistcens(censdata, "weibull")
      tinf_fit_weibull_median <- qweibull(0.5, shape = tinf_fit_weibull$estimate["shape"]
                                             , scale = tinf_fit_weibull$estimate["scale"])
      tinf_fit_weibull_loglik <- as.numeric(logLik(tinf_fit_weibull) )
      tinf_fit_weibull_aic <- as.numeric(AIC(tinf_fit_weibull))
      tinf_fit_weibull_bic <- as.numeric(BIC(tinf_fit_weibull))
      tinf_fit_weibull_p1 <- as.numeric(tinf_fit_weibull$estimate["shape"][[1]])
      tinf_fit_weibull_p2 <- as.numeric(tinf_fit_weibull$estimate["scale"][[1]])
      # Add results to a row
      c("tinf","Weibull", tinf_fit_weibull_median, tinf_fit_weibull_loglik, tinf_fit_weibull_aic
        , tinf_fit_weibull_bic, tinf_fit_weibull_p1, tinf_fit_weibull_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      tinf_weibull_results_row <- c("tinf","Weibull", NA, NA, NA, NA, NA, NA)
    })
    # Add to tinf results df
    tinf_fit_results_df <- rbind(tinf_fit_results_df, tinf_weibull_results_row)
  
    # Lognormal distribution fit
    tinf_lnorm_results_row <- tryCatch({
      # Compute fit to data
      tinf_fit_lnorm <- fitdistcens(censdata, "lnorm")
      tinf_fit_lnorm_median <- exp(tinf_fit_lnorm$estimate["meanlog"])
      tinf_fit_lnorm_loglik <- as.numeric( logLik(tinf_fit_lnorm) )
      tinf_fit_lnorm_aic <- as.numeric( AIC(tinf_fit_lnorm) )
      tinf_fit_lnorm_bic <- as.numeric( BIC(tinf_fit_lnorm) )
      tinf_fit_lnorm_p1 <- as.numeric( tinf_fit_lnorm$estimate["shape"][[1]] )
      tinf_fit_lnorm_p2 <- as.numeric( tinf_fit_lnorm$estimate["rate"][[1]] )
      # Add results to a row to be added to results df
      c("tinf","Lognormal", tinf_fit_lnorm_median, tinf_fit_lnorm_loglik, tinf_fit_lnorm_aic
        , tinf_fit_lnorm_bic, tinf_fit_lnorm_p1
                                  , tinf_fit_lnorm_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      tinf_lnorm_results_row <- c("tinf","Lognormal", NA, NA, NA, NA, NA, NA)
    })
    # Add to tinf results df
    tinf_fit_results_df <- rbind(tinf_fit_results_df, tinf_lnorm_results_row)
    
    #summary(tinf_fit_norm)
    #summary(tinf_fit_gamma)
    #summary(tinf_fit_weibull)
    #summary(tinf_fit_lnorm)
    
    
    # Determine the best fitting distribution
    tinf_best_fit_ix <- which(as.numeric(tinf_fit_results_df$loglik) == min(tinf_fit_results_df$loglik, na.rm = TRUE))
    tinf_fit_results_df[tinf_best_fit_ix,"best_fit"] = 1
    
    # Write data to file for each sim rep
    #write.csv(  tinf_fit_results_df
    #          , file = paste0(output_folder,results_type,"_",td_type,"_tinf_fit_results_table.csv")
    #          , row.names = FALSE)
    
    
    
    
    ####
    # Fit time period between tinf and treport (TD)
    # (Then add the median of this period to the median fitted tinf to get an estimate
    # for the median treport)
    
    # Define time period
    treport_tinf_period <- td - tinf
    treport_tinf_period_ex_inf <- treport_tinf_period[ !is.nan( treport_tinf_period ) ]
    
    # Compute medians of time period
    treport_tinf_period_median <- median(td - tinf)
    treport_tinf_period_median_finite <- median( td[is.finite(td)] - tinf[is.finite(tinf)] )
    
    # Create df to store results of fitting
    treport_tinf_fit_results_df = data.frame( data_type = c("treport_tinf","treport_tinf")
                                      ,source = c("treport_tinf_all", "treport_tinf_ex_inf")
                                      ,median = c( treport_tinf_period_median
                                                  ,  treport_tinf_period_median_finite)
                                      ,loglik = c(NA, NA)
                                      ,AIC =    c(NA, NA)
                                      ,BIC =    c(NA, NA)
                                      ,param_1 =c(NA, NA)
                                      ,param_2 =c(NA, NA)
                                      #, best_fit =c("","")
    )
    
    ## Fit various distributions to data (NOTE that the treport - tinf period is NOT censored)
    
    # Normal distribution fit
    treport_tinf_norm_results_row <- tryCatch({
      # Compute fit to data
      treport_tinf_fit_norm <- fitdist(treport_tinf_period_ex_inf, "norm")
      # Mean
      treport_tinf_fit_norm_median <- as.numeric(treport_tinf_fit_norm$estimate["mean"]) # qnorm(0.5, mean=fit_norm$estimate["mean"], sd = fit_norm$estimate["sd"]) # Normal median equals the mean
      # loglik
      treport_tinf_fit_norm_loglik <- as.numeric(logLik(treport_tinf_fit_norm) )
      # AIC
      treport_tinf_fit_norm_aic <- as.numeric(AIC(treport_tinf_fit_norm))
      #BIC#
      treport_tinf_fit_norm_bic <- as.numeric(BIC(treport_tinf_fit_norm))
      # param 1
      treport_tinf_fit_norm_p1 <- as.numeric(treport_tinf_fit_norm$estimate["mean"][[1]])
      # param 2
      treport_tinf_fit_norm_p2 <- as.numeric(treport_tinf_fit_norm$estimate["sd"][[1]])
      # Add results to a row
      c("treport_tinf","Normal", treport_tinf_fit_norm_median, treport_tinf_fit_norm_loglik
        , treport_tinf_fit_norm_aic, treport_tinf_fit_norm_bic, treport_tinf_fit_norm_p1
        , treport_tinf_fit_norm_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      treport_tinf_norm_results_row <- c("treport_tinf","Normal", NA, NA, NA, NA, NA, NA)
    })
    
    # Add to treport_tinf results df
    treport_tinf_fit_results_df <- rbind(treport_tinf_fit_results_df, treport_tinf_norm_results_row)
    
    # Gamma distribution fit
    treport_tinf_gamma_results_row <- tryCatch({
      treport_tinf_fit_gamma <- fitdist(treport_tinf_period_ex_inf, "gamma")
      treport_tinf_fit_gamma_median <- qgamma(0.5, shape = treport_tinf_fit_gamma$estimate["shape"]
                                      , rate = treport_tinf_fit_gamma$estimate["rate"])
      treport_tinf_fit_gamma_loglik <- as.numeric(logLik(treport_tinf_fit_gamma) )
      treport_tinf_fit_gamma_aic <- as.numeric(AIC(treport_tinf_fit_gamma)) 
      treport_tinf_fit_gamma_bic <- as.numeric(BIC(treport_tinf_fit_gamma))
      treport_tinf_fit_gamma_p1 <- as.numeric(treport_tinf_fit_gamma$estimate["shape"][[1]])
      treport_tinf_fit_gamma_p2 <- as.numeric(treport_tinf_fit_gamma$estimate["rate"][[1]])
      # Add results to a row
      c("treport_tinf","Gamma", treport_tinf_fit_gamma_median, treport_tinf_fit_gamma_loglik
        , treport_tinf_fit_gamma_aic, treport_tinf_fit_gamma_bic, treport_tinf_fit_gamma_p1
        , treport_tinf_fit_gamma_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      treport_tinf_gamma_results_row <- c("treport_tinf","Gamma", NA, NA, NA, NA, NA, NA)
    })
    
    # Add to tinf results df
    treport_tinf_fit_results_df <- rbind(treport_tinf_fit_results_df, treport_tinf_gamma_results_row)
    
    #tinf_fit_results_df<-tinf_fit_results_df[1:2,]
    
    # Weibull distribution fit
    treport_tinf_weibull_results_row <- tryCatch({
      # Compute fit to data
      treport_tinf_fit_weibull <- fitdist(treport_tinf_period_ex_inf, "weibull")
      treport_tinf_fit_weibull_median <- qweibull(0.5, shape = treport_tinf_fit_weibull$estimate["shape"]
                                          , scale = treport_tinf_fit_weibull$estimate["scale"])
      treport_tinf_fit_weibull_loglik <- as.numeric(logLik(treport_tinf_fit_weibull) )
      treport_tinf_fit_weibull_aic <- as.numeric(AIC(treport_tinf_fit_weibull))
      treport_tinf_fit_weibull_bic <- as.numeric(BIC(treport_tinf_fit_weibull))
      treport_tinf_fit_weibull_p1 <- as.numeric(treport_tinf_fit_weibull$estimate["shape"][[1]])
      treport_tinf_fit_weibull_p2 <- as.numeric(treport_tinf_fit_weibull$estimate["scale"][[1]])
      # Add results to a row
      c("treport_tinf","Weibull", treport_tinf_fit_weibull_median, treport_tinf_fit_weibull_loglik
        , treport_tinf_fit_weibull_aic, treport_tinf_fit_weibull_bic, treport_tinf_fit_weibull_p1
        , treport_tinf_fit_weibull_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      treport_tinf_weibull_results_row <- c("treport_tinf","Weibull", NA, NA, NA, NA, NA, NA)
    })
    # Add to tinf results df
    treport_tinf_fit_results_df <- rbind(treport_tinf_fit_results_df, treport_tinf_weibull_results_row)
    
    # Lognormal distribution fit
    treport_tinf_lnorm_results_row <- tryCatch({
      # Compute fit to data
      treport_tinf_fit_lnorm <- fitdist(treport_tinf_period_ex_inf, "lnorm")
      treport_tinf_fit_lnorm_median <- exp(treport_tinf_fit_lnorm$estimate["meanlog"])
      treport_tinf_fit_lnorm_loglik <- as.numeric( logLik(treport_tinf_fit_lnorm) )
      treport_tinf_fit_lnorm_aic <- as.numeric( AIC(treport_tinf_fit_lnorm) )
      treport_tinf_fit_lnorm_bic <- as.numeric( BIC(treport_tinf_fit_lnorm) )
      treport_tinf_fit_lnorm_p1 <- as.numeric( treport_tinf_fit_lnorm$estimate["meanlog"][[1]] )
      treport_tinf_fit_lnorm_p2 <- as.numeric( treport_tinf_fit_lnorm$estimate["sdlog"][[1]] )
      # Add results to a row to be added to results df
      c("treport_tinf","Lognormal", treport_tinf_fit_lnorm_median, treport_tinf_fit_lnorm_loglik
        , treport_tinf_fit_lnorm_aic, treport_tinf_fit_lnorm_bic, treport_tinf_fit_lnorm_p1
        , treport_tinf_fit_lnorm_p2)
    }, error = function(e) {
      # If the Weibull distribution fit fails then return an empty row
      treport_tinf_lnorm_results_row <- c("treport_tinf","Lognormal", NA, NA, NA, NA, NA, NA)
    })
    # Add to treport_tinf results df
    treport_tinf_fit_results_df <- rbind(treport_tinf_fit_results_df, treport_tinf_lnorm_results_row)
    
    #summary(treport_tinf_fit_norm)
    #summary(treport_tinf_fit_gamma)
    #summary(treport_tinf_fit_weibull)
    #summary(treport_tinf_fit_lnorm)
    
    # Determine the best fitting distribution
    treport_tinf_best_fit_ix <- which(as.numeric(treport_tinf_fit_results_df$loglik) == min(treport_tinf_fit_results_df$loglik, na.rm = TRUE))
    treport_tinf_fit_results_df[treport_tinf_best_fit_ix,"best_fit"] = 1
    
    # Merge results dataframes for tinf and (treport-tinf) and then save file
    final_fit_results_df <- rbind(tinf_fit_results_df, treport_tinf_fit_results_df)
    write.csv( final_fit_results_df
               , file = paste0(output_folder,results_type,"_",td_type,"_fit_results.csv")
               , row.names = FALSE )
    
    ##### Output samples from fitted distribution
    
    ## Estimate TD (and 3TD) from fitted data
    # Draw samples from best fit distribution for tinf and best fit distribution for
    # (treport-tinf). Then add together to get fitted distribution for treport (= TD or 3TD)
    tinf_fit_samples <- switch( tinf_best_fit_ix # Row index of best-fit distribution in results df
                               , NULL # rows 1 and 2 in dfare simulated data and not best-fit distribution
                               , NULL
                               , rnorm( 10000, mean = tinf_fit_norm$estimate["mean"], sd = tinf_fit_norm$estimate["sd"] )
                               , rgamma( 10000, shape = tinf_fit_gamma$estimate["shape"], rate = tinf_fit_gamma$estimate["rate"] )
                               , rweibull( 10000, shape = tinf_fit_weibull$estimate["shape"], scale = tinf_fit_weibull$estimate["scale"] )
                               , rlnorm( 10000, meanlog = tinf_fit_lnorm$estimate["meanlog"], sdlog = tinf_fit_lnorm$estimate["sdlog"] )
    )
     
    treport_tinf_fit_samples <- switch( treport_tinf_best_fit_ix # Row index of best-fit distribution in results df
                                       , NULL # rows 1 and 2 in dfare simulated data and not best-fit distribution
                                       , NULL
                                       , rnorm( 10000, mean = treport_tinf_fit_norm$estimate["mean"], sd = treport_tinf_fit_norm$estimate["sd"] )
                                       , rgamma( 10000, shape = treport_tinf_fit_gamma$estimate["shape"], rate = treport_tinf_fit_gamma$estimate["rate"] )
                                       , rweibull( 10000, shape = treport_tinf_fit_weibull$estimate["shape"], scale = treport_tinf_fit_weibull$estimate["scale"] )
                                       , rlnorm( 10000, meanlog = treport_tinf_fit_lnorm$estimate["meanlog"], sdlog = treport_tinf_fit_lnorm$estimate["sdlog"] )
    ) 
      
    treport_fit_samples <- tinf_fit_samples + treport_tinf_fit_samples
    
    
    # Output fitted TD values
    df <- data.frame( c1 = treport_fit_samples )
    col.name <- paste0(sample_type,"_",td_type)
    names(df) <- col.name
    write.csv( df
               , file = paste0(output_folder,results_type,"_",td_type
                               #,"_"
                               #, tinf_fit_results_df[ tinf_best_fit_ix, "source"],"_"
                               #, treport_tinf_fit_results_df[ treport_tinf_best_fit_ix,"source"]
                               ,"_fit_samples.csv")
               , row.names = FALSE )
  
  
    
    #### Plots
    
    # Before plotting fit results, save the current graphical parameters
    old_par <- par(no.readonly = TRUE)
    
  
    
    # Save to file histograms and fits of:
    # (1) time of infection data
    # (2) Period between time of report (TD) and time of infection
    # (3) time of report (detection)
    
    jpeg(paste0(output_folder,results_type,"_",td_type,"_fits.jpg")
         #, width=900, height=900
         , units = "in", width = 6, height = 9
         , res = 200, quality=100)
    
    # Split plot area into 3 vertical sections (3 rows, 1 column) so can add all 
    # to same output file
    par(mfcol = c(3, 1))
    
    # Save to file a histogram of time of infection data on density scale
    #jpeg(paste0(output_folder,results_type,"_",td_type,"_tinf_fit.jpg")
    #     , width=480, height=480, quality=150)
    
    hist(tinf, breaks = 30, freq = FALSE, col = "lightblue"
         , main = paste0("tinf censored at ",censoring_point," days with best fit distribution")
         , xlab = "Days since first importation", ylab = "Density"
         , xlim = c(0,200)
         , ylim = c(0,0.1)
    )
    
    # Create sequence of x values for fitted line
    #xfit <- seq(min(x[is.finite(x)]), max(x[is.finite(x)]), length.out = 200)
    tinf_xfit <- seq(0, 200, length.out = 200)
    
    # Compute fitted density values using estimated parameters
    # y-values from best-fit distribution
    tinf_yfit <- switch( tinf_best_fit_ix # Row index of best-fit distribution in results df
                              , NULL # rows 1 and 2 in df are simulated data and not best-fit distribution
                              , NULL
                              , dnorm( tinf_xfit, mean = tinf_fit_norm$estimate["mean"], sd = tinf_fit_norm$estimate["sd"] )
                              , dgamma( tinf_xfit, shape = tinf_fit_gamma$estimate["shape"], rate = tinf_fit_gamma$estimate["rate"] )
                              , dweibull( tinf_xfit, shape = tinf_fit_weibull$estimate["shape"], scale = tinf_fit_weibull$estimate["scale"] )
                              , dlnorm( tinf_xfit, meanlog = tinf_fit_lnorm$estimate["meanlog"], sdlog = tinf_fit_lnorm$estimate["sdlog"] )
    )
    
    # Add fitted density line to histogram
    lines(tinf_xfit, tinf_yfit, col = "green", lwd = 2)
    # and fitted median
    abline(v = tinf_fit_results_df[tinf_best_fit_ix,"median"], col = "green",  lty = 2, lwd = 2)
    # Add median from simulated (right-censored) data
    abline(v = median(tinf), col = "blue",  lty = 2, lwd = 2)
    
    # Add legend
    legend("topright"
           ,legend = c("tinf"
                       , "tinf median"
                       , paste0(tinf_fit_results_df[tinf_best_fit_ix,"source"]," fit")
                       , paste0(tinf_fit_results_df[tinf_best_fit_ix,"source"]," fit median")
           )
           ,fill = c("lightblue",NA,NA,NA)
           ,border = c("black",NA,NA,NA)
           ,col = c(NA, "blue","green","green")
           ,lty = c(NA,2,1,2)
           ,bty = "n"
           ,pt.cex = 2
           ,cex = 0.8)
    #dev.off()
    
    
    # Save histogram of treport minus tinf period to file, along with fit and median values
    #jpeg(paste0(output_folder,results_type,"_",td_type,"_treport_tinf_fit.jpg")
    #     , width=480, height=480, quality=150)
    
    hist(treport_tinf_period, breaks = 30, freq = FALSE, col = "lightgreen"
         , main = "treport minus tinf with best fit distribution"
         , xlab = "Period between tinf and treport", ylab = "Density"
         , xlim = c(0,50)
         , ylim = c(0,0.2)
         #, add = TRUE
    )
    
    # Add median from simulated (right-censored) data
    abline(v = median( treport_tinf_period_ex_inf ), col = "darkgreen",  lty = 2, lwd = 2)
    
    # Compute fitted density values using estimated parameters
    treport_tinf_xfit <- seq(0, 50, length.out = 200)
    
    # y-values from best-fit distribution
    treport_tinf_yfit <- switch( treport_tinf_best_fit_ix # Row index of best-fit distribution in results df
                                 , NULL # rows 1 and 2 in dfare simulated data and not best-fit distribution
                                 , NULL
                                 , dnorm( treport_tinf_xfit, mean = treport_tinf_fit_norm$estimate["mean"], sd = treport_tinf_fit_norm$estimate["sd"] )
                                 , dgamma( treport_tinf_xfit, shape = treport_tinf_fit_gamma$estimate["shape"], rate = treport_tinf_fit_gamma$estimate["rate"] )
                                 , dweibull( treport_tinf_xfit, shape = treport_tinf_fit_weibull$estimate["shape"], scale = treport_tinf_fit_weibull$estimate["scale"] )
                                 , dlnorm( treport_tinf_xfit, meanlog = treport_tinf_fit_lnorm$estimate["meanlog"], sdlog = treport_tinf_fit_lnorm$estimate["sdlog"] )
    )
    
    # Add fitted density line to histogram
    lines(treport_tinf_xfit, treport_tinf_yfit, col = "red", lwd = 2)
    # and fitted median
    abline(v = treport_tinf_fit_results_df[ treport_tinf_best_fit_ix, "median" ]
           , col = "red",  lty = 2, lwd = 2)
    
    # Add legend
    legend("topright"
           ,legend = c("Period between treport and tinf"
                       , "Period between treport and tinf median"
                       , paste0(treport_tinf_fit_results_df[treport_tinf_best_fit_ix,"source"]," fit")
                       , paste0(treport_tinf_fit_results_df[treport_tinf_best_fit_ix,"source"]," fit median")
           )
           ,fill = c("lightgreen",NA,NA,NA)
           ,border = c("black",NA,NA,NA)
           ,col = c(NA, "darkgreen","red","red")
           ,lty = c(NA,2,1,2)
           , lwd = 2
           ,bty = "n"
           ,pt.cex = 2
           ,cex = 0.8)
    #dev.off()
    
    # Save plot of distribution fitted to treport (TD)
    #jpeg(paste0(output_folder,results_type,"_",td_type,"_td_fit.jpg")
    #     , width=480, height=480, quality=150)
    
    # Fitted distribution
    hist( treport_fit_samples, breaks = 30, freq = FALSE
          , col = rgb(1, 0, 0, alpha = 0.3)
          , main = paste0("Simulated TDs and fitted TDs adjusted for right censoring at ",censoring_point," days")
          , xlab = "Days since first importation", ylab = "Density"
          , xlim = c(0,200)
          , ylim = c(0,0.1)
          #,alpha=0.3
          #, add = TRUE
    )
    # Simulated right-censored data
    hist(td, breaks = 30, freq = FALSE
         , rgb(0, 0, 1, alpha = 0.3)
         , main = ""
         #, xlab = "Days since first importation", ylab = "Density"
         , xlim = c(0,200)
         , ylim = c(0,0.1)
         , add = TRUE
    )
    # Add medians
    abline(v = median(treport_fit_samples),  col = "red",    lty = 2, lwd = 2)
    abline(v = median(td), col = "blue",  lty = 2, lwd = 2)
    
    legend("topright"
           , legend = c( paste0("Simulation TDs right censored at ",censoring_point," days (",simreps," simreps)")
                         , paste0("Adjusted TDs using "
                                  , tinf_fit_results_df[ tinf_best_fit_ix, "source"]
                                  ," fit for tinf and "
                                  , treport_tinf_fit_results_df[ treport_tinf_best_fit_ix,"source"]
                                  ," fit for (treport-tinf) (10000 samples)")
                         , "Simulated TDs median"
                         , "Adjusted TD fitted median")
           , fill = c(rgb(0, 0,1, alpha = 0.3), rgb(1, 0, 0, alpha = 0.3),NA,NA)
           , border = c("black","black",NA,NA)
           , col = c(NA,NA, "blue","red")
           , lty = c(NA, NA,2,2)
           , bty = "n"
           , pt.cex = 2
           , cex = 0.8 )
    dev.off()
    
    # Return to previous plot parameters
    par(old_par)
  }  
}
