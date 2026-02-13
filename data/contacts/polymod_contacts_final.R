# 12 June 2025

# Generate contact rate matrix (symmetrical) and fit contact number distributions
# disaggregated by age and setting (home, work/school, travel, other)

#### Components in this script #################################################

# (1) Load POLYMOD Survey data
# (2) Functions
#   - fit_contact_dist() fits statistical distributions to contact count data
#   - assess_model_fits_df() collects fit summary information across multiple
#     contact datasets and enables comparison of different statistical distribution
#     fits
# (3) Build contact dataset starting with list of participants rather than list
#     of contacts. The difference is that this will include participants with zero
#     contacts in all or some settings
# (4) Fit statistical distributions and record results
# (5) Output
# (6) Check that fitted stat distributions are similar to sampled POLYMOD data
# (7) Check can draw a number of contacts from each distribution
# (8) Attempt at fitting zero-inflated negative binomial
# (9) Output contact matrices arranged in an array

################################################################################

#### (1) Load POLYMOD Survey data ##############################################

# socialmixr package
# https://cran.r-project.org/web//packages//socialmixr/vignettes/socialmixr.html
install.packages("socialmixr")
library(socialmixr)
View(socialmixr::list_surveys())
# Want to use a pre-Covid survey for UK
list_surveys()[c(1,7)]
#Key: <date_added>
#  date_added                       title             creator                                    url
#<char>                      <char>              <char>                                 <char>
#  1: 2017-11-07 POLYMOD social contact data        Joël Mossong https://doi.org/10.5281/zenodo.3874557
#  2: 2018-09-05  Social contact data for UK Albert Jan van Hoek https://doi.org/10.5281/zenodo.3874717

# 2 Albert Jan van Hoek survey only for infants so not relevant

#### Load Polymod survey
data("polymod")
get_citation(polymod)
# Mossong J, Hens N, Jit M, Beutels P, Auranen K, Mikolajczyk R, Massari M, Salmaso S, Tomba GS, Wallinga
# J, Heijne J, Sadkowska-Todys M, Rosinska M, Edmunds WJ (2017). “POLYMOD social contact data.”
# doi:10.5281/zenodo.1157934 <https://doi.org/10.5281/zenodo.1157934>, Version 1.1.
saveRDS(polymod, "polymod.rds")
polymod <- readRDS("~/mSCAPE/1_epidemic_modelling/contacts/polymod/polymod.rds")
### Filter for UK participants only
polymod_uk_part <- subset( polymod$participants, polymod$participants$country == "United Kingdom")
polymod_uk_cont <- subset( polymod$contacts, polymod$contacts$part_id %in% polymod_uk_part$part_id )

#### (2) Functions #############################################################

library(fitdistrplus)

# Fitting contact counts to a statistical distributions: negative binomial, exponential, gamma

fit_contact_dist <- function( contact_data, distribution = "nbinom" ){

  # Check if selected distribution is included in this function
  '%ni%' <- Negate('%in%')
  if( distribution %ni% c("nbinom","exp","gamma") ){
    print("Unable to fit. This function only fits to either negative binomial (nbinom), exponential (exp) or gamma distributions")
  } else {

    # Fit distribution
    #cont_dist_fit <- fitdist( contact_data, distribution )
    cont_dist_fit <- tryCatch({
      fitdist(contact_data, distribution)
    }, error = function(e) {
      message("Error in fitdist: ", e$message)
      return(NULL)
    })

    # Display a summary of the fitted distribution
    if( !is.null( cont_dist_fit) ){
      print( summary( cont_dist_fit ) )
      # Define x-axis range and include 0, which may be absent from the data
      x_range <- 0:max( contact_data )

      # Different method for different distributions
      if (distribution == "nbinom"){
        # Extract r and mu negative binomial parameters from fit
        r <- unname( cont_dist_fit$estimate[1] )
        mu <- unname( cont_dist_fit$estimate[2] )
        # Define p in terms of r and mu
        p <- r / (r + mu)
        # Plot histogram of data with fit as line
        hist( contact_data, prob = TRUE)#, breaks = 100)
        lines(x_range, dnbinom(x_range, size = r, prob = p), col = "red", lwd = 2)

      } else if ( distribution == "exp"){
        # Extract r and mu negative binomial parameters from fit
        lambda <- unname( cont_dist_fit$estimate[1] )
        # Plot histogram of data with fit as line
        hist( contact_data, prob = TRUE)#, breaks = 100)
        lines(x_range, dexp(x_range, rate = lambda), col = "blue", lwd = 2)
      } else if ( distribution == "gamma" ){
        # Extract shape and rate parameter estimates for gamma distribution fit
        shape <- unname( cont_dist_fit$estimate[1] )
        rate <- unname( cont_dist_fit$estimate[2] )
        # Plot histogram of data with fit as line
        hist( contact_data, prob = TRUE)#, breaks = 100)
        lines(x_range, dgamma(x_range, shape = shape, scale = 1 / rate), col = "green", lwd = 2)
      }
      # Pause to view hist/line plot
      #readline(prompt = "Press [enter] to continue")

      # Plot diagnostic plots for the fitted distribution
      plot( cont_dist_fit )

      return( cont_dist_fit )

    } else {
      return( NA )
    }
  }
}

### Select statistical model to maximise LL and minimise AIC and BIC

# Cycle through the settings and record the goodness of fit measures
# for the three statistical models (negative binomial, exponential and gamma)
assess_model_fits_df <- function( contact_setting_list, age_group ){

  cont_selection_model_fit <- data.frame( contact_setting = rep(NA, length(contact_setting_list) )
                                          , n_participants = rep(NA, length(contact_setting_list) )
                                          , age_group = rep(NA, length(contact_setting_list) )
                                          , nbinom_r = rep(NA, length(contact_setting_list) )
                                          , nbinom_p = rep(NA, length(contact_setting_list) )
                                          , exp_rate = rep(NA, length(contact_setting_list) )
                                          , gamma_shape = rep(NA, length(contact_setting_list) )
                                          , gamma_scale = rep(NA, length(contact_setting_list) )
                                          , loglik_nbinom = rep(NA, length(contact_setting_list) )
                                          , loglik_exp = rep(NA, length(contact_setting_list) )
                                          , loglik_gamma = rep(NA, length(contact_setting_list) )
                                          , AIC_nbinom = rep(NA, length(contact_setting_list) )
                                          , AIC_exp = rep(NA, length(contact_setting_list) )
                                          , AIC_gamma = rep(NA, length(contact_setting_list) )
                                          , BIC_nbinom = rep(NA, length(contact_setting_list) )
                                          , BIC_exp = rep(NA, length(contact_setting_list) )
                                          , BIC_gamma = rep(NA, length(contact_setting_list) )
                                          , best_loglik = rep(NA, length(contact_setting_list) )
                                          , best_AIC = rep(NA, length(contact_setting_list) )
                                          , best_BIC = rep(NA, length(contact_setting_list) )
                                          , auto_recommend_model = rep(NA, length(contact_setting_list) )
  )

  for( i in 1:length(contact_setting_list)){
    print(paste0("k=",k,", i=",i))
    # Fit distribution models to contact data
    fit_nb <- fit_contact_dist(contact_data = unlist(contact_setting_list[ i ]), distribution = "nbinom")
    r  <- ifelse( unique( is.na( fit_nb ) ), NA, unname( fit_nb$estimate[1] ) )
    mu <- ifelse( unique( is.na( fit_nb ) ), NA, unname( fit_nb$estimate[2] ) )
    p  <- ifelse( unique( is.na( fit_nb ) ), NA, r / (r + mu) )  # Define p in terms of r and mu
    fit_exp <- fit_contact_dist(contact_data = unlist(contact_setting_list[ i ]), distribution = "exp")
    fit_gamma <- fit_contact_dist(contact_data = unlist(contact_setting_list[ i ]), distribution = "gamma")
    # Fill df with model fit parameters
    cont_selection_model_fit[i,"contact_setting"] <- names(contact_setting_list)[ i ]
    cont_selection_model_fit[i,"age_group"] <- age_group #"all"
    cont_selection_model_fit[i,"n_participants"] <- length(unlist(contact_setting_list[1]))
    cont_selection_model_fit[i,"nbinom_r"] <- r
    cont_selection_model_fit[i,"nbinom_p"] <- p
    cont_selection_model_fit[i,"exp_rate"]      <- ifelse( unique( is.na( fit_exp ))  , NA, unname( fit_exp$estimate ) )
    cont_selection_model_fit[i,"gamma_shape"]   <- ifelse( unique( is.na( fit_gamma )), NA, unname( fit_gamma$estimate[1] ) )
    cont_selection_model_fit[i,"gamma_scale"]   <- ifelse( unique( is.na( fit_gamma )), NA, 1 / unname( fit_gamma$estimate[2] ) )
    cont_selection_model_fit[i,"loglik_nbinom"] <- ifelse( unique( is.na( fit_exp ))  , NA, fit_nb$loglik )
    cont_selection_model_fit[i,"loglik_exp"]    <- ifelse( unique( is.na( fit_exp ))  , NA, fit_exp$loglik )
    cont_selection_model_fit[i,"loglik_gamma"]  <- ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$loglik )
    cont_selection_model_fit[i,"AIC_nbinom"]    <- ifelse( unique( is.na( fit_nb ))   , NA, fit_nb$aic )
    cont_selection_model_fit[i,"AIC_exp"]       <- ifelse( unique( is.na( fit_exp ))  , NA, fit_exp$aic )
    cont_selection_model_fit[i,"AIC_gamma"]     <- ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$aic )
    cont_selection_model_fit[i,"BIC_nbinom"]    <- ifelse( unique( is.na( fit_nb ))   , NA, fit_nb$bic )
    cont_selection_model_fit[i,"BIC_exp"]       <- ifelse( unique( is.na( fit_exp ))  , NA, fit_exp$bic )
    cont_selection_model_fit[i,"BIC_gamma"]     <- ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$bic )
    # Goodness of fit parameters for comparison
    logliks <- c( ifelse( unique( is.na( fit_nb )), NA, fit_nb$loglik )
                  , ifelse( unique( is.na( fit_exp )), NA, fit_exp$loglik )
                  , ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$loglik )
    )
    best_loglik <- c("nbinom","exp","gamma")[which( logliks == max( logliks ))]
    #aics <- c( fit_nb$aic, fit_exp$aic, fit_gamma$aic )
    aics <- c( ifelse( unique( is.na( fit_nb )), NA, fit_nb$aic )
               , ifelse( unique( is.na( fit_exp )), NA, fit_exp$aic )
               , ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$aic )
    )
    best_aic <- c("nbinom","exp","gamma")[which( aics == min( aics ))]
    #bics <- c( fit_nb$bic, fit_exp$bic, fit_gamma$bic )
    bics <- c( ifelse( unique( is.na( fit_nb )), NA, fit_nb$bic )
               , ifelse( unique( is.na( fit_exp )), NA, fit_exp$bic )
               , ifelse( unique( is.na( fit_gamma )), NA, fit_gamma$bic )
    )
    best_bic <- c("nbinom","exp","gamma")[which( bics == min( bics ))]
    # Continue filling df
    cont_selection_model_fit[i,"best_loglik"] <- ifelse( length( best_loglik ) == 0, NA, best_loglik )
    cont_selection_model_fit[i,"best_AIC"] <- ifelse( length( best_aic ) == 0, NA, best_aic )
    cont_selection_model_fit[i,"best_BIC"] <- ifelse( length( best_bic ) == 0, NA, best_bic )
    #cont_selection_model_fit[i,"auto_recommend_model"] <- tryCatch({ names( which.max( table( c(best_loglik, best_aic, best_bic) ) ) )[ max( table(c(best_loglik, best_aic, best_bic) ) ) >= 2 ] }
    #                                                               , error = function(e) { message("Error in auto_recommend_model: ", e$message) return(NULL) })
  }
  return( cont_selection_model_fit )
}

#### (3) Build contact dataset starting with list of participants rather than list
#     of contacts. #############################################################

### Build empirical contact distributions from individual participants

# List of participants
part_ids <- unique( polymod_uk_part$part_id )

# Initialise df to store number of contacts for each participant

part_cont <- data.frame( part_id        = rep(NA, length(part_ids))
                       , part_age       = rep(NA, length(part_ids))
                       , cnt_home       = rep(NA, length(part_ids))
                       , cnt_work       = rep(NA, length(part_ids))
                       , cnt_school     = rep(NA, length(part_ids))
                       , cnt_transport  = rep(NA, length(part_ids))
                       , cnt_leisure    = rep(NA, length(part_ids))
                       , cnt_otherplace = rep(NA, length(part_ids))
)

# Loop through participants and aggregate their contacts by setting
#** NOTE**: Some contacts have multiple settings - these are counted in multiple settings
#           and some have no setting recorded - these are therefore not included
for(i in 1:length( part_ids ) ){

  # Filter contact df for individual survey participant
  part_filtered <- subset( polymod_uk_part, polymod_uk_part$part_id == part_ids[ i ] )
  cont_filtered <- subset( polymod_uk_cont, polymod_uk_cont$part_id == part_ids[ i ] )

  # Fill df
  part_cont[i,"part_id"] <- part_ids[ i ]
  part_cont[i,"part_age"] <- part_filtered$part_age
  part_cont[i,"cnt_home"] <- sum( cont_filtered$cnt_home )
  part_cont[i,"cnt_work"] <- sum( cont_filtered$cnt_work )
  part_cont[i,"cnt_school"] <- sum( cont_filtered$cnt_school )
  part_cont[i,"cnt_transport"] <- sum( cont_filtered$cnt_transport )
  part_cont[i,"cnt_leisure"] <- sum( cont_filtered$cnt_leisure )
  part_cont[i,"cnt_otherplace"] <- sum( cont_filtered$cnt_otherplace )

}

View( part_cont )

# Define a list of contact settings, which combines some
contact_setting_list_via_part <- list( home = part_cont$cnt_home
                                       , work_school = part_cont$cnt_work + part_cont$cnt_school
                                       , transport = part_cont$cnt_transport
                                       , leisure = part_cont$cnt_leisure
                                       , other = part_cont$cnt_otherplace
)

# Define age groups
age_groups <- list( all = seq(0,100), "0_3" = seq(0,3), "4_10" = seq(4,10), "11_15" = seq(11,15)
                    , "16_17" = seq(16,17), "18_20" = seq(18,20), "21_24" = seq(21,24)
                    , "25_29" = seq(25,29), "30_34" = seq(30,34), "35_39" = seq(35,39)
                    , "40_44" = seq(40,44), "45_49" = seq(45,49), "50_54" = seq(50,54)
                    , "55_59" = seq(55,59), "60_64" = seq(60,64), "65_69" = seq(65,69)
                    , "70_74" = seq(70,74), "75plus" = seq(75,100))

# Initialise list to hold dataframes containing participant-contact details
# disaggragated by age group
part_cont_by_age_list <- list()

# Loop through age groups including all ages
for (j in 1:length( age_groups )){

  # List of participant IDs
  part_ids <- unique( subset( polymod_uk_part, polymod_uk_part$part_age %in% unname( unlist( age_groups[ j ] ) ) )$part_id )

  # Initialise blank df for participant contacts
  part_cont_temp <- data.frame(   part_id        = rep(NA, length(part_ids))
                                , part_age       = rep(NA, length(part_ids))
                                , cnt_home       = rep(NA, length(part_ids))
                                , cnt_work       = rep(NA, length(part_ids))
                                , cnt_school     = rep(NA, length(part_ids))
                                , cnt_transport  = rep(NA, length(part_ids))
                                , cnt_leisure    = rep(NA, length(part_ids))
                                , cnt_otherplace = rep(NA, length(part_ids))
  )

  # Loop through participants, gathering info on their contacts

  for(i in 1:length( part_ids ) ){

    # Filter contact df for individual survey participant
    part_filtered <- subset( polymod_uk_part, ( polymod_uk_part$part_id == part_ids[ i ] ) )
    cont_filtered <- subset( polymod_uk_cont, ( polymod_uk_cont$part_id == part_ids[ i ] ) )

    # Fill df
    part_cont_temp[i,"part_id"]        <- part_ids[ i ]
    part_cont_temp[i,"part_age"]       <- part_filtered$part_age
    part_cont_temp[i,"cnt_home"]       <- sum( cont_filtered$cnt_home )
    part_cont_temp[i,"cnt_work"]       <- sum( cont_filtered$cnt_work )
    part_cont_temp[i,"cnt_school"]     <- sum( cont_filtered$cnt_school )
    part_cont_temp[i,"cnt_transport"]  <- sum( cont_filtered$cnt_transport )
    part_cont_temp[i,"cnt_leisure"]    <- sum( cont_filtered$cnt_leisure )
    part_cont_temp[i,"cnt_otherplace"] <- sum( cont_filtered$cnt_otherplace )
  }

  # Add to participant-contacts split by setting and age group to the list and name with age group
  part_cont_by_age_list[[ j ]] <- part_cont_temp
  names( part_cont_by_age_list )[ j ] <- names( age_groups )[ j ]
}
rm(i,j)

# Check that this is working as expected
a <- sum( part_cont_by_age_list$all$cnt_home == part_cont[,"cnt_home"] )
b <- length( part_cont[,"cnt_home"] )
a == b #[1] TRUE
rm(a,b)


#### (4) Fit statistical distributions and record results ######################

# Fit statistical distributions to contact data for each age group and measure
# goodness of fit so can decide which stat dist to use in epidemic modelling
#max_cont <- 0
for (k in 1:length( part_cont_by_age_list ) ) {

  print(paste0("k=",k))

  # Age group
  ag <- names( part_cont_by_age_list )[ k ]

  # Contact numbers
  temp_cnt_home <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_home"] ) )
  temp_cnt_work <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_work"] ) )
  temp_cnt_school <- unlist( unname( part_cont_by_age_list[[ag]]["cnt_school"] ) )
  temp_cnt_transport <- unlist( unname( part_cont_by_age_list[[ag]]["cnt_transport"] ) )
  temp_cnt_leisure <- unlist( unname( part_cont_by_age_list[[ag]]["cnt_leisure"] ) )
  temp_cnt_otherplace <- unlist( unname( part_cont_by_age_list[[ag]]["cnt_otherplace"] ) )

  cont_by_setting <- list( home = temp_cnt_home, work = temp_cnt_work, school = temp_cnt_school
                           , transport = temp_cnt_transport, leisure = temp_cnt_leisure
                           , otherplace = temp_cnt_otherplace
                           # Below are composites
                           , work_school = temp_cnt_work + temp_cnt_school
                           , other_t_l_o = temp_cnt_transport + temp_cnt_leisure + temp_cnt_otherplace
                           , other_l_o = temp_cnt_leisure + temp_cnt_otherplace
  )
  #max_cont <- max( max_cont, max(unlist(cont_by_setting) )

  # Fit three statistical models (negative binomial, exponential and gamma) to
  # contact numbers and record fit parameters and goodness of fit measures in df
  temp_fit_assessment_df <- assess_model_fits_df( contact_setting_list = cont_by_setting
                                                , age_group = ag )

  # Add the model fit assessment for individual age group data to a main df to
  # hold all fit data together
  if(k == 1){ fit_assessment_df <- temp_fit_assessment_df }
  else{ fit_assessment_df <- rbind( fit_assessment_df, temp_fit_assessment_df) }

}

saveRDS( fit_assessment_df, "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/fit_assessment_df.rds")
fit_assessment_df <- readRDS( "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/fit_assessment_df.rds" )

# Plot goodness of fit measurements for the diferent statistical distributions
# loglik
plot( fit_assessment_df$loglik_nbinom, typ="l")
lines( fit_assessment_df$loglik_exp, col = "blue" )
lines( fit_assessment_df$loglik_gamma, col = "red" )

# BIC
plot( fit_assessment_df$BIC_nbinom, typ="l")
lines( fit_assessment_df$BIC_exp, col = "blue" )
lines( fit_assessment_df$BIC_gamma, col = "red" )

# AIC
plot( fit_assessment_df$BIC_nbinom, typ="l")
lines( fit_assessment_df$BIC_exp, col = "blue" )
lines( fit_assessment_df$BIC_gamma, col = "red" )


## Plot stat dist fits by setting, by age and by stat dist

library(glue)
# Define unique combinations for plotting
#names(age_groups)
#stat_dists <- c("nbinom","exp","gamma")
contact_settings <- unique( fit_assessment_df$contact_setting )

# Set axis limits
x_limits <- c(0,50)
x_range <- 0:max(x_limits)
y_limits <- c(0,1)#c(0,500)

### Plot by age groups

# Loop through age groups
for(ag in names( part_cont_by_age_list )){

  # Filter data, by age group, to be used in individual plot
  temp_fit_assessment_df <- subset( fit_assessment_df, fit_assessment_df$age_group == ag )

  # Compile contact numbers by contact setting or combinations thereof
  temp_cnt_home <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_home"] ) )
  temp_cnt_work <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_work"] ) )
  temp_cnt_school <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_school"] ) )
  temp_cnt_transport <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_transport"] ) )
  temp_cnt_leisure <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_leisure"] ) )
  temp_cnt_otherplace <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_otherplace"] ) )

  cont_by_setting <- list( home = temp_cnt_home
                         , work = temp_cnt_work
                         , school = temp_cnt_school
                         , transport = temp_cnt_transport
                         , leisure = temp_cnt_leisure
                         , otherplace = temp_cnt_otherplace
                         # Below are composites
                         , work_school = temp_cnt_work + temp_cnt_school
                         , other_t_l_o = temp_cnt_transport + temp_cnt_leisure + temp_cnt_otherplace
                         , other_l_o = temp_cnt_leisure + temp_cnt_otherplace
  )

  # Plots will be set out with all 9 contact types/settings on a single plot view,
  # with a separate plot view for different age groups.
  # Individual plots contain all three statistical distribution fits

  png( glue("~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/statistical_dist_fit_plots/combined_by_age/stat_dist_fits_age_group_{ag}.png")
       , width = 900, height = 600)

  # Set up the plotting grid
  par(mfrow = c(2, 5), mar = c(2, 2, 2, 0.5), oma = c(4, 4, 2, 6)) # oma reserves space for outer margin (legend)[1][3][7]

  # Loop through the contact settings and plot all on a single view
  for( s in 1: ( length(cont_by_setting) + 1) ){

    plot_idx <- s

    # There are 9 plots to be made, but the common legend for all will be placed
    # in the empty 10th position in the grid
    if (plot_idx < 10){

      ylab <- ifelse( plot_idx %in% c(1,6) , "Density" , "" )
      xlab <- ifelse( plot_idx %in% c(6,7,8,9) , "Contacts" , "" )

      hist( cont_by_setting[[ plot_idx ]]
          , prob = TRUE
          , xlim = x_limits, ylim = y_limits
          , xlab = xlab, ylab = ylab
          , main = paste( paste0("age=", ag), names( cont_by_setting )[ plot_idx ], sep = " | ")
          , axes = TRUE
      )
      lines(x_range, dnbinom(x_range, size = temp_fit_assessment_df$nbinom_r[ plot_idx ], prob = temp_fit_assessment_df$nbinom_p[ plot_idx ]), col = "red", lwd = 2)
      lines(x_range, dexp(x_range, rate = temp_fit_assessment_df$exp_rate[ plot_idx ]), col = "blue", lwd = 2)
      lines(x_range, dgamma(x_range, shape = temp_fit_assessment_df$gamma_shape[ plot_idx ], scale = temp_fit_assessment_df$gamma_scale[ plot_idx ] ), col = "green", lwd = 2)
      points( x = names( table( cont_by_setting[[ plot_idx ]] ))
              , y = unname(table( cont_by_setting[[ plot_idx ]] )) / sum( unname( table( cont_by_setting[[ plot_idx ]] )))
              ,pch = 16)

      box()
    } else if ( plot_idx == 10){
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Blank")
      legend("topright", legend = c("survey data - value"
                                    ,"survey data - histogram"
                                    ,"negative binomial"
                                    ,"exponential"
                                    ,"gamma")
             , lty=c(NA,NA,1,1,1)
             , pch = c(16,15,NA,NA,NA)
             , col = c("black","grey","red","blue","green")
             , lwd = 3
      )
    }


  }

  dev.off()
}

### For inspecting individual fits
ag <- "all"
setting <- "otherplace"
cnt_data <- as.data.frame( part_cont_by_age_list[ ag ])[,paste0(ag,".","cnt_",setting)]
# Try excluding zero contacts
cnt_data <- subset( cnt_data, cnt_data != 0 )

library(fitdistrplus)

# Fit negative binomial
fit_nb <- fitdist( cnt_data, "nbinom" )
# nbinom parameter estimates
r <- unname( fit_nb$estimate[1] )
mu <- unname( fit_nb$estimate[2] )
# Define p in terms of r and mu
p <- r / (r + mu)

# Fit exponential distribution
fit_exp <- fitdist( cnt_data, "exp" )
# Extract r and mu negative binomial parameters from fit
lambda <- unname( fit_exp$estimate[1] )

# Fit gamma distribution
fit_gamma <- fitdist( cnt_data, "gamma" )
# Extract shape and rate parameter estimates for gamma distribution fit
shape <- unname( fit_gamma$estimate[1] )
rate <- unname( fit_gamma$estimate[2] )

# Plot histogram of data with fit as line
par(mfrow = c(1, 1))
xrange <- c(0,50)
hist( cnt_data, prob = TRUE, ylim = c(0,1), xlim = xrange)#, breaks = 100)
points( x = as.numeric( names( table( cnt_data ) ) )
        ,y = unname( table( cnt_data ) ) / sum( unname( table( cnt_data ) ) )
        , pch = 16
)
# plot nbinom, exp and gamma distributions
lines(x_range, dnbinom(x_range, size = r, prob = p), col = "red", lwd = 2)
lines(x_range, dexp(x_range, rate = lambda), col = "blue", lwd = 2)
lines(x_range, dgamma(x_range, shape = shape, scale = 1 / rate), col = "green", lwd = 2)

plot(fit_nb)
plot(fit_exp)
plot(fit_gamma)

summary( fit_nb )
summary(fit_exp)
summary(fit_gamma)

### Make separate plots for each combination of age and contact setting
# Loop through age groups
for(ag in names( part_cont_by_age_list )){

  # Filter data, by age group, to be used in individual plot
  temp_fit_assessment_df <- subset( fit_assessment_df, fit_assessment_df$age_group == ag )

  # Compile contact numbers by contact setting or combinations thereof
  temp_cnt_home <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_home"] ) )
  temp_cnt_work <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_work"] ) )
  temp_cnt_school <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_school"] ) )
  temp_cnt_transport <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_transport"] ) )
  temp_cnt_leisure <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_leisure"] ) )
  temp_cnt_otherplace <- unlist( unname( part_cont_by_age_list[[ ag ]]["cnt_otherplace"] ) )

  cont_by_setting <- list( home = temp_cnt_home
                           , work = temp_cnt_work
                           , school = temp_cnt_school
                           , transport = temp_cnt_transport
                           , leisure = temp_cnt_leisure
                           , otherplace = temp_cnt_otherplace
                           # Below are composites
                           , work_school = temp_cnt_work + temp_cnt_school
                           , other_t_l_o = temp_cnt_transport + temp_cnt_leisure + temp_cnt_otherplace
                           , other_l_o = temp_cnt_leisure + temp_cnt_otherplace
  )


  # Loop through the contact settings and plot all on a single view
  for( s in 1: length(cont_by_setting) ){

    # Initialise image file
    png( glue("~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/statistical_dist_fit_plots/separate_by_age/age_group_{ag}_setting_{names( cont_by_setting )[ s ]}.png")
         , width = 1200, height = 800)

    # Plot
    hist( cont_by_setting[[ s ]]
        , prob = TRUE
        , xlim = x_limits, ylim = y_limits
        , xlab = "Number of contacts", ylab = "Density"
        , main = paste( paste0("age=", ag), names( cont_by_setting )[ s ], sep = " | ")
        , axes = TRUE
      )
    lines(x_range, dnbinom(x_range, size = temp_fit_assessment_df$nbinom_r[ s ], prob = temp_fit_assessment_df$nbinom_p[ s ]), col = "red", lwd = 2)
    lines(x_range, dexp(x_range, rate = temp_fit_assessment_df$exp_rate[ s ]), col = "blue", lwd = 2)
    lines(x_range, dgamma(x_range, shape = temp_fit_assessment_df$gamma_shape[ s ], scale = temp_fit_assessment_df$gamma_scale[ s ] ), col = "green", lwd = 2)
    points( x = names( table( cont_by_setting[[ s ]] ))
          , y = unname(table( cont_by_setting[[ s ]] )) / sum( unname( table( cont_by_setting[[ s ]] )))
          , pch = 16)
    box()
    legend("topright", legend = c("survey data - value"
                                 ,"survey data - histogram"
                                 ,"negative binomial"
                                 ,"exponential"
                                 ,"gamma")
             , lty=c(NA,NA,1,1,1)
             , pch = c(16,15,NA,NA,NA)
             , col = c("black","grey","red","blue","green")
             , lwd = 3
    )
  dev.off()
  }
}

#### (5) Output age, setting, selected distribution name, distribution parameters

# Initialise output df
cols <- c("single_year_age","age_group", "setting", "nbinom_r","nbinom_p","gamma_shape","gamma_scale")
polymod_contact_age_setting_distribution <- as.data.frame( matrix(data = NA, nrow = 909, ncol = 7 ) )
names( polymod_contact_age_setting_distribution ) <- cols

row_n <- 0

# Loop through contact settings, age groups and single year ages and fill output df
for( i in 1:length( unique( fit_assessment_df$contact_setting ) ) ){

  c_setting <- unique( fit_assessment_df$contact_setting )[ i ]

  for(j in 1:length( unique( fit_assessment_df$age_group ) ) ){

    ag <- unique( fit_assessment_df$age_group )[ j ]

    # Exclude the 'all' age groups as an option
    if( ag  == "all" ){ next }

    # Convert age groups to constituent single year ages

    # Check if age group string contains "plus" instead of an upper bound
    if( grepl("plus$", ag) ){
      lower <- as.numeric(sub("plus$", "", ag))
      age_group_bounds <- c(lower,100)
    } else {
      age_group_bounds <- as.numeric( unlist( strsplit( ag , "_" ) ) )
    }

    age_group_constituents <- seq( age_group_bounds[1], age_group_bounds[2], 1 )

    for( k in age_group_constituents ){

      row_n <- row_n + 1

      polymod_contact_age_setting_distribution$single_year_age[ row_n ] <- k
      polymod_contact_age_setting_distribution$age_group[ row_n ] <- ag
      polymod_contact_age_setting_distribution$setting[ row_n ] <- c_setting

      fit_assessment_df_subset <- subset( fit_assessment_df, age_group == ag & contact_setting == c_setting )

      polymod_contact_age_setting_distribution$nbinom_r[ row_n ] <- fit_assessment_df_subset$nbinom_r
      polymod_contact_age_setting_distribution$nbinom_p[ row_n ] <- fit_assessment_df_subset$nbinom_p
      polymod_contact_age_setting_distribution$gamma_shape[ row_n ] <- fit_assessment_df_subset$gamma_shape
      polymod_contact_age_setting_distribution$gamma_scale[ row_n ] <- fit_assessment_df_subset$gamma_scale

    }
  }
}

# File converted to single year age
saveRDS( polymod_contact_age_setting_distribution , "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_age_setting_distribution.rds")
polymod_contact_age_setting_distribution <- readRDS( "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_age_setting_distribution.rds")

# File in age groups
contact_setting_age_group_distributions <- fit_assessment_df[ , 1:8]
saveRDS( contact_setting_age_group_distributions , "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/contact_setting_age_group_distributions.rds")

################################################################################

#### (6) Check that fitted stat distributions are similar to sampled POLYMOD data ####################


test_age <- "45_49"
test_setting <- "cnt_work"
dat <- part_cont_by_age_list[[ test_age ]]$cnt_work

samples_nbinom <- rnbinom(n = 10000
                   , size = subset( fit_assessment_df
                                    , contact_setting == "work" &
                                      age_group == "45_49")$nbinom_r
                   , prob = subset( fit_assessment_df
                                    , contact_setting == "work" &
                                      age_group == "45_49")$nbinom_p)

samples_gamma <- rgamma(n = 10000
                          , shape = subset( fit_assessment_df
                                           , contact_setting == "work" &
                                             age_group == "45_49")$gamma_shape
                          , scale = subset( fit_assessment_df
                                           , contact_setting == "work" &
                                             age_group == "45_49")$gamma_scale)

# histogram of POLYMOD data
hist(dat, col=rgb(1,0,0,0.25)
     , xlim=c(0, 50)
     , ylim = c(0,1)
     , probability = TRUE
     , breaks = seq(0,50,1) #c(0,5,10,15,20,25,30,35,40,45,50)
     , main='POLYMOD data and samples from fitted distribution'
     , xlab='Number of contacts'
     , ylab = "Density")

# Add data sampled from fitted gamma distribution
hist(samples_nbinom
     , col=rgb(0,1,0,0.25)
     , add=TRUE
    , xlim=c(0, 50)
    , ylim = c(0,1)
    , probability = TRUE
    , breaks = seq(0,10000,1)#c(0,5,10,15,20,25,30,35,40,45,50,1000)
    )

# Add data sampled from fitted gamma distribution
hist(samples_gamma
     , col=rgb(0,1,1,0.25)
     , add=TRUE
     , xlim=c(0, 50)
     , ylim = c(0,1)
     , probability = TRUE
     , breaks = seq(0,10000,1) #c(0,5,10,15,20,25,30,35,40,45,50,1000)
)

points( x = names( table(dat) ), y = unname(table(dat))/sum(unname(table(dat))), pch = 16)

# Add legend
legend('topright'
       , c('POLYMOD data - values','POLYMOD data - histogram', 'Samples from fitted nbinom distribution', 'Samples from fitted gamma distribution')
       , fill=c(NA,rgb(1,0,0,0.25), rgb(0,1,0,0.25), rgb(0,0,1,0.25))
       , col = c("black",NA,NA,NA)
       , pch = c(16,NA,NA,NA)
       )

################################################################################
#### (7) Check can draw a number of contacts from each distribution ############

View( contact_setting_age_group_distributions )

# Initialise df to record results
test_dist_results <- as.data.frame( matrix( data = NA
                                            , nrow = nrow( contact_setting_age_group_distributions )
                                            , ncol = 6
                                              )
                                    )
colnames( test_dist_results ) <- c("contact_setting","n_participants","age_group","contact_sample_nbinom","contact_sample_exp","contact_sample_gamma")

# Loop through contact settings, age groups
# and generate a sample number of contacts from the three fitted distributions
for( i in 1:nrow(contact_setting_age_group_distributions)){
  test_dist_results[ i, 1:3 ] <- contact_setting_age_group_distributions[ i, 1:3 ]
  test_dist_results[i,"contact_sample_nbinom"] <- rnbinom( n = 1
                                                         , size = contact_setting_age_group_distributions[i,"nbinom_r"]
                                                         , prob = contact_setting_age_group_distributions[i,"nbinom_p"])
  test_dist_results[i,"contact_sample_exp"] <- rexp( n = 1
                                                           , rate = contact_setting_age_group_distributions[i,"exp_rate"])
  test_dist_results[i,"contact_sample_gamma"] <- rgamma( n = 1
                                                        , shape = contact_setting_age_group_distributions[i,"gamma_shape"]
                                                        , scale = contact_setting_age_group_distributions[i,"gamma_scale"])

}

# Currently only using home, work/school and other (transport, leisure, otherplace)
View( subset( test_dist_results, contact_setting %in% c("home","work_school","other_t_l_o") ) )
# but may use more later
# The only combinations of setting, age group and distribution that don't generate
# a number of contacts is work_school, 70_74 and 75plus, with gamma distribution.
# The other distributions generate 0 (or close to 0) for these parameters.

################################################################################
#### (8) Attempt at fitting zero-inflated negative binomial ####################

# Try fitting a zero-inflated negative binomial
install.packages("pscl")    # Run only if not already installed
library(pscl)

df <- cnt_data #data.frame(contacts = contacts)

# Fit a zero-inflated negative binomial model with no predictors (intercept only)
zinb_model <- zeroinfl(cnt_data ~ 1 | 1, data = as.data.frame(df), dist = "negbin")

# View model summary
summary(zinb_model)

# Plot model and data
xrange <- c(0,50)
hist( cnt_data, prob = TRUE, ylim = c(0,1), xlim = xrange)#, breaks = 100)
points( x = as.numeric( names( table( cnt_data ) ) )
        ,y = unname( table( cnt_data ) ) / sum( unname( table( cnt_data ) ) )
        , pch = 16
)
# plot zero-inflated negative binomial model
lines(x_range, dnbinom(x_range, size = r, prob = p), col = "red", lwd = 2)
lines(x_range, dexp(x_range, rate = lambda), col = "blue", lwd = 2)
lines(x_range, dgamma(x_range, shape = shape, scale = 1 / rate), col = "green", lwd = 2)

# Extract estimated parameters
# Count model (negative binomial)
mu <- exp(coef(zinb_model)[["count_(Intercept)"]][1])           # mean of NB (intercept-only model)
theta <- zinb_model$theta                           # NB dispersion parameter

# Zero-inflation model (logit)
pi <- plogis(coef(zinb_model)[["zero_(Intercept)"]][1])         # probability of extra zeros

# Define range of counts to plot (e.g., 0 to max observed + 5)
max_count <- max(cnt_data)
counts <- 0:(max_count + 5)

# Compute ZINB probabilities for each count
zinb_pmf <- function(y, mu, theta, pi) {
  # NB probability for y
  nb_prob <- dnbinom(y, size = theta, mu = mu)
  # ZINB probability
  ifelse(y == 0,
         pi + (1 - pi) * nb_prob,
         (1 - pi) * nb_prob)
}

probs <- zinb_pmf(counts, mu, theta, pi)

# Plot the PMF
barplot(probs, names.arg = counts,
        xlab = "Number of contacts",
        ylab = "Probability",
        main = "Zero-Inflated Negative Binomial PMF")



################################################################################


#### (9) Output contact matrices arranged in an array ##########################

vignette("socialmixr", package = "socialmixr")

# Age groups set at pre-school, primary school, secondary school
# , sixth-form/apprenticeship, young-adult/university, then 5yr age groups
# Values are lower limits of age groups
age_limits = c(0,4,11,16,18,21,25,30,35,40,45,50,55,60,65,70,75)#,80,85,90,95,100)

# Matrix without setting filer
contact_matrix_uk <- contact_matrix( survey = polymod
                                     , countries = "United Kingdom"
                                     , age.limits = age_limits
                                     , estimated.contact.age = "sample" #"mean"
                                     , missing.contact.age = "remove"
                                     , symmetric = TRUE
                                     #, filter = list(cnt_work = 0,cnt_home = 0)#cnt_home = 1,cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                     , weigh.dayofweek = TRUE
                                     , weigh.age = TRUE
                                     , return.part.weights = TRUE
                                     , per.capita = TRUE
                                     #, split = TRUE
)
View(contact_matrix_uk$matrix)
#lines(contact_matrix_uk$matrix[,17])
heatmap(contact_matrix_uk$matrix)

# Contact matrix for HOME contact setting
contact_matrix_uk_home <- contact_matrix( survey = polymod
                                        , countries = "United Kingdom"
                                        , age.limits = age_limits
                                        , estimated.contact.age = "sample" #"mean"
                                        , missing.contact.age = "remove"
                                        , symmetric = TRUE
                                        , filter = list(cnt_home = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                        , weigh.dayofweek = TRUE
                                        , weigh.age = TRUE
                                        , return.part.weights = TRUE
                                        , per.capita = TRUE
                                        #, split = TRUE
)
View(contact_matrix_uk_home$matrix)
heatmap(contact_matrix_uk_home$matrix)

mat <- contact_matrix_uk_home$matrix

### Function for plotting heatmaps
plot_heatmap <- function(mat, title_text){
  # Function to extract the start of the interval
  extract_start <- function(label) as.numeric(sub("\\[(.*),.*\\)", "\\1", label))

  # Sort row and column labels by extracted start
  row_order <- order(sapply(rownames(mat), extract_start))
  col_order <- order(sapply(colnames(mat), extract_start))
  mat_sorted <- mat[row_order, col_order]

  # Plot heatmap without clustering
  heatmap(mat_sorted, Rowv=NA, Colv=NA, scale="none",
          labRow=rownames(mat_sorted), labCol=colnames(mat_sorted),
          xlab="age groups", ylab="age groups", main=title_text)
}

plot_heatmap( mat = contact_matrix_uk_home$matrix
              , title_text = "All settings")


# Contact matrix for WORK contact setting
contact_matrix_uk_work <- contact_matrix( survey = polymod
                                        , countries = "United Kingdom"
                                        , age.limits = age_limits
                                        , estimated.contact.age = "sample" #"mean"
                                        , missing.contact.age = "remove"
                                        , symmetric = TRUE
                                        , filter = list(cnt_work = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                        , weigh.dayofweek = TRUE
                                        , weigh.age = TRUE
                                        , return.part.weights = TRUE
                                        , per.capita = TRUE
                                        #, split = TRUE
)
View(contact_matrix_uk_work$matrix)
heatmap(contact_matrix_uk_work$matrix)
plot_heatmap( mat = contact_matrix_uk_work$matrix
              , title_text = "Work contacts")


# Contact matrix for SCHOOL contact setting
contact_matrix_uk_school <- contact_matrix( survey = polymod
                                          , countries = "United Kingdom"
                                          , age.limits = age_limits
                                          , estimated.contact.age = "sample" #"mean"
                                          , missing.contact.age = "remove"
                                          , symmetric = TRUE
                                          , filter = list(cnt_school = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                          , weigh.dayofweek = TRUE
                                          , weigh.age = TRUE
                                          , return.part.weights = TRUE
                                          , per.capita = TRUE
                                          #, split = TRUE
)
View(contact_matrix_uk_school$matrix)
heatmap(contact_matrix_uk_school$matrix)
plot_heatmap( mat = contact_matrix_uk_school$matrix
              , title_text = "school contacts")


# Contact matrix for WORK & SCHOOL contact setting
# Need to ensure no overlap when summing the two together
# - school could be an adult's work but a child unlikely to have work contacts.
# So set filter for work to cnt_school = 1, cnt_work = 0.

# Contact matrix for SCHOOL NOT WORK contact setting
contact_matrix_uk_school_not_work <- contact_matrix( survey = polymod
                                            , countries = "United Kingdom"
                                            , age.limits = age_limits
                                            , estimated.contact.age = "sample" #"mean"
                                            , missing.contact.age = "remove"
                                            , symmetric = TRUE
                                            , filter = list(cnt_school = 1, cnt_work = 0)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                            , weigh.dayofweek = TRUE
                                            , weigh.age = TRUE
                                            , return.part.weights = TRUE
                                            , per.capita = TRUE
                                            #, split = TRUE
)
View( contact_matrix_uk_school_not_work$matrix )
heatmap(contact_matrix_uk_school_not_work$matrix)
sum( contact_matrix_uk_school_not_work$matrix ) #[1] 27.79072 (very close to but not the same as contact_matrix_uk_school$matrix)

# Contact matrix for WORK and NOT SCHOOL contact setting
contact_matrix_uk_work_not_school <- contact_matrix( survey = polymod
                                                     , countries = "United Kingdom"
                                                     , age.limits = age_limits
                                                     , estimated.contact.age = "sample" #"mean"
                                                     , missing.contact.age = "remove"
                                                     , symmetric = TRUE
                                                     , filter = list(cnt_work = 1, cnt_school = 0)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                                     , weigh.dayofweek = TRUE
                                                     , weigh.age = TRUE
                                                     , return.part.weights = TRUE
                                                     , per.capita = TRUE
                                                     #, split = TRUE
)
View(contact_matrix_uk_work_not_school$matrix)
heatmap(contact_matrix_uk_work_not_school$matrix)
sum(contact_matrix_uk_work_not_school$matrix) #[1] 32.01691 (similar to work (=32.15017) but not exactly the same)

# Contact matrix for WORK AND SCHOOL contact setting
contact_matrix_uk_work_and_school <- contact_matrix( survey = polymod
                                                     , countries = "United Kingdom"
                                                     , age.limits = age_limits
                                                     , estimated.contact.age = "sample" #"mean"
                                                     , missing.contact.age = "remove"
                                                     , symmetric = TRUE
                                                     , filter = list(cnt_work = 1, cnt_school = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                                     , weigh.dayofweek = TRUE
                                                     , weigh.age = TRUE
                                                     , return.part.weights = TRUE
                                                     , per.capita = TRUE
                                                     #, split = TRUE
)
View( contact_matrix_uk_work_and_school$matrix )
heatmap( contact_matrix_uk_work_and_school$matrix )
sum( contact_matrix_uk_work_and_school$matrix ) #[1] 0.1300395
plot_heatmap( mat = contact_matrix_uk_work_and_school$matrix
              , title_text = "work + school contacts")

# Sum SCHOOL and WORK contact matrices together
contact_matrix_uk_school_work <- contact_matrix_uk_school$matrix + contact_matrix_uk_work$matrix

View(contact_matrix_uk_school_work)
heatmap(contact_matrix_uk_school_work)
sum( contact_matrix_uk_school_work ) #[1] 60.06809
sum( contact_matrix_uk_school$matrix ) #[1] 27.91792
sum( contact_matrix_uk_work$matrix ) #[1] 32.15017
# > 27.91792+32.15017
# [1] 60.06809

# Sum SCHOOL (not work) and WORK (not school) contact matrices together
contact_matrix_uk_school_work_ex <- contact_matrix_uk_school_not_work$matrix + contact_matrix_uk_work_not_school$matrix

View(contact_matrix_uk_school_work_ex)
heatmap(contact_matrix_uk_school_work_ex)
sum( contact_matrix_uk_school_work_ex ) #[1] 59.80763
sum( contact_matrix_uk_school_not_work$matrix ) #[1] 27.79072
sum( contact_matrix_uk_work_not_school$matrix ) #[1] 32.01691
# > 27.79072 + 32.01691
# [1] 59.80763

# Sum SCHOOL (not work), WORK (not school), and SCHOOL AND WORK contact matrices together
contact_matrix_uk_school_work_sandw <- contact_matrix_uk_school_not_work$matrix +
                                        contact_matrix_uk_work_not_school$matrix +
                                            contact_matrix_uk_work_and_school$matrix

View( contact_matrix_uk_school_work_sandw )
heatmap( contact_matrix_uk_school_work_sandw )
sum( contact_matrix_uk_school_work_sandw ) #[1] 59.93767
sum( contact_matrix_uk_school_not_work$matrix ) #[1] 27.79072
sum( contact_matrix_uk_work_not_school$matrix ) #[1] 32.01691
sum( contact_matrix_uk_work_and_school$matrix ) #[1] 0.1300395
# > 27.79072 + 32.01691 + 0.1300395
# [1] 59.93767
plot_heatmap( mat = contact_matrix_uk_school_work_sandw
              , title_text = "school + work")


# TRANSPORT
contact_matrix_uk_transport <- contact_matrix( survey = polymod
                                                     , countries = "United Kingdom"
                                                     , age.limits = age_limits
                                                     , estimated.contact.age = "sample" #"mean"
                                                     , missing.contact.age = "remove"
                                                     , symmetric = TRUE
                                                     , filter = list(cnt_transport = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                                     , weigh.dayofweek = TRUE
                                                     , weigh.age = TRUE
                                                     , return.part.weights = TRUE
                                                     , per.capita = TRUE
                                                     #, split = TRUE
)
View(contact_matrix_uk_transport$matrix)
heatmap(contact_matrix_uk_transport$matrix)
sum(contact_matrix_uk_transport$matrix) #[1] 7.483837

# LEISURE
contact_matrix_uk_leisure <- contact_matrix( survey = polymod
                                               , countries = "United Kingdom"
                                               , age.limits = age_limits
                                               , estimated.contact.age = "sample" #"mean"
                                               , missing.contact.age = "remove"
                                               , symmetric = TRUE
                                               , filter = list(cnt_leisure = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                               , weigh.dayofweek = TRUE
                                               , weigh.age = TRUE
                                               , return.part.weights = TRUE
                                               , per.capita = TRUE
                                               #, split = TRUE
)
View( contact_matrix_uk_leisure$matrix )
heatmap( contact_matrix_uk_leisure$matrix )
sum( contact_matrix_uk_leisure$matrix ) #[1] 34.1587

# OTHERPLACE
contact_matrix_uk_otherplace <- contact_matrix( survey = polymod
                                             , countries = "United Kingdom"
                                             , age.limits = age_limits
                                             , estimated.contact.age = "sample" #"mean"
                                             , missing.contact.age = "remove"
                                             , symmetric = TRUE
                                             , filter = list(cnt_otherplace = 1)#cnt_home = 1,cnt_work = 0, cnt_school = 1,cnt_transport = 1,cnt_leisure = 1,cnt_otherplace = 1)
                                             , weigh.dayofweek = TRUE
                                             , weigh.age = TRUE
                                             , return.part.weights = TRUE
                                             , per.capita = TRUE
                                             #, split = TRUE
)
View( contact_matrix_uk_otherplace$matrix )
heatmap( contact_matrix_uk_otherplace$matrix )
sum( contact_matrix_uk_otherplace$matrix ) #[1] 34.41013

# OTHER = sum( TRANSPORT, LEISURE, OTHERPLACE)
# NOTE: May include double counting
contact_matrix_uk_other <- contact_matrix_uk_transport$matrix + contact_matrix_uk_leisure$matrix + contact_matrix_uk_otherplace$matrix
View( contact_matrix_uk_other )
heatmap( contact_matrix_uk_other )
sum( contact_matrix_uk_other ) #[1] 76.05267
sum( contact_matrix_uk_otherplace$matrix ) #[1] 34.41013
sum( contact_matrix_uk_leisure$matrix ) #[1] 34.1587
sum( contact_matrix_uk_transport$matrix ) #[1] 7.483837
34.41013 + 34.1587 + 7.483837 # [1] 76.05267
rownames(contact_matrix_uk_other)
colnames(contact_matrix_uk_other)

#### Save contact matrices to files to be imported into Julia for model

# Need to reformat from matrix to df so col and row names (age groups) are not
# lost when imported into Julia.
# Also adjust age groups so they are exclusive and in same format as contact distribution files

# Reformat age groups
matrix_age_group <- rownames( contact_matrix_uk_home_df )
df_age_group <- c()
for(i in 1:length(matrix_age_group)){
  temp <- unlist( regmatches( matrix_age_group[ i ], gregexpr("[0-9]+", matrix_age_group[ i ] ) ) )
  temp[2] <- ifelse( is.na( temp[2] ) , 101 , temp[2] ) # Applies to highest age category "75plus"
  df_age_group[ i ] <- paste0(temp[1],":",as.character( as.numeric( temp[2] ) - 1 ) )
}

saveRDS( df_age_group , "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_matrix_age_groups.rds" )
# [1] "0:3"    "4:10"   "11:15"  "16:17"  "18:20"  "21:24"  "25:29"  "30:34"  "35:39"  "40:44"  "45:49"
# [12] "50:54"  "55:59"  "60:64"  "65:69"  "70:74"  "75:100"

# HOME
contact_matrix_uk_home_df <- as.data.frame( contact_matrix_uk_home$matrix )
colnames( contact_matrix_uk_home_df ) <- df_age_group
rownames( contact_matrix_uk_home_df ) <- df_age_group
saveRDS( contact_matrix_uk_home_df, "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_matrix_home.rds")
plot_heatmap( mat = contact_matrix_uk_home$matrix, title_text = "home contacts")

# WORK/SCHOOL (avoids double counting by summing 3 contact matrices: w not s, s not w, s and w)
contact_matrix_uk_school_work_sandw_df <- as.data.frame( contact_matrix_uk_school_work_sandw )
colnames( contact_matrix_uk_school_work_sandw_df ) <- df_age_group
rownames( contact_matrix_uk_school_work_sandw_df ) <- df_age_group
heatmap( as.matrix(contact_matrix_uk_school_work_sandw_df) )
saveRDS( contact_matrix_uk_school_work_sandw_df, "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_matrix_school_work.rds")
plot_heatmap( mat = contact_matrix_uk_school_work_sandw, title_text = "work+school contacts")


# Other
contact_matrix_uk_other_df <- as.data.frame( contact_matrix_uk_other )
colnames( contact_matrix_uk_other_df ) <- df_age_group
rownames( contact_matrix_uk_other_df ) <- df_age_group
saveRDS( contact_matrix_uk_other_df, "~/mSCAPE/1_epidemic_modelling/contacts/polymod/2025_06/polymod_contact_matrix_other.rds")
plot_heatmap( mat = contact_matrix_uk_other, title_text = "other contacts")

#** NOTE **
# (1) Some contacts are categorised in multiple settings - so there is some double
#     counting of contacts
# (2) Some contacts have no recorded setting - so they are removed

#** Potential further refinements **
# (A) Teachers in schools - rather than a small probability of adult contact with
#     children across whole population, it should really be close to zero for most
#     and a large number of child contacts for adult teachers
# (B) Adjust work and school contacts by day of week? Otherwise will have larger
#     proportion of zero contacts than should
