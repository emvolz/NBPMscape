"""
Configuration management for NBPMscape package

Functions:

- load_config(config_path::String):                  Load configuration (model parameters) from a YAML file and return as a nested dictionary.

- validate_config(config::Dict):                     Validate configuration values and return any warnings or errors.

- update_configurable_parameters(P, config::Dict):   Create a new parameter NamedTuple with updated values from config.


"""

using YAML

"""
Function        load_config(config_path::String)

Description     Load configuration (model parameters) from a YAML file and return as a nested dictionary.

Arguments       config_path::String     Path and filename for .yaml file containing model parameters.

Example         load_config("config/outbreak_params_covid19_like.yaml")
"""
function load_config(config_path::String)
    if !isfile(config_path)
        error("Configuration file not found: $config_path")
    end
    
    try
        return YAML.load_file(config_path)
    catch e
        error("Failed to parse YAML configuration file: $e")
    end
end


"""
Function        validate_config(config::Dict)

Description     Validate configuration values and return any warnings or errors.
"""
function validate_config(config::Dict) # config=config_data
    
    # Initialize warnings and errors
    warnings = String[]
    errors = String[]
    
    # Check that probabilities and proportions are between 0 and 1
    prob_fields = [
        "parameters.prop_severe_hosp_short_stay",
        "parameters.prop_severe_hosp_long_stay", 
        "parameters.prop_moderate_ED",
        "parameters.prop_moderate_GP",
        "parameters.prop_mild",
        "parameters.p_sampled_icu",
        "parameters.sample_target_prob_icu",
        "parameters.proportion_hosp_swabbed",
        "parameters.sensitivity_mg_virus",
        "parameters.sensitivity_mg_bacteria",
        "parameters.sensitivity_mg_fungi",
        "parameters.rho_hosp",
        "parameters.rho_asymptomatic",
        "parameters.swab_proportion_at_48h"
    ]
    
    for field in prob_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(0 <= value <= 1)
                push!(errors, "Probability or proportion field $field must be between 0 and 1, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    # Check that rates are greater than zero
    rate_fields = [
        "parameters.gp_only_rate",
        "parameters.ed_direct_rate",
        "parameters.gp_before_hosp_rate",
        "parameters.ed_from_gp_rate",
        "parameters.hosp_admit_direct_rate",
        "parameters.hosp_admit_from_gp_rate",
        "parameters.hosp_recovery_rate",
        "parameters.hosp_short_stay_recovery_rate",
        "parameters.hosp_long_stay_recovery_rate",
        "parameters.hosp_death_rate",
        "parameters.triage_icu_rate",
        "parameters.icu_to_death_rate",
        "parameters.icu_to_stepdown_leading_to_recovery_rate",
        "parameters.stepdown_to_recovery_after_icu_rate",
        "parameters.infectivity",
        "parameters.frate",
        "parameters.grate",
        "parameters.commuterate",
        "parameters.importrate"
    ]
    
    for field in rate_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value < 0 #<= 0
                push!(errors, "Rate field $field must be positive, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    # Check that relative contact rates are not negative
    relative_contact_rate_fields = [
                                    "parameters.fcont"
                                    ,"parameters.gcont"
                                    ,"parameters.oocont"
                                    ]
    
    for field in relative_contact_rate_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value < 0 #<= 0
                push!(errors, "Relative contact rate field $field must not be negative, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that hospital upper time limits
    hosp_time_ul_fields = ["parameters.tdischarge_ed_upper_limit"
                          ,"parameters.tdischarge_hosp_short_stay_upper_limit"
                          ]
    
    for field in hosp_time_ul_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value <= 0 #<= 0
                push!(errors, "Hospital time upper limit field $field must be greater than zero, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that number of imports is greater than zero
    nimports_fields = [ "parameters.nimports" ]
    
    for field in nimports_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value <= 0 #<= 0
                push!(errors, "Number of imports field ($field) must be greater than zero, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check parameter values are not negative
    other_non_neg_fields = [ "parameters.icu_swab_lag_max"
                            , "parameters.icu_ari_admissions"
                            , "parameters.gp_practices_swab"
                            ,"parameters.gp_swabs_mg"
                            ,"parameters.gp_ari_swabs"
                            ,"parameters.n_hosp_samples_per_week"
                            ,"parameters.swab_time_mode" ]
    
    for field in other_non_neg_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value < 0 #<= 0
                push!(errors, "$field must be not be negative, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that number of imports is greater than zero
    evo_fields = [ "parameters.μ"
                 , "parameters.ω" ]
    
    for field in nimports_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value <= 0 #<= 0
                push!(errors, "$field must be greater than zero, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check parameter values are greater than zero
    other_fields = [ "parameters.n_icu_samples_per_week"
                 , "parameters.icu_ari_admissions"
                 , "parameters.gp_practices_total"
                 , "parameters.pop_eng"
                 , "parameters.gp_ari_consults"
                 , "parameters.hariss_courier_to_analysis"
                 , "parameters.hosp_ari_admissions"
                  ]
    
    for field in nimports_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value <= 0 #<= 0
                push!(errors, "$field must be greater than zero, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that number of imports is greater than zero
    icu_sample_type_fields = [ "parameters.icu_sample_type" ]
    
    for field in icu_sample_type_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["regional","fixed"])
                push!(errors, "$field must be either regional or fixed, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that number of imports is greater than zero
    icu_site_stage_fields = [ "parameters.icu_site_stage" ]
    
    for field in icu_site_stage_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["current","engagement","longlist"])
                push!(errors, "$field must be either current, engagement or longlist, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that number of imports is greater than zero
    sample_icu_cases_version_fields = [ "parameters.sample_icu_cases_version" ]
    
    for field in sample_icu_cases_version_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["number","proportion"])
                push!(errors, "$field must be either number or proportion, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that pathogen type is a valid option
    pathogen_type_fields = ["parameters.pathogen_type"]
    
    for field in pathogen_type_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["virus","bacteria","fungi"]) # value = "virus"
                push!(errors, "$field must be either virus, bacteria or fungi, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check turnaround times
    turnaround_time_fields = ["parameters.turnaroundtime_icu"
                             ,"parameters.turnaroundtime_hariss"
                             ,"parameters.turnaroundtime_rcgp"
                             ]
    #value_vec = []
    #warnings=[]
    #errors=[]
    for field in turnaround_time_fields
        keys = split(field, ".") #field=turnaround_time_fields[1]
        value = config
        try
            for key in keys
                value = value[key]
                println(value) #push!(value_vec,value)
            end
            # Check that the lower limit is less than the upper limit
            if length(value) != 2 # value=[2,2]
                push!(errors, "$field must be a two element vector, got $(value)")
            end
            if value[1] >= value[2]
                push!(errors, "The lower limit (1st element in vector) of $field must be less than the upper limit, got lower limit = $(value[1]) and upper limit = $(value[2])")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Warning if turnaround times are different in the different care settings
    tt_icu = config["parameters"]["turnaroundtime_icu"]
    tt_rcgp = config["parameters"]["turnaroundtime_rcgp"]
    tt_hariss = config["parameters"]["turnaroundtime_hariss"]   

    if !(tt_icu == tt_rcgp == tt_hariss)
        push!(warnings, "Warning: turnaround times are not all the same! turnaroundtime_icu = $(tt_icu), turnaroundtime_rcgp = $(tt_rcgp), and turnaroundtime_hariss = $(tt_hariss)")
    end 
    

    # Check that sample allocation method is a valid option
    sample_allocation_fields = ["parameters.sample_allocation"]
    
    for field in sample_allocation_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["equal","weighted"]) 
                push!(errors, "$field must be either equal or weighted, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that sample weighting is a valid option
    weight_samples_by_fields = ["parameters.weight_samples_by"]
    
    for field in weight_samples_by_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in ["ae_mean","catchment_pop"]) 
                push!(errors, "$field must be either ae_mean or catchment_pop, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that parameter values are either true or false
    boolean_fields = ["parameters.icu_only_sample_before_death"
                     ,"parameters.hariss_only_sample_before_death"]
    
    for field in boolean_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in [true,false])
                push!(errors, "$field must be either true or false, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check that initial day of week is within 1 to 7 range
    initial_dow_fields = ["parameters.initial_dow"]
    
    for field in initial_dow_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(value in [1,2,3,4,5,6,7])
                push!(errors, "$field must be between 1 and 7 (Sunday=1,..., Saturday=7), got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    # Check collection days are all within 1 to 7 range
    phl_collection_dow_fields = ["parameters.phl_collection_dow"]
    
    for field in phl_collection_dow_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !( all(x -> 1 <= x <= 7, value) )
                push!(errors, "All elements in $field must be between 1 and 7 (Sunday=1,..., Saturday=7), got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    # Check that the proportion of samples that are required to be adults is either "free" or a Float
    sample_proportion_adult_fields = ["parameters.sample_proportion_adult"]
    
    for field in sample_proportion_adult_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !( value == "free" || isa(value, Float64) )
                push!(errors, "$field must be free or of type Float between 0.0 and 1.0, got $value")
            end
            if isa(value, Float64)
                if !( 0 <= value <= 1 )
                    push!(errors, "$field must be between 0.0 and 1.0, got $value")
                end
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # Check day of the week relative contact rates are all non-negative
    dowcont_fields = ["parameters.dowcont"]
    
    for field in dowcont_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !( all(x -> 0 <= x , value) )
                push!(errors, "All elements in $field must be non-negative, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end


    # Check hospital and ICU ARI admissions adult and child proportions are non-negative and sum to 1
    hosp_ari_admissions_fields = ["parameters.hosp_ari_admissions_adult_p"
                                 ,"parameters.hosp_ari_admissions_child_p"]
    icu_ari_admissions_fields = ["parameters.icu_ari_admissions_adult_p"
                                ,"parameters.icu_ari_admissions_child_p"]
    
    # First, check sums to one
    # hosp_ari_admissions
    try
        hosp_ari_adult = config["parameters"]["hosp_ari_admissions_adult_p"]
        hosp_ari_child = config["parameters"]["hosp_ari_admissions_child_p"]

        if ( hosp_ari_adult + hosp_ari_child ) != 1
            push!(errors, "Values for hosp_ari_admissions_adult_p and hosp_ari_admissions_child_p must sum to one, got $(hosp_ari_adult) and $(hosp_ari_child)")
        end
    catch
        push!(warnings, "Could not validate fields hosp_ari_admissions_adult_p and hosp_ari_admissions_child_p")
    end

    # icu_ari_admissions
    try
        icu_ari_adult = config["parameters"]["icu_ari_admissions_adult_p"]
        icu_ari_child = config["parameters"]["icu_ari_admissions_child_p"]

        if ( icu_ari_adult + icu_ari_child ) != 1
            push!(errors, "Values for icu_ari_admissions_adult_p and icu_ari_admissions_child_p must sum to one, got $(icu_ari_adult) and $(icu_ari_child)")
        end
    catch
        push!(warnings, "Could not validate fields icu_ari_admissions_adult_p and icu_ari_admissions_child_p")
    end

    # Second, check non-negative
    # hosp_ari_admissions
    for field in hosp_ari_admissions_fields
        keys = split(field, ".")
        value = config
        
        try
            for key in keys
                value = value[key]
                #println(value)
            end
            if !( all(x -> 0 <= x , value) )
                push!(errors, "Value for $field must be non-negative, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    # icu_ari_admissions
    for field in icu_ari_admissions_fields
        keys = split(field, ".")
        value = config
        
        try
            for key in keys
                value = value[key]
            end
            if !( all(x -> 0 <= x , value) )
                push!(errors, "Value for $field must be non-negative, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end

    return warnings, errors
end



"""
Function    update_configurable_parameters(P, config::Dict)

Description     Create a new parameter NamedTuple with updated values from config.

Arguments       P::NamedTuple   Contains default model parameters
                config::Dict    Contains parameters loaded from configuration file

"""
function update_configurable_parameters(P, config::Dict)
    # Convert NamedTuple P to dictionary
    P_dict = Dict(pairs(P))
    
    # Get parameters from config
    if haskey(config, "parameters")
        params = config["parameters"]
        
        for (key, value) in params
            # Convert special names
            param_key = if key == "rho_hosp"
                :ρ_hosp
            elseif key == "rho_asymptomatic"
                :ρ_asymptomatic
            elseif key == "mu"
                :μ
            elseif key == "omega"
                :ω
            else
                Symbol(key)
            end
            
            # Print confirmation that parameter has been updated from configuration file...
            if haskey(P_dict, param_key)
                P_dict[param_key] = value
                @info "Updated $param_key: $value"
            else
                # ... or not
                @warn "Unknown parameter: $key"
            end
        end
    end
    
    return NamedTuple(P_dict)
end

## OLD VERSION
#function update_configurable_parameters(P, config::Dict)
#    flat_config = flatten_dict(config)
#    P_dict = Dict(pairs(P))
#    param_mapping = get_parameter_mapping()
#    
#    updated_params = []
#    for (config_key, param_key) in param_mapping
#        if haskey(flat_config, config_key)
#            # Handle tuple conversion for dowcont
#            if param_key == :dowcont && isa(flat_config[config_key], Vector)
#                P_dict[param_key] = tuple(flat_config[config_key]...)
#            else
#                P_dict[param_key] = flat_config[config_key]
#            end
#            push!(updated_params, param_key)
#        end
#    end
#    
#    if !isempty(updated_params)
#        @info "Updated $(length(updated_params)) parameters: $(join(updated_params, ", "))"
#    end
#    
#    return NamedTuple(P_dict)
#end

"""
NOT IN USE
Function        get_config_path(config_name::String="default_config.yaml")

Description     Get the full path to a configuration file in the config directory.
"""
#function get_config_path(config_name::String="default_config.yaml")
#    pkg_dir = dirname(dirname(@__FILE__))
#    config_dir = joinpath(pkg_dir, "config")
#    
#    if !isdir(config_dir)
#        mkpath(config_dir)
#    end
#    
#    return joinpath(config_dir, config_name)
#end

"""
NOT IN USE
Function        flatten_dict(d::Dict, parent_key::String="", sep::String="_")

Description     Flatten a nested dictionary for easier parameter mapping.
"""
#function flatten_dict(d::Dict, parent_key::String="", sep::String="_")
#    items = []
#    for (k, v) in d
#        new_key = parent_key == "" ? string(k) : "$(parent_key)$(sep)$(k)"
#        if isa(v, Dict)
#            append!(items, flatten_dict(v, new_key, sep))
#        else
#            push!(items, (Symbol(new_key), v))
#        end
#    end
#    return Dict(items)
#end


"""
NOT IN USE
Function        get_parameter_mapping()

Description     Get the mapping dictionary from config keys to parameter names.
"""
#function get_parameter_mapping()
#    return Dict(
#        ## Pathogen specific parameters ##
#        # Transmission parameters
#        :transmission_infectivity => :infectivity,
#        :transmission_infectivity_shape => :infectivity_shape,
#        :transmission_infectivity_scale => :infectivity_scale,
#        # Disease progression parameters
#        :disease_progression_latent_shape => :latent_shape,
#        :disease_progression_latent_scale => :latent_scale,
#        :disease_progression_infectious_shape => :infectious_shape,
#        :disease_progression_infectious_scale => :infectious_scale,
#        # Transmission reduction factors
#        :transmission_reduction_rho_hosp => :ρ_hosp,
#        :transmission_reduction_rho_asymptomatic => :ρ_asymptomatic,
#        # Infection severity probabilities (age-independent)
#        :severity_probabilities_prop_severe_hosp_short_stay => :prop_severe_hosp_short_stay,
#        :severity_probabilities_prop_severe_hosp_long_stay => :prop_severe_hosp_long_stay,
#        :severity_probabilities_prop_moderate_ED => :prop_moderate_ED,
#        :severity_probabilities_prop_moderate_GP => :prop_moderate_GP,
#        :severity_probabilities_prop_mild => :prop_mild,
#        # Healthcare progression rates (1/duration in days)
#        :healthcare_rates_gp_only_rate => :gp_only_rate,
#        :healthcare_rates_ed_direct_rate => :ed_direct_rate,
#        :healthcare_rates_gp_before_hosp_rate => :gp_before_hosp_rate,
#        :healthcare_rates_ed_from_gp_rate => :ed_from_gp_rate,
#        :healthcare_rates_hosp_admit_direct_rate => :hosp_admit_direct_rate,
#        :healthcare_rates_hosp_admit_from_gp_rate => :hosp_admit_from_gp_rate,
#        :healthcare_rates_hosp_recovery_rate => :hosp_recovery_rate,
#        :healthcare_rates_hosp_short_stay_recovery_rate => :hosp_short_stay_recovery_rate,
#        :healthcare_rates_hosp_long_stay_recovery_rate => :hosp_long_stay_recovery_rate,
#        :healthcare_rates_hosp_death_rate => :hosp_death_rate,
#       :healthcare_rates_triage_icu_rate => :triage_icu_rate,
#        :healthcare_rates_icu_to_death_rate => :icu_to_death_rate,
#        :healthcare_rates_icu_to_stepdown_leading_to_recovery_rate => :icu_to_stepdown_leading_to_recovery_rate,
#        :healthcare_rates_stepdown_to_recovery_after_icu_rate => :stepdown_to_recovery_after_icu_rate,
#        # Import parameters
#        :imports_importrate => :importrate,
#        :imports_nimports => :nimports,
#        :imports_import_t_df => :import_t_df,
#        :imports_import_t_s => :import_t_s,
#        # Molecular evolution parameters
#        :evolution_mu => :μ,
#        :evolution_omega => :ω,
#        
#        ## Non-pathogen specific parameters ##
#        # Relative contact rates
#        :contact_rates_fcont => :fcont,
#        :contact_rates_gcont => :gcont,
#        :contact_rates_oocont => :oocont,
#        # Scale adjustment for contact rates by day of the week (Sun-Sat)
#       :day_of_week_scaling_dowcont => :dowcont,
#        # Network dynamics
#        :network_dynamics_frate => :frate,
#        :network_dynamics_grate => :grate,
#        # Sampling and testing parameters
#        :sampling_psampled => :psampled,
#        :sampling_sample_target_prob => :sample_target_prob,
#        :sampling_turnaroundtime_icu => :turnaroundtime_icu,
#        :sampling_turnaroundtime_rcgp => :turnaroundtime_rcgp,
#        :sampling_turnaroundtime_hariss => :turnaroundtime_hariss,
#        :sampling_icu_swab_lag_max => :icu_swab_lag_max,
#        # Discharge timing
#        :discharge_timing_tdischarge_ed_upper_limit => :tdischarge_ed_upper_limit,
#        :discharge_timing_tdischarge_hosp_short_stay_upper_limit => :tdischarge_hosp_short_stay_upper_limit,
#        # Commute parameters
#        :commute_commuterate => :commuterate,
#        # Metagenomic testing sensitivity
#        :metagenomic_sensitivity_sensitivity_mg_virus => :sensitivity_mg_virus,
#        :metagenomic_sensitivity_sensitivity_mg_bacteria => :sensitivity_mg_bacteria,
#        :metagenomic_sensitivity_sensitivity_mg_fungi => :sensitivity_mg_fungi
#    )
#end