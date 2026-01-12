"""
Configuration management for NBPMscape package

Functions:

- load_config(config_path::String):                  Load configuration (model parameters) from a YAML file and return as a nested dictionary.

- validate_config(config::Dict):                     Validate configuration values and return any warnings or errors.

- update_configurable_parameters(P, config::Dict):   Create a new parameter NamedTuple with updated values from config.

- print_changes(P::NamedTuple, config_data, default_configurable_params )   Prints a dataframe to display the default parameters, the parameters in 
                                                                            the configuration file, the parameters after initialization from the 
                                                                            configuration file, and whether each parameter has been updated

- convert_params_to_dfs(P::NamedTuple):              Uses some of the parameter values to populate dataframes which are required within the model


"""

using YAML

"""
Function        load_config(config_path::String)

Description     Load configuration (model parameters) from a YAML file and return as a nested dictionary.

Arguments       config_path::String     Path and filename for .yaml file containing model parameters.

Example         config_data = load_config("config/outbreak_params_covid19_like.yaml")
                config_data = load_config("non_existent_config.yaml")
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

Arguments       config::Dict    Configuration data loaded using load_config function

Returns         Returns vectors of warnings and errors

Example         warnings, errors = validate_config(config_data)
"""
function validate_config(config::Dict) # config=config_data # config=load_config("C:\\Users\\kdrake\\AppData\\Local\\Temp\\jl_7OHsHe7obS.yaml")
    
    # Initialize warnings and errors
    warnings = String[];
    errors = String[];
    
    # Helper function to safely get nested value
    function safe_get_value(config, field_path)
        keys = split(field_path, ".")
        value = config
        try
            for key in keys #key = keys[1]
                if haskey(value, key)
                    value = value[key]
                else
                    return nothing, "Field $field_path not found"
                end
            end
            return value, nothing
        catch e
            return nothing, "Error accessing $field_path: $e"
        end
    end


    # Check that probabilities and proportions are between 0 and 1
    prob_fields = [
        "parameters.prop_severe_hosp_short_stay", # field_path = "parameters.prop_severe_hosp_short_stay"
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
        "parameters.rho_hosp", # Note that ho is used in the yaml file instead of ρ which is used in P
        "parameters.rho_asymptomatic", # Note that ho is used in the yaml file instead of ρ which is used in P
        "parameters.swab_proportion_at_48h"
    ];
    
    for field in prob_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(0 <= value <= 1)
            push!(errors, "Probability or proportion field $field must be between 0 and 1, got $value")
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
    ];
    
    for field in rate_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value < 0
            push!(errors, "Rate field $field must be positive, got $value")
        end
    end
  
    # Check that relative contact rates are not negative
    relative_contact_rate_fields = [
                                    "parameters.fcont"
                                    ,"parameters.gcont"
                                    ,"parameters.oocont"
                                    ];
    
    for field in relative_contact_rate_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value < 0
            push!(errors, "Relative contact rate field $field must not be negative, got $value")
        end
    end
    
    # Check that hospital upper time limits are greater than zero
    hosp_time_ul_fields = ["parameters.tdischarge_ed_upper_limit"
                          ,"parameters.tdischarge_hosp_short_stay_upper_limit"
                          ];
    
    for field in hosp_time_ul_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value <= 0
            push!(errors, "Hospital time upper limit field $field must be greater than zero, got $value")
        end
    end

    # Check that number of imports is greater than zero
    nimports_fields = [ "parameters.nimports" ];
    
    for field in nimports_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value <= 0
            push!(errors, "Number of imports field ($field) must be greater than zero, got $value")
        end
    end

    # Check parameter values are not negative
    other_non_neg_fields = [ "parameters.icu_swab_lag_max"
                            , "parameters.icu_ari_admissions"
                            , "parameters.gp_practices_swab"
                            ,"parameters.gp_swabs_mg"
                            ,"parameters.gp_ari_swabs"
                            ,"parameters.n_hosp_samples_per_week"
                            ,"parameters.swab_time_mode" ];
    
    for field in other_non_neg_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value < 0
            push!(errors, "$field must be not be negative, got $value")
        end
    end

    # Check that Evolutionary parameters are greater than zero
    evo_fields = [ "parameters.mu" #"parameters.μ" mu is used in the yaml file instead of μ
                 , "parameters.omega" ]; #, "parameters.ω" ]; omega is used in the yaml file instead of ω

    for field in evo_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value <= 0
            push!(errors, "$field must be greater than zero, got $value")
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
                  ];
    
    for field in other_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && value < 0
            push!(errors, "$field must not be negative, got $value")
        end
    end

    # Check that ICU sample type is a valid option
    icu_sample_type_fields = [ "parameters.icu_sample_type" ];
    
    for field in icu_sample_type_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(value in ["regional","fixed"])
            push!(errors, "$field must be either regional or fixed, got $value")
        end
    end

    # Check that ICU site stages parameter value is a valid option
    icu_site_stage_fields = [ "parameters.icu_site_stage" ];
    
    for field in icu_site_stage_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(value in ["current","engagement","longlist"])
            push!(errors, "$field must be either current, engagement or longlist, got $value")
        end
    end

    # Check that ICU sampling method is a valid option
    sample_icu_cases_version_fields = [ "parameters.sample_icu_cases_version" ];
    
    for field in sample_icu_cases_version_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(value in ["number","proportion"])
            push!(errors, "$field must be either number or proportion, got $value")
        end
    end

    # Check that pathogen type is a valid option
    pathogen_type_fields = ["parameters.pathogen_type"];
    
    for field in pathogen_type_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(value in ["virus","bacteria","fungi"])
            push!(errors, "$field must be either virus, bacteria or fungi, got $value")
        end
    end

    # Check turnaround times
    turnaround_time_fields = ["parameters.turnaroundtime_icu"
                             ,"parameters.turnaroundtime_hariss"
                             ,"parameters.turnaroundtime_rcgp"
                             ];

    for field in turnaround_time_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing
            # Check that it's a two-element vector
            if !isa(value, Vector) || length(value) != 2
                push!(errors, "$field must be a two element vector, got $(value)")
            elseif value[1] >= value[2]
                push!(errors, "The lower limit (1st element in vector) of $field must be less than the upper limit, got lower limit = $(value[1]) and upper limit = $(value[2])")
            end
        end
    end

    # Warning if turnaround times are different in the different care settings
    tt_icu, error_msg = safe_get_value(config,"parameters.turnaroundtime_icu")
    tt_rcgp, error_msg = safe_get_value(config,"parameters.turnaroundtime_rcgp")
    tt_hariss, error_msg = safe_get_value(config,"parameters.turnaroundtime_hariss")

    if !(tt_icu == tt_rcgp == tt_hariss)
        push!(warnings, "Warning: turnaround times are not all the same! turnaroundtime_icu = $(tt_icu), turnaroundtime_rcgp = $(tt_rcgp), and turnaroundtime_hariss = $(tt_hariss)");
    end 
    
    # Check that sample allocation method is a valid option
    sample_allocation_value, error_msg = safe_get_value(config, "parameters.sample_allocation")
    if error_msg !== nothing
        push!(warnings, "Could not validate sample_allocation: $error_msg")
    elseif sample_allocation_value !== nothing && !(sample_allocation_value in ["equal","weighted"])
        push!(errors, "sample_allocation must be either 'equal' or 'weighted', got $sample_allocation_value")
    end
    
    # Check that sample weighting is a valid option
    weight_samples_by_value, error_msg = safe_get_value(config, "parameters.weight_samples_by")
    if error_msg !== nothing
        push!(warnings, "Could not validate weight_samples_by: $error_msg");
    elseif weight_samples_by_value !== nothing && !(weight_samples_by_value in ["ae_mean","catchment_pop"]) 
        push!(errors, "weight samples by must be either ae_mean or catchment_pop, got $weight_samples_by_value");
    end

    # Check that parameter values are either true or false
    boolean_fields = ["parameters.icu_only_sample_before_death"
                     ,"parameters.hariss_only_sample_before_death"];
    
    for field in boolean_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !(value in [true,false])
            push!(errors, "$field must be either true or false, got $value")
        end
    end

    # Check that initial day of week is within 1 to 7 range
    initial_dow_value, error_msg = safe_get_value(config, "parameters.initial_dow")
    if error_msg !== nothing
        push!(warnings, "Could not validate initial_dow: $error_msg");
    elseif initial_dow_value !== nothing && !(initial_dow_value in [1,2,3,4,5,6,7])
        push!(errors, "initial_dow must be between 1 and 7 (Sunday=1,..., Saturday=7), got $initial_dow_value");
    end

    # Check collection days are all within 1 to 7 range
    phl_collection_dow_value, error_msg = safe_get_value(config, "parameters.phl_collection_dow")
    if error_msg !== nothing
        push!(warnings, "Could not validate phl_collection_dow: $error_msg");
    elseif phl_collection_dow_value !== nothing && !( all(x -> 1 <= x <= 7, phl_collection_dow_value) )
        push!(errors, "All elements in phl_collection_dow must be between 1 and 7 (Sunday=1,..., Saturday=7), got $phl_collection_dow_value");
    end
   
    # Check that the proportion of samples that are required to be adults is either "free" or a Float
    sample_proportion_adult_value, error_msg = safe_get_value(config, "parameters.sample_proportion_adult")
    if error_msg !== nothing
        push!(warnings, "Could not validate sample_proportion_adult: $error_msg");
    elseif sample_proportion_adult_value !== nothing && !( sample_proportion_adult_value == "free" || (isa(sample_proportion_adult_value, Float64) && ( 0 <= sample_proportion_adult_value <= 1 )))
        push!(errors, "sample_proportion_adult must be 'free' or of type Float between 0.0 and 1.0, got $sample_proportion_adult_value");
    end

    # Check day of the week relative contact rates are all non-negative
    dowcont_value, error_msg = safe_get_value(config, "parameters.dowcont")
    if error_msg !== nothing
        push!(warnings, "Could not validate dowcont: $error_msg");
    elseif dowcont_value !== nothing && !( all(x -> 0 <= x , dowcont_value) )
        push!(errors, "All elements in dowcont must be non-negative, got $dowcont_value");
    end

    # Check hospital and ICU ARI admissions adult and child proportions are non-negative and sum to 1
    hosp_ari_admissions_fields = ["parameters.hosp_ari_admissions_adult_p"
                                 ,"parameters.hosp_ari_admissions_child_p"];
    icu_ari_admissions_fields = ["parameters.icu_ari_admissions_adult_p"
                                ,"parameters.icu_ari_admissions_child_p"];
    
    ari_admissions_fields = append!(hosp_ari_admissions_fields, icu_ari_admissions_fields);
       
    # First, check values are between 0 and 1
    for field in ari_admissions_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg");
        elseif value !== nothing && !( all(x -> 0 <= x  <=1 , value) )
            push!(errors, "Value for $field must be non-negative, got $value");
        end
    end

    # Second, check sum to one
    # hosp_ari_admissions
    hosp_ari_adult, error_msg_adult = safe_get_value(config, "parameters.hosp_ari_admissions_adult_p")
    hosp_ari_child, error_msg_child = safe_get_value(config, "parameters.hosp_ari_admissions_child_p")
    try
        if ( hosp_ari_adult + hosp_ari_child ) != 1
            push!(errors, "Values for hosp_ari_admissions_adult_p and hosp_ari_admissions_child_p must sum to one, got $(hosp_ari_adult) and $(hosp_ari_child)")
        end
    catch # Errors in the parameter values are already captured above
    end
    
    # icu_ari_admissions
    icu_ari_adult, error_msg_adult = safe_get_value(config,"parameters.icu_ari_admissions_adult_p")
    icu_ari_child, error_msg_child = safe_get_value(config,"parameters.icu_ari_admissions_child_p")
    try
        if ( icu_ari_adult + icu_ari_child ) != 1
            push!(errors, "Values for icu_ari_admissions_adult_p and icu_ari_admissions_child_p must sum to one, got $(icu_ari_adult) and $(icu_ari_child)")
        end
    catch # Errors in the parameter values are already captured above
    end

    # Check values for destination proportions after attendance at the Emergency Department. All values should be between 0 and 1 and they should sum to 1.
    ed_ari_dest_adult_fields = ["parameters.ed_ari_destinations_adult_p_discharged"
                                 ,"parameters.ed_ari_destinations_adult_p_short_stay"
                                 ,"parameters.ed_ari_destinations_adult_p_longer_stay"];
    ed_ari_dest_child_fields = ["parameters.ed_ari_destinations_child_p_discharged"
                                 ,"parameters.ed_ari_destinations_child_p_short_stay"
                                 ,"parameters.ed_ari_destinations_child_p_longer_stay"];
    ed_ari_dest_fields = append!(ed_ari_dest_adult_fields, ed_ari_dest_child_fields)
    
    # First, check values are between 0 and 1.
    for field in ed_ari_dest_fields
        value, error_msg = safe_get_value(config, field)
        if error_msg !== nothing
            push!(warnings, "Could not validate field $field: $error_msg")
        elseif value !== nothing && !( all(x -> 0 <= x  <=1 , value) )
            push!(errors, "Value for $field must be non-negative, got $value")
        end
    end

    # Second, check sums to one
    # adult
    ed_ari_dest_adult_p_discharged, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_adult_p_discharged")
    ed_ari_dest_adult_p_short_stay, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_adult_p_short_stay")
    ed_ari_dest_adult_p_longer_stay, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_adult_p_longer_stay")

    try
        if ( ed_ari_dest_adult_p_discharged + ed_ari_dest_adult_p_short_stay + ed_ari_dest_adult_p_longer_stay ) != 1
            push!(errors, "Values for ed_ari_destinations_adult_p_discharged, ed_ari_destinations_adult_p_short_stay and ed_ari_destinations_adult_p_longer_stay must sum to one, got $(ed_ari_dest_adult_p_discharged),  $(ed_ari_dest_adult_p_short_stay) and $(ed_ari_dest_adult_p_longer_stay)")
        end
    catch # Errors in the parameter values are already captured above
    end
    
    # child
    ed_ari_dest_child_p_discharged, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_child_p_discharged")
    ed_ari_dest_child_p_short_stay, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_child_p_short_stay")
    ed_ari_dest_child_p_longer_stay, error_msg = safe_get_value(config,"parameters.ed_ari_destinations_child_p_longer_stay")

    try
        if ( ed_ari_dest_child_p_discharged + ed_ari_dest_child_p_short_stay + ed_ari_dest_child_p_longer_stay ) != 1
            push!(errors, "Values for ed_ari_destinations_child_p_discharged, ed_ari_destinations_child_p_short_stay and ed_ari_destinations_child_p_longer_stay must sum to one, got $(ed_ari_dest_child_p_discharged),  $(ed_ari_dest_child_p_short_stay) and $(ed_ari_dest_child_p_longer_stay)")
        end
    catch # Errors in the parameter values are already captured above
    end
        
    return warnings, errors
end



"""
Function    update_configurable_parameters(P, config::Dict)

Description     Create a new parameter NamedTuple with updated values from config.

Arguments       P::NamedTuple   Contains default model parameters
                config::Dict    Contains parameters loaded from configuration file
                default_configurable_params::NamedTuple Contains default configurable parameters for comparison against 

Returns     NamedTuple(P_dict)

Example     P = update_configurable_parameters(P, config_data)
            P = update_configurable_parameters(P, config_data, default_configurable_params)
"""
#function update_configurable_parameters(P, config::Dict) # P = NBPMscape.P # config = config_data
function update_configurable_parameters(P, config::Dict, default_configurable_params::NamedTuple)
    # Convert NamedTuple P to dictionary
    P_dict = Dict(pairs(P)); # println(P_dict)
    
    # Get parameters from config
    if haskey(config, "parameters")
        params = config["parameters"]
        
        for (key, value) in params # println(params)
        #    println(key," ",value)
        #end
            # Convert special names
            param_key = if key == "rho_hosp"
                :ρ_hosp
            elseif key == "rho_asymptomatic"
                :ρ_asymptomatic
            elseif key == "mu"
                :μ
            elseif key == "omega"
                :ω
            #elseif key in ["ed_ari_destinations_adult_p_discharged","ed_ari_destinations_adult_p_short_stay","ed_ari_destinations_adult_p_longer_stay" ] 
             #   :ed_ari_destinations_adult
            #elseif key in ["ed_ari_destinations_child_p_discharged","ed_ari_destinations_child_p_short_stay","ed_ari_destinations_child_p_longer_stay"]
            #    :ed_ari_destinations_child
            else
                Symbol(key)
            end
        
            # Print confirmation that parameter has been updated from configuration file...
            if haskey(P_dict, param_key)
                
                # Parameters in P that are dataframes need to filled with values from the .yaml (config)
                #if key == :ed_ari_destinations_adult_p_discharged
                #    #P.ed_ari_destinations_adult[Ldischarged] = params[ed_ari_destinations_adult_discharged] #config["parameters"]
                #    row_number_discharged_adult = findfirst(P_dict[:ed_ari_destinations_adult].destination .== :discharged)
                #    P_dict[:ed_ari_destinations_adult][row_number_discharged_adult,:proportion_of_attendances] = value #params[ed_ari_destinations_adult_discharged] #config["parameters"]
                #    @info "Updated ed_ari_destinations_adult_p_discharged: $value"
                #elseif key == :ed_ari_destinations_adult_p_short_stay
                #    row_number_short_stay_adult = findfirst(P_dict[:ed_ari_destinations_adult].destination .== :short_stay)
                #    P_dict[:ed_ari_destinations_adult][row_number_short_stay_adult,:proportion_of_attendances] = value
                #    @info "Updated ed_ari_destinations_adult_p_short_stay: $value"
                #elseif key == :ed_ari_destinations_adult_p_longer_stay
                #    row_number_longer_stay_adult = findfirst(P_dict[:ed_ari_destinations_adult].destination .== :longer_stay)
                #    P_dict[:ed_ari_destinations_adult][row_number_longer_stay_adult,:proportion_of_attendances] = value
                #    @info "Updated ed_ari_destinations_adult_p_longer_stay: $value"
                #elseif key == :ed_ari_destinations_child_p_discharged
                #    row_number_longer_stay_child = findfirst(P_dict[:ed_ari_destinations_child].destination .== :discharged)
                #    P_dict[:ed_ari_destinations_child][row_number_longer_stay_child, :proportion_of_attendances] = value
                #    @info "Updated ed_ari_destinations_child_p_discharged: $value"
                #elseif key == :ed_ari_destinations_child_p_short_stay
                #    row_number_short_stay_child = findfirst(P_dict[:ed_ari_destinations_child].destination .== :short_stay)
                #    P_dict[:ed_ari_destinations_child][row_number_short_stay_child,:proportion_of_attendances] = value
                #    @info "Updated ed_ari_destinations_child_p_short_stay: $value"
                #elseif key == :ed_ari_destinations_child_p_longer_stay
                #    row_number_longer_stay_child = findfirst(P_dict[:ed_ari_destinations_child].destination .== :longer_stay)
                #    P_dict[:ed_ari_destinations_child][row_number_longer_stay_child,:proportion_of_attendances] = value
                #    @info "Updated ed_ari_destinations_child_p_longer_stay: $value"
                ## Parameters that are dataframes read from csv files must be done here
                #elseif key == :icu_nhs_trust_sampling_sites
                #    P_dict[:icu_nhs_trust_sampling_sites] = CSV.read(value, DataFrame)
                #    @info "Updated $param_key: $value"
                #elseif key == :hariss_nhs_trust_sampling_sites
                #    P_dict[:hariss_nhs_trust_sampling_sites] = CSV.read(value, DataFrame)
                #    @info "Updated $param_key: $value"
                ## everything else
                #else 

                    # Change type of dowcont value - type in default parameters in core.jl is Tuple but only Vector format allowed in .yaml config file
                    if param_key == :dowcont
                        value = tuple(value...)
                    end

                    P_dict[param_key] = value;
                    @info "Updated $param_key: $value"
                #end
            else
                # ... or not
                # unless parameter is part of a set that will be used to fill a dataframe
                #if key in [ed_ari_destinations_adult_p_discharged,ed_ari_destinations_adult_p_short_stay,ed_ari_destinations_adult_p_longer_stay
                #                 ,ed_ari_destinations_child_p_discharged,ed_ari_destinations_child_p_short_stay,ed_ari_destinations_child_p_longer_stay]
                #else
                    #@warn "Unknown parameter: $key"
                    @warn "Unknown parameter: $param_key"
                #end
            end
        end
    end
    
    # Generate warning if config file does not contain any parameters that are expected and report which parameters
    #keys(P_dict);
    #keys(config["parameters"]);

    return NamedTuple(P_dict)
    #return ( P = NamedTuple(P_dict), 
end




"""
Function        print_changes(P::NamedTuple, config_data, default_configurable_params )

Description     Prints a dataframe to display the default parameters, the parameters in 
                the configuration file, the parameters after initialization from the 
                configuration file, and whether each parameter has been updated.
                Handles cases where config_data has more or fewer parameters than expected.


Arguments   P::NamedTuple   Parameters after updating from the configuration file
            config_data     Parameters from the configuration file
            default_configurable_params     Only the parameters that can be updated
                                            from the configuration file with values
                                            prior to initialization from the 
                                            configuration file

Returns     Prints a dataframe to screen showing values and updates to parameters

"""
function print_changes(P::NamedTuple, config_data, default_configurable_params )
    
    # Create table of default and updated parameters with warnings when required
    params_change_df = DataFrame( model_parameter = String[]
                                , default_value = Any[]
                                , config_file_value = Any[]
                                , value_after_initialization = Any[]
                                , updated = Bool[]
                                , status = String[]
                                );

    # Helper function to safely get config value
    function get_config_value(config_data, param_name_config)
        try
            if haskey(config_data, "parameters") && haskey(config_data["parameters"], String(param_name_config))
                return config_data["parameters"][String(param_name_config)]
            else
                return "Not found in config file"
            end
        catch e
            return "Error accessing config: $e"
        end
    end

    # Helper function to convert parameter names for config lookup
    function get_config_param_name(param_name)
        if param_name == :ρ_hosp
            return "rho_hosp"
        elseif param_name == :ρ_asymptomatic
            return "rho_asymptomatic"
        elseif param_name == :μ
            return "mu"
        elseif param_name == :ω
            return "omega"
        else
            return String(param_name)
        end
    end

    # Helper function to safely compare values
    function values_equal(val1, val2, tolerance=1e-10)
        try
            if (val1 isa Real) && (val2 isa Real)
                return abs(val1 - val2) < tolerance
            else
                return val1 == val2
            end
        catch
            return false
        end
    end

    # Process all parameters from default_configurable_params
    if default_configurable_params !== nothing
        param_names = collect(propertynames(default_configurable_params));
    
        for param_name in param_names #p in 1:length(default_configurable_params) # param_name = param_names[1]
            try
                # Get parameter name for config lookup
                param_name_config = get_config_param_name(param_name);

                # Get values
                default_value = getproperty(default_configurable_params, param_name)
                config_file_value = get_config_value(config_data, param_name_config)
                
                # Get value after initialization (safely)
                value_after_initialization = try
                                                getproperty(P, param_name)
                                            catch
                                                "Not found in P"
                                            end
                                       
                # Determine if updated
                # (includes a tolerance of 10 decimal places when comparing scalar values, but no tolerance for tuples or vectors etc)
                updated_boolean = !values_equal(default_value, value_after_initialization);
                
                # Determine status
                status = if config_file_value == "Not found in config file"
                            "Missing from config"
                        elseif value_after_initialization == "Not found in P"
                            "Missing from final params"
                        elseif updated_boolean
                            "Updated"
                        else
                            "Unchanged"
                        end
                
                # Add row to dataframe for printing
                params_change_df_row = [
                                        String(param_name),
                                        default_value,
                                        config_file_value,
                                        value_after_initialization,
                                        updated_boolean,
                                        status
                                        ];
                push!(params_change_df, params_change_df_row);
            
            catch e

                # Handle any errors for individual parameters
                push!(params_change_df, [
                                        String(param_name),
                                        "Error: $e",
                                        "Error: $e",
                                        "Error: $e",
                                        false,
                                        "Error processing"
                                        ]
                    )
            end
        end
    end

    # Check for extra parameters in config_data that aren't in default_configurable_params
    if haskey(config_data, "parameters")
        config_param_names = Set(keys(config_data["parameters"]))
        expected_param_names = Set()
        
        # Build set of expected parameter names (including name conversions)
        if default_configurable_params !== nothing
            for param_name in propertynames(default_configurable_params)
                push!(expected_param_names, get_config_param_name(param_name))
            end
        end
        
        # Find extra parameters
        extra_params = setdiff(config_param_names, expected_param_names)
        
        for extra_param in extra_params
            try
                config_value = config_data["parameters"][extra_param]
                push!(params_change_df, [
                    extra_param,
                    "Not in defaults",
                    config_value,
                    "Unknown",
                    false,
                    "Extra in config"
                ])
            catch e
                push!(params_change_df, [
                    extra_param,
                    "Error: $e",
                    "Error: $e",
                    "Error: $e",
                    false,
                    "Error processing extra"
                ])
            end
        end
    end

    # Sort by status and parameter name for better readability
    sort!(params_change_df, [:status, :model_parameter]);
    
    # Print summary statistics
    println("\n" * "="^80)
    println("PARAMETER INITIALIZATION SUMMARY")
    println("="^80)
    
    if nrow(params_change_df) > 0
        status_counts = combine(groupby(params_change_df, :status), nrow => :count);
        for row in eachrow(status_counts)
            println("$(row.status): $(row.count) parameters")
        end
        
        println("\nDETAILED PARAMETER CHANGES:")
        println("-"^80)
        show(params_change_df, allrows=true, allcols=true)
        println("\n" * "="^80)
        
        # Print warnings for problematic parameters
        error_rows = filter(row -> occursin("Error", row.status), params_change_df);
        if nrow(error_rows) > 0
            println("\n⚠️  WARNINGS - Parameters with errors:")
            for row in eachrow(error_rows)
                println("  • $(row.model_parameter): $(row.status)")
            end
        end
        
        extra_rows = filter(row -> row.status == "Extra in config", params_change_df);
        if nrow(extra_rows) > 0
            println("\n📋 INFO - Extra parameters found in config file:")
            for row in eachrow(extra_rows)
                println("  • $(row.model_parameter): $(row.config_file_value)")
            end
        end
        
        missing_rows = filter(row -> row.status == "Missing from config", params_change_df);
        if nrow(missing_rows) > 0
            println("\n📝 INFO - Parameters using default values (not in config):")
            for row in eachrow(missing_rows)
                println("  • $(row.model_parameter): $(row.default_value)")
            end
        end
        
    else
        println("No parameters to display")
    end
    
    println("="^80)
end


"""
Function        convert_params_to_dfs(P::NamedTuple)

Description     Uses some of the parameter values to populate dataframes which are required within the model

Arguments   P::NamedTuple   Parameter names and values (defaults updated from configuration file)

Returns     P::NamedTuple with the new dataframe added as additional parameters

"""
function convert_params_to_dfs(P::NamedTuple)
    icu_nhs_trust_sampling_sites = CSV.read( P.icu_nhs_trust_sampling_sites_file, DataFrame );
    P = (; P..., icu_nhs_trust_sampling_sites = icu_nhs_trust_sampling_sites); # P.icu_nhs_trust_sampling_sites
    hariss_nhs_trust_sampling_sites = CSV.read( P.hariss_nhs_trust_sampling_sites_file, DataFrame);
    P = (; P..., hariss_nhs_trust_sampling_sites = hariss_nhs_trust_sampling_sites);
    ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                         , proportion_of_attendances = [P.ed_ari_destinations_adult_p_discharged
                                                                      , P.ed_ari_destinations_adult_p_short_stay
                                                                     , P.ed_ari_destinations_adult_p_longer_stay]
                                         );
    P = (; P..., ed_ari_destinations_adult = ed_ari_destinations_adult);
    ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                         , proportion_of_attendances = [P.ed_ari_destinations_child_p_discharged
                                                                       ,P.ed_ari_destinations_child_p_short_stay
                                                                       ,P.ed_ari_destinations_child_p_longer_stay]
                                          );
    P = (; P..., ed_ari_destinations_child = ed_ari_destinations_child);
    # And remove the keys that have been converted to dataframes	
    #	drop = Set([:icu_nhs_trust_sampling_sites_file, :hariss_nhs_trust_sampling_sites_file
    #				,:ed_ari_destinations_adult_p_discharged, ed_ari_destinations_adult_p_short_stay, ed_ari_destinations_adult_p_longer_stay
    #				,:ed_ari_destinations_child_p_discharged,:ed_ari_destinations_child_p_short_stay,:ed_ari_destinations_child_p_longer_stay
    #				])
    #	P = (; (k => v for (k, v) in pairs(P) if !(k in drop))...)
    return P;
end