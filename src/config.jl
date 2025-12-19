"""
Configuration management for NBPMscape package
"""

using YAML

"""
    load_config(config_path::String)

Load configuration from a YAML file and return as a nested dictionary.
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
    get_config_path(config_name::String="default_config.yaml")

Get the full path to a configuration file in the config directory.
"""
function get_config_path(config_name::String="default_config.yaml")
    pkg_dir = dirname(dirname(@__FILE__))
    config_dir = joinpath(pkg_dir, "config")
    
    if !isdir(config_dir)
        mkpath(config_dir)
    end
    
    return joinpath(config_dir, config_name)
end

"""
    flatten_dict(d::Dict, parent_key::String="", sep::String="_")

Flatten a nested dictionary for easier parameter mapping.
"""
function flatten_dict(d::Dict, parent_key::String="", sep::String="_")
    items = []
    for (k, v) in d
        new_key = parent_key == "" ? string(k) : "$(parent_key)$(sep)$(k)"
        if isa(v, Dict)
            append!(items, flatten_dict(v, new_key, sep))
        else
            push!(items, (Symbol(new_key), v))
        end
    end
    return Dict(items)
end

"""
    validate_config(config::Dict)

Validate configuration values and return any warnings or errors.
"""
function validate_config(config::Dict)
    warnings = String[]
    errors = String[]
    
    # Check that probabilities are between 0 and 1
    prob_fields = [
        "severity_probabilities.prop_severe_hosp_short_stay",
        "severity_probabilities.prop_severe_hosp_long_stay", 
        "severity_probabilities.prop_moderate_ED",
        "severity_probabilities.prop_moderate_GP",
        "severity_probabilities.prop_mild",
        "sampling.psampled",
        "sampling.sample_target_prob",
        "metagenomic_sensitivity.sensitivity_mg_virus",
        "metagenomic_sensitivity.sensitivity_mg_bacteria",
        "metagenomic_sensitivity.sensitivity_mg_fungi"
    ]
    
    for field in prob_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if !(0 <= value <= 1)
                push!(errors, "Probability field $field must be between 0 and 1, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    # Check that rates are positive
    rate_fields = [
        "healthcare_rates.gp_only_rate",
        "healthcare_rates.ed_direct_rate",
        "transmission.infectivity"
    ]
    
    for field in rate_fields
        keys = split(field, ".")
        value = config
        try
            for key in keys
                value = value[key]
            end
            if value <= 0
                push!(errors, "Rate field $field must be positive, got $value")
            end
        catch
            push!(warnings, "Could not validate field $field")
        end
    end
    
    return warnings, errors
end

"""
    get_parameter_mapping()

Get the mapping dictionary from config keys to parameter names.
"""
function get_parameter_mapping()
    return Dict(
        :contact_rates_fcont => :fcont,
        :contact_rates_gcont => :gcont,
        :contact_rates_oocont => :oocont,
        :day_of_week_scaling_dowcont => :dowcont,
        :transmission_infectivity => :infectivity,
        :transmission_infectivity_shape => :infectivity_shape,
        :transmission_infectivity_scale => :infectivity_scale,
        :disease_progression_latent_shape => :latent_shape,
        :disease_progression_latent_scale => :latent_scale,
        :disease_progression_infectious_shape => :infectious_shape,
        :disease_progression_infectious_scale => :infectious_scale,
        :transmission_reduction_rho_hosp => :ρ_hosp,
        :transmission_reduction_rho_asymptomatic => :ρ_asymptomatic,
        :network_dynamics_frate => :frate,
        :network_dynamics_grate => :grate,
        :severity_probabilities_prop_severe_hosp_short_stay => :prop_severe_hosp_short_stay,
        :severity_probabilities_prop_severe_hosp_long_stay => :prop_severe_hosp_long_stay,
        :severity_probabilities_prop_moderate_ED => :prop_moderate_ED,
        :severity_probabilities_prop_moderate_GP => :prop_moderate_GP,
        :severity_probabilities_prop_mild => :prop_mild,
        :healthcare_rates_gp_only_rate => :gp_only_rate,
        :healthcare_rates_ed_direct_rate => :ed_direct_rate,
        :healthcare_rates_gp_before_hosp_rate => :gp_before_hosp_rate,
        :healthcare_rates_ed_from_gp_rate => :ed_from_gp_rate,
        :healthcare_rates_hosp_admit_direct_rate => :hosp_admit_direct_rate,
        :healthcare_rates_hosp_admit_from_gp_rate => :hosp_admit_from_gp_rate,
        :healthcare_rates_hosp_recovery_rate => :hosp_recovery_rate,
        :healthcare_rates_hosp_short_stay_recovery_rate => :hosp_short_stay_recovery_rate,
        :healthcare_rates_hosp_long_stay_recovery_rate => :hosp_long_stay_recovery_rate,
        :healthcare_rates_hosp_death_rate => :hosp_death_rate,
        :healthcare_rates_triage_icu_rate => :triage_icu_rate,
        :healthcare_rates_icu_to_death_rate => :icu_to_death_rate,
        :healthcare_rates_icu_to_stepdown_leading_to_recovery_rate => :icu_to_stepdown_leading_to_recovery_rate,
        :healthcare_rates_stepdown_to_recovery_after_icu_rate => :stepdown_to_recovery_after_icu_rate,
        :sampling_psampled => :psampled,
        :sampling_sample_target_prob => :sample_target_prob,
        :sampling_turnaroundtime_icu => :turnaroundtime_icu,
        :sampling_turnaroundtime_rcgp => :turnaroundtime_rcgp,
        :sampling_turnaroundtime_hariss => :turnaroundtime_hariss,
        :sampling_icu_swab_lag_max => :icu_swab_lag_max,
        :discharge_timing_tdischarge_ed_upper_limit => :tdischarge_ed_upper_limit,
        :discharge_timing_tdischarge_hosp_short_stay_upper_limit => :tdischarge_hosp_short_stay_upper_limit,
        :imports_commuterate => :commuterate,
        :imports_importrate => :importrate,
        :imports_nimports => :nimports,
        :imports_import_t_df => :import_t_df,
        :imports_import_t_s => :import_t_s,
        :evolution_mu => :μ,
        :evolution_omega => :ω,
        :metagenomic_sensitivity_sensitivity_mg_virus => :sensitivity_mg_virus,
        :metagenomic_sensitivity_sensitivity_mg_bacteria => :sensitivity_mg_bacteria,
        :metagenomic_sensitivity_sensitivity_mg_fungi => :sensitivity_mg_fungi
    )
end

"""
    update_configurable_parameters(P, config::Dict)

Create a new parameter NamedTuple with updated values from config.
"""
function update_configurable_parameters(P, config::Dict)
    flat_config = flatten_dict(config)
    P_dict = Dict(pairs(P))
    param_mapping = get_parameter_mapping()
    
    updated_params = []
    for (config_key, param_key) in param_mapping
        if haskey(flat_config, config_key)
            # Handle tuple conversion for dowcont
            if param_key == :dowcont && isa(flat_config[config_key], Vector)
                P_dict[param_key] = tuple(flat_config[config_key]...)
            else
                P_dict[param_key] = flat_config[config_key]
            end
            push!(updated_params, param_key)
        end
    end
    
    if !isempty(updated_params)
        @info "Updated $(length(updated_params)) parameters: $(join(updated_params, ", "))"
    end
    
    return NamedTuple(P_dict)
end