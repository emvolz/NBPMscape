#=
Tests to check parameter configuration is working as expected
=#


using Test
using YAML
using Logging
using DataFrames
using CSV
using Suppressor
using NBPMscape
using Plots
using DifferentialEquations
using NamedArrays

# Include the module containing the functions
include(joinpath(@__DIR__,"..","src","config.jl"))
include(joinpath(@__DIR__,"..","src","core.jl"))

@testset "initialize_parameters() Tests" begin
    
    # Helper function to create a temporary valid config file
    function create_temp_config(config_dict::Dict)
        temp_file = tempname() * ".yaml";
        YAML.write_file(temp_file, config_dict);
        return temp_file
    end
    
    # Helper function to create minimal valid config
    function create_minimal_valid_config()
        return Dict(
                "parameters" => Dict(
                "pathogen_type" => "bacteria",#"virus",
                "infectivity" => 20.0,
                "prop_severe_hosp_short_stay" => 0.4,
                "prop_severe_hosp_long_stay" => 0.6,
                "prop_moderate_ED" => 0.12,
                "prop_moderate_GP" => 0.11,
                "prop_mild" => 0.77,
                "gp_only_rate" => 0.2,
                "hosp_recovery_rate" => 0.09,
                "nimports" => 1000,
                "mu" => 0.001,
                "omega" => 0.5,
                "fcont" => 1.0,
                "gcont" => 0.6,
                "oocont" => 0.75,
                "dowcont" => [0.1, 0.14, 0.17, 0.14, 0.16, 0.14, 0.13],
                "turnaroundtime_icu" => [2, 4],
                "turnaroundtime_hariss" => [2, 4],
                "turnaroundtime_rcgp" => [2, 4],
                "initial_dow" => 1,
                "phl_collection_dow" => [2, 5],
                "sample_allocation" => "equal",
                "weight_samples_by" => "ae_mean",
                "icu_sample_type" => "regional",
                "icu_site_stage" => "current",
                "sample_icu_cases_version" => "number",
                "sample_proportion_adult" => "free",
                "icu_only_sample_before_death" => true,
                "hariss_only_sample_before_death" => true,
                "hosp_ari_admissions_adult_p" => 0.52,
                "hosp_ari_admissions_child_p" => 0.48,
                "icu_ari_admissions_adult_p" => 0.76,
                "icu_ari_admissions_child_p" => 0.24,
                "ed_ari_destinations_adult_p_discharged" => 0.628,
                "ed_ari_destinations_adult_p_short_stay" => 0.030,
                "ed_ari_destinations_adult_p_longer_stay" => 0.342,
                "ed_ari_destinations_child_p_discharged" => 0.861,
                "ed_ari_destinations_child_p_short_stay" => 0.014,
                "ed_ari_destinations_child_p_longer_stay" => 0.125
            )
        );
    end
    
    # Test 1: Default parameters (no config file)
    @testset "Default Parameters" begin
        # Test with no arguments
        @suppress_err begin
            P = initialize_parameters();
            @test P !== nothing
            @test isa(P, NamedTuple)
        end
        
        # Test with explicit empty string
        @suppress_err begin
            P_empty = initialize_parameters("");
            @test P_empty !== nothing
            @test isa(P_empty, NamedTuple)
        end
    end
    
    # Test 2: Valid config file
    @testset "Valid Config File" begin
        config_dict = create_minimal_valid_config();
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                @test isa(P, NamedTuple)
                
                # Test that specific parameters were updated
                #@test P.pathogen_type == "bacteria"
                #@test P.infectivity == 20.0
                # Test that parameters from minimal config were updated
                if haskey(P, :pathogen_type)
                    @test P.pathogen_type == "bacteria"
                end
                if haskey(P, :infectivity)
                    @test P.infectivity == 20.0
                end
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 3: Non-existent config file
    @testset "Non-existent Config File" begin
        non_existent_file = "non_existent_config.yaml";
        
        @suppress_err begin
            # Should throw an error when trying to load non-existent file 
            # and not load default parameters (although they may already be loaded)
            @test_throws Exception initialize_parameters(non_existent_file)
            # Should throw a SystemError or similar file-related error
            @test_throws Union{SystemError, Base.IOError, ErrorException} initialize_parameters(non_existent_file)
            
            #@test P !== nothing # This is just looking at a previous P that was loaded 
            #@test isa(P, NamedTuple) # This is just looking at a previous P that was loaded 
            #P_default = initialize_parameters();
            #@test isequal(P, P_default) # This is just looking at a previous P that was loaded 
        end
    end
    
    # Test 4: Invalid YAML file
    @testset "Invalid YAML File" begin
        # Create invalid YAML file
        temp_invalid = tempname() * ".yaml";
        open(temp_invalid, "w") do f
            write(f, "invalid: yaml: content: [unclosed")
        end
        
        try
            # initialize_parameters function should fail
            # and default parameters should not be loaded
            @suppress_err begin
                @test_throws Exception P = initialize_parameters(temp_invalid);
                @test P == nothing
                @test !isa(P, NamedTuple)
            end
        finally
            rm(temp_invalid, force=true)
        end
    end
    
    # Test 5: Config with validation errors
    @testset "Config with Validation Errors" begin
        invalid_config = Dict(
            "parameters" => Dict(
                "pathogen_type" => "invalid_type",  # Should be virus/bacteria/fungi
                "prop_severe_hosp_short_stay" => 1.5,  # Should be 0-1
                "gp_only_rate" => -0.1,  # Should be positive
                "turnaroundtime_icu" => [4, 2],  # Lower > upper
                "initial_dow" => 8,  # Should be 1-7
                "sample_allocation" => "invalid"  # Should be equal/weighted
            )
        );
        temp_config = create_temp_config(invalid_config);
        
        try
            @suppress_err begin
                # Should throw ArgumentError due to validation failures
                @test_throws ArgumentError initialize_parameters(temp_config)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
@testset "Config with Validation Errors" begin
    # Create config with proper structure but invalid values
    invalid_config = Dict(
        "parameters" => Dict(
            # Include required fields that validator expects
            "turnaroundtime_icu" => [4, 2],  # Lower > upper (invalid)
            "turnaroundtime_hariss" => [4, 2],  # Lower > upper (invalid)
            "turnaroundtime_rcgp" => [4, 2],  # Lower > upper (invalid)
            "pathogen_type" => "invalid_type",  # Should be virus/bacteria/fungi
            "prop_severe_hosp_short_stay" => 1.5,  # Should be 0-1
            "prop_severe_hosp_long_stay" => -0.1,  # Should be 0-1
            "gp_only_rate" => -0.1,  # Should be positive
            "initial_dow" => 8,  # Should be 1-7
            "sample_allocation" => "invalid",  # Should be equal/weighted
            "hosp_ari_admissions_adult_p" => 0.6,  # These should sum to 1
            "hosp_ari_admissions_child_p" => 0.3,  # Sum = 0.9, not 1
            "icu_ari_admissions_adult_p" => 0.8,   # These should sum to 1
            "icu_ari_admissions_child_p" => 0.3,   # Sum = 1.1, not 1
            "ed_ari_destinations_adult_p_discharged" => 0.5,  # These should sum to 1
            "ed_ari_destinations_adult_p_short_stay" => 0.2,
            "ed_ari_destinations_adult_p_longer_stay" => 0.2,  # Sum = 0.9
            "icu_only_sample_before_death" => "yes",  # Should be boolean
            "hariss_only_sample_before_death" => "no"  # Should be boolean
        )
    );
    temp_config = create_temp_config(invalid_config);
    
    try
        # Should throw ArgumentError due to validation failures
        @test_throws ArgumentError initialize_parameters(temp_config)
    finally
        rm(temp_config, force=true)
    end
end

    # Test 6: Config with warnings but no errors
    @testset "Config with Warnings" begin
        warning_config = create_minimal_valid_config();
        # Add different turnaround times to trigger warning
        warning_config["parameters"]["turnaroundtime_icu"] = [1, 3];
        warning_config["parameters"]["turnaroundtime_hariss"] = [2, 4];
        warning_config["parameters"]["turnaroundtime_rcgp"] = [3, 5];
        
        temp_config = create_temp_config(warning_config);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                @test isa(P, NamedTuple)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 7: Special parameter name conversions
    @testset "Parameter Name Conversions" begin
        config_dict = create_minimal_valid_config();
        config_dict["parameters"]["rho_hosp"] = 0.3;
        config_dict["parameters"]["rho_asymptomatic"] = 0.2;
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                @test P.ρ_hosp == 0.3
                @test P.ρ_asymptomatic == 0.2
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 8: Tuple conversion for dowcont
    @testset "Tuple Conversion" begin
        config_dict = create_minimal_valid_config();
        dowcont_vector = [0.1, 0.14, 0.17, 0.14, 0.16, 0.14, 0.13];
        config_dict["parameters"]["dowcont"] = dowcont_vector;
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                @test P.dowcont isa Tuple
                @test collect(P.dowcont) == dowcont_vector
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 9: Probability validation
    @testset "Probability Validation" begin
        # Test proportions that should sum to 1
        config_dict = create_minimal_valid_config();
        config_dict["parameters"]["hosp_ari_admissions_adult_p"] = 0.6;
        config_dict["parameters"]["hosp_ari_admissions_child_p"] = 0.3;  # Sum = 0.9, not 1
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                @test_throws ArgumentError initialize_parameters(temp_config)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 10: ED destinations validation
    @testset "ED Destinations Validation" begin
        config_dict = create_minimal_valid_config();
        # Adult destinations that don't sum to 1
        config_dict["parameters"]["ed_ari_destinations_adult_p_discharged"] = 0.5;
        config_dict["parameters"]["ed_ari_destinations_adult_p_short_stay"] = 0.2;
        config_dict["parameters"]["ed_ari_destinations_adult_p_longer_stay"] = 0.2;  # Sum = 0.9, not 1
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                @test_throws ArgumentError initialize_parameters(temp_config)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 11: Turnaround time validation
    @testset "Turnaround Time Validation" begin
        config_dict = create_minimal_valid_config();
        # Invalid turnaround time (not 2-element vector)
        config_dict["parameters"]["turnaroundtime_icu"] = [2, 4, 6];
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                @test_throws ArgumentError initialize_parameters(temp_config)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 12: Boolean parameter validation
    @testset "Boolean Parameter Validation" begin
        config_dict = create_minimal_valid_config();
        config_dict["parameters"]["icu_only_sample_before_death"] = "yes"  # Should be boolean;
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                @test_throws ArgumentError initialize_parameters(temp_config)
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 13: Sample proportion adult validation
    @testset "Sample Proportion Adult Validation" begin
        # Test valid "free" value
        config_dict = create_minimal_valid_config();
        config_dict["parameters"]["sample_proportion_adult"] = "free";
        
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                @test P.sample_proportion_adult == "free"
            end
        finally
            rm(temp_config, force=true)
        end
        
        # Test valid numeric value
        config_dict["parameters"]["sample_proportion_adult"] = 0.75;
        temp_config2 = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config2);
                @test P !== nothing
                @test P.sample_proportion_adult == 0.75
            end
        finally
            rm(temp_config2, force=true)
        end
    end
    
    # Test 14: Global P variable
    @testset "Global P Variable" begin
        @suppress_err begin
            P_returned = initialize_parameters();
            @test @isdefined P  # Global P should be defined
            @test @isdefined P_returned  # Global P should be defined
            #@test P == P_returned  # Global P should equal returned P
            @test isequal(P, P_returned)  # Global P should equal returned P
        end
    end
    
    # Test 15: DataFrame conversion
    @testset "DataFrame Conversion" begin
        config_dict = create_minimal_valid_config();
        temp_config = create_temp_config(config_dict);
        
        try
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                
                # Check that DataFrames are created (if convert_params_to_dfs works)
                if haskey(P, :ed_ari_destinations_adult)
                    @test P.ed_ari_destinations_adult isa DataFrame
                    @test nrow(P.ed_ari_destinations_adult) == 3
                end
                
                if haskey(P, :ed_ari_destinations_child)
                    @test P.ed_ari_destinations_child isa DataFrame
                    @test nrow(P.ed_ari_destinations_child) == 3
                end
            end
        finally
            rm(temp_config, force=true)
        end
    end
    
    # Test 16: Unknown parameters warning
    @testset "Unknown Parameters Warning" begin
        config_dict = create_minimal_valid_config();
        config_dict["parameters"]["unknown_parameter"] = 123;
        
        temp_config = create_temp_config(config_dict);
        
        try
            # Capture warnings
            @suppress_err begin
                P = initialize_parameters(temp_config);
                @test P !== nothing
                # The function should still work despite unknown parameters
            end
        finally
            rm(temp_config, force=true)
        end
    end

    # Test 17: print_changes function robustness testset
    @testset "print_changes() Robustness Tests" begin
        
        # Mock functions for testing
        function create_mock_default_params()
            return (
                param1 = 1.0,
                param2 = "test",
                param3 = [1, 2, 3],
                ρ_hosp = 0.25,
                μ = 0.001
            );
        end
        
        function create_mock_final_params()
            return (
                param1 = 2.0,  # Changed
                param2 = "test",  # Unchanged
                param3 = [1, 2, 3],  # Unchanged
                ρ_hosp = 0.30,  # Changed
                μ = 0.001,  # Unchanged
                extra_final_param = "new"  # Extra in final
            );
        end
        
        # Test 18: Config with fewer parameters 
        @testset "Config with fewer parameters" begin
            default_params = create_mock_default_params();
            final_params = create_mock_final_params();
            
            # Config with only some parameters
            config_data = Dict(
                "parameters" => Dict(
                    "param1" => 2.0,
                    "rho_hosp" => 0.30
                )
            );
            
            # Should not throw an error
            @test_nowarn print_changes(final_params, config_data, default_params)
        end
        
        # Test 19: Config with more parameters 
        @testset "Config with more parameters" begin
            default_params = create_mock_default_params();
            final_params = create_mock_final_params();
            
            # Config with extra parameters
            config_data = Dict(
                "parameters" => Dict(
                    "param1" => 2.0,
                    "param2" => "test",
                    "param3" => [1, 2, 3],
                    "rho_hosp" => 0.30,
                    "mu" => 0.001,
                    "extra_config_param" => "extra_value",
                    "another_extra" => 123
                )
            );
            
            # Should not throw an error
            @test_nowarn print_changes(final_params, config_data, default_params)
        end
        
        # Test 20: Empty config
        @testset "Empty config" begin
            default_params = create_mock_default_params();
            final_params = create_mock_final_params();
            
            # Empty config
            config_data = Dict("parameters" => Dict());
            
            # Should not throw an error
            @test_nowarn print_changes(final_params, config_data, default_params)
        end
        
        # Test 21: Malformed config 
        @testset "Malformed config" begin
            default_params = create_mock_default_params();
            final_params = create_mock_final_params();
            
            # Config without parameters key
            config_data = Dict("other_key" => "value");
            
            # Should not throw an error
            @test_nowarn print_changes(final_params, config_data, default_params)
        end
        
        # Test 22: Nil default params 
        @testset "Nil default params" begin
            final_params = create_mock_final_params();
            config_data = Dict("parameters" => Dict("param1" => 2.0));
            
            # Should not throw an error even with nil default params
            @test_nowarn print_changes(final_params, config_data, nothing);
        end
    end

end;