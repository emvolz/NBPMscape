#- Simulate multiple replicates using default covid-like parameters 
#- Run for final analysis for research article
#-- after changes made following code review 
#-- importmodel=:GLEAM_outbreak_country_not_specified in simforest()
#-- covid-like parameters
#-- updated jump rates to account for competing pathways to ED/hosp (direct or via GP)

#using Pkg
#Pkg.add("JLD2")
#Pkg.add("DataFrames")
#Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="536099580c3b6868fdc02c1c2cb64b639681a077") # 'master' branch committed at 16:29 08 Oct 2026 (merge pull request #8)
using Pkg
using JLD2
using DataFrames
using NBPMscape

println("Julia version: ", VERSION)
println("Project:")
# To check that using the local Project.toml and Manifest.toml rather than the version associated with the Julia version (which has different package versions that are incompatible)
println("Active project: ", Base.active_project())
Pkg.status("NBPMscape")

println("Starting analysis...")

initialize_parameters("config/outbreak_params_covid19_like_symptomatic_base_minus_20_20261008.yaml"); # Updated jump rates for care pathways

NREPS = 10
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=70.0, maxgenerations=100, max_cases=100000000, importmodel=:GLEAM_outbreak_country_not_specified), 1:NREPS )
@save "covidlike-bm20-1.6.1-sims-nrep10.jld2" sims

### Also save versions filtered for G and cases accessing healthcare, i.e. GP, ED, hospital and ICU (reduces file size and more likely to be able to reload)

# Create vectors to store filtered dataframes: G filtered for ICU and GP rows
#sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]
#sims_G_gp_filter = [DataFrame() for _ in 1:length(sims)]
# Saving as a single file avoids some overlap of cases, for example where the individual visits the GP before ICU
sims_G_filtered = [DataFrame() for _ in 1:length(sims)]
# Create df to save the max_cases_df output - useful for checking the impact on the number of infections
sims_max_cases_df = [DataFrame() for _ in 1:length(sims)]
# Create df to save the simid and tinf values for each simulation replicate
sims_simid_tinf_df = [DataFrame() for _ in 1:length(sims)]
# Create df to save the tinf values for each simulation replicate
sims_tinf_df = [Vector() for _ in 1:length(sims)]

# Loop through sim replicates
for s in 1:length(sims)	

	try
		fo = sims[s]
		mc = fo.max_cases_df
		
		G = fo.G

		if size(mc,1) > 0
			sims_max_cases_df[s] = mc
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])	
		else
			sims_max_cases_df[s] = missing
			sims_simid_tinf_df[s] = missing
		end

		if size(G,1) > 0
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])
		else
			sims_simid_tinf_df[s] = missing
		end
		
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end

	try

		fo = sims[s]

		# Filter for G, which is df containing infection information
		# and only retain cases (rows) that progressed to healthcareICU (i.e. capable of detection under primary care, secondary care or ICU sampling methodology)
		G_filtered = fo.G[ isfinite.(fo.G.tgp) .| isfinite.(fo.G.ted) .| isfinite.(fo.G.thospital) .| isfinite.(fo.G.ticu), : ]
		# Trim columns to save on storage space
		G_filtered = G_filtered[:,["pid"
						,"tinf", "tgp", "ted", "thospital", "ticu", "tstepdown", "tdischarge", "trecovered", "tdeceased"
						, "severity", "fatal", "iscommuter"
						, "homeregion"
						#, "commuteregion"
						#, "generation"
						#, "F", "G", "H"
						, "infector_age", "infectee_age"
						, "importedinfection"
						, "simid"
						]]
		
		
		if size(G_filtered,1) > 0
			sims_G_filtered[s] = G_filtered
		else
			sims_G_filtered[s] = missing
		end
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end
	
	# Save times of infection directly from G. 
	# sims_simid_tinf_df only includes unique combinations of simid and tinf so if more than one infection happens at the 
	# same time in the same simtree then only one will be recorded
	try
		fo = sims[s]
		tinfs = fo.G[:,:tinf]
		if length(tinfs) > 0 
			sims_tinf_df[s] = tinfs
		else
			sims_tinf_df[s] = missing
		end
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end
	
end

#Save filtered data to files
#@save "covidlike-1.3.4-sims-nrep1000_filtered_G_icu.jld2" sims_G_icu_filter
#@save "covidlike-1.3.4-sims-nrep1000_filtered_G_gp.jld2" sims_G_gp_filter
@save "covidlike-bm20-1.6.1-sims-nrep10_filtered_G.jld2" sims_G_filtered
@save "covidlike-bm20-1.6.1-sims-nrep10_max_cases_df.jld2" sims_max_cases_df
@save "covidlike-bm20-1.6.1-sims-nrep10_sims_simid_tinf_df.jld2" sims_simid_tinf_df
@save "covidlike-bm20-1.6.1-sims-nrep10_sims_tinf_df.jld2" sims_tinf_df

println("Script finished successfully.")
