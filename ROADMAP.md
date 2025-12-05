# Version 0.2.0 (implemented)

## Age structure 

- Polymod (or Warwick Uni Social Contact Survey) mixing matrix
    * distinct mixing for work/school, home, other, (and possibly travel) 
- Scale rates of :severe infection with age 
- Stratify hospital and ICU admission rates, and death rates

## Update model compartments
- Split 'Removed' compartments into 'Death' and 'Recovered'



# Version 0.3.0 

## Adapt care pathways
- Add new infection severity level to enable modelling of secondary care sampling

## Configuration files 
- For editing pathogen parameter input values and sampling parameters

## Standardise analaysis of sampling strategies
- Adapt analysis scripts to standard format enabling ease of use

## Incorporate analysis in R into Julia
- In particular the fitting of statistical distributions to adjust for right-censoring of simulation results
- This should make analysis easier

## Fit joint distributions 
- Enables adjustment for right-censoring when combining simulation results for sampling from more than one healthcare setting
- Possibly using a copula

## More realistic importation dynamics 
- Short-duration contacts at port-of-entry (similar to commuter dynamics)
- Home region should be drawn prop-to population density (or possibly derived from ONS International Passenger Survey)
    * simulate 1-way migration from port-of-entry to home region

## Long distance migration 
- Scalable migration matrix for occasional long-distance travel 

## Non-metagenomic surveillance model
- 'Death' compartment to be used as pathway for detections in non-metagenomic surveillance model
- Incorporate imported cases as potential trigger for detection in non-metagenomic surveillance model



# Version 0.4.0 

## Sample parameters 

- Defined ranges for covid-like and pandemic-flu-like 
- Pipeline for experiments: sampling parameters and simulating forests 




# Version 0.5.0 

## Phylogenetic tree simulation 

- Include 'nodes' in Infection for tracking MRCAs 
- Compute chronological and genetic distances between sample units 
- Tree figures using Phylo.jl 

## Maps 

- Density of imports, infections, clusters etc. 
