<!-- 
library(rmarkdown) 
render("parameters.md", output_format = "html_document")
render("~/mSCAPE/presentations/nbpmscape-slides.md", output_format = "html_document")

-->

# NBPMscape.jl
## Disease Transmission and Cluster Growth Simulation

---

## Parameter values
Here we describe the parameter values and estimates used in simulating disease transmission and cluster growth of a SARS-CoV-2-like respiratory pathogen using *NBMscape.jl*.

### Importation probability
We assume that the initial case of infection is imported from outside the UK and that imports continue after the initial case. We assume a daily importation rate of 0.5 and 1000 imports as per the analysis by du Plessis et al (2020) of the establishment of the SARS-CoV-2 epidemic in the UK. 

We also include t-distribution parameters from an analysis of lineages studied in Volz et al (2020).$\color{red}{\text{[HOW USED IN MODEL?]}}$

The location of the imported cases in the UK is determined using a probability based on international air passenger arrivals and airport location. International air passenger traffic is sourced from the [UK Civil Aviation Authority Table 12.1](https://www.caa.co.uk/data-and-analysis/uk-aviation-market/airports/uk-airport-data/) for 2024. [Accessed: 31 October 2024 and 15 March 2025]. Geographic area definitions from the UK Office for National Statistics (ONS) were used to determine the International Territory Level 2 (ITL2) in which airports in England are located.

```julia
importrate = 0.5
nimports = 1000
import_t_df = 2.48 $\color{red}{\text{[TO DO]}}$
import_t_s = 8.75 $\color{red}{\text{[TO DO]}}$
$\color{red}{\text{[TO DO - CAA import data]}}$
```

### Infectivity
During the first wave of SARS-CoV-2 in the UK, infectivity was related to viral load [$\color{red}{\text{[ref]}}$]. This varies with time during infection and can be modelled with a gamma distribution. We use the gamma distribution parameters estimated by Metcalf et al (2021): Gamma(shape = 2.2, scale = 2.5).
```julia
infectivity = 6.0        # Transmission rate scale
infectivity shape = 2.2  # Shape parameter for gamma distribution
infectivity_scale = 2.5  # Scale parameter for gamma distribution
```

### Latent period duration
The latent infection period is the duration between the time of infection and the time that an individual becomes infectious. This differs from the incubation period which is the duration between the time of infection and the onset of symptoms. These two periods may be different, which can result in pre-symptomatic transmission. The latent infection period for an individual is assumed to be drawn from a gamma distribution. We use the values estimated by Galmiche et al (2023): Gamma(shape = 4.24, scale = 1.08). $\color{red}{\text{[We use gamma distribution parameters (shape = 3.2638, scale = 0.979) that we estimate from the mean latent period estimate of 3.3 days (95% CI: 0.2, 7.9) by Zhao et al (2020).]}}$

```julia
# Galmiche et al (2023)
latent_shape = 4.24
latent_scale = 1.08

# Zhao et al (2020)
latent_shape = 3.2638
latent_scale = 0.979
```

### Infectious period duration
Individuals remain infectious for different periods of time. The length of the infectious period follows a gamma distribution. We use the parameters estimated by Verity et al (2020) which have a mean of 24 days: Gamma(shape = 8.16, scale = 3.03).

```julia
infectious_shape = 8.16
infectious_scale = 3.03

```

### Number of contacts
The number of contacts varies by individual. We analysed the POLYMOD UK survey [Massong et al, (2008)] to produce an empirical distribution of the number of contacts in various settings: workplace / school, other (non-home). We analysed the UK ONS 2021 census data for household size ([ONS nomis 2021 UK census TS017 Household size data](https://www.nomisweb.co.uk/query/construct/summary.asp?mode=construct&version=0&dataset=2037) [Accessed: 18 March 2025]) to produce an empirical distribution for the number of household contacts (= household size - 1). These were fit to statistical distributions using the *fitdistrplus* R package [Delignette-Muller & Dutang (2015)]. 

Household: Poisson(λ = 1.358975) [$\color{red}{\text{or Negative Binomial(size = 4.447099, p = 0.7659633}}$]

Workplace/school: Negative Binomial(r = 1.44, p = 0.1366)

Casual/other: Gamma(shape = 1.42, scale = 6.27746)


```julia
fnegbinomr = 4.447099   # Household contact negative binomial distribution r value
fnegbinomp = 0.7659633  # Household contact negative binomial distribution p value
# fpoissonrate = 1.358975 # Household contact Poisson distribution rate value
gnegbinomr = 1.44       # Workplace / school contact negative binomial distribution r value
gnegbinomp = 0.1366     # Workplace / school contact negative binomial distribution p value
oorateshape = 1.42      # Casual / other contact gamma distribution shape
ooratescale = 6.27746   # Casual / other contact gamma distribution scale
oorateshape1 = 2.42     # [$\color{red}{\text{ADD DETAILS}}$]
ooratescale1 = 6.27746  # [$\color{red}{\text{ADD DETAILS}}$]
```

### Changes to contact networks
The contact network has a dynamic structure and allows for contacts to be made and lost. The rate at which contacts are made and lost varies by the type of contact. [$\color{red}{\text{ref}}$]

```julia
frate = 0.0       # Rate of gaining household contacts
grate = 1 / 30.0  # Rate of gaining workplace / school contacts
```

### Commuting
The model incorporates movement of individuals between regions. The commuting rates are based on an analysis of 'origin-destination' data for England and Wales which is part of the ONS UK 2021 census and available through the ([nomis website](https://www.nomisweb.co.uk/sources/census_2021_od) [Accessed: 14 March 2025]). In particular, the 'ODWP01EW - Location of usual residence and place of work' dataset was used.
```julia
commuterate = 2.0 # [$\color{red}{\text{TO UPDATE?}}$]
```

### Relative probability of transmission
The probability of transmission is correlated with the distance between contacts and the duration spent with the contact [Ferreti et al, (2023)]. The POLYMOD UK survey [Massong et al, (2008)] contains information on the time spent with a contact and whether there was physical contact or not. This study also contains information on the frequency of contact. We accessed the survey data via the *socialmixr* R package [$\color{red}{\text{ref}}$]. We derived relative transmission probabilities based on this data [$\color{red}{\text{[TO DO]}}$].
```julia
fcont = 1.00 TODO # Household contacts
gcont = 0.50      # Workplace / school contacts
oocont = 0.50     # Casual / other contacts
```

### Day of the week
We include some temporal variation in contact rates. We analsyed contact patterns in the POLYMOD UK survey [Massong et al, (2008)] by the day of the week. These factors are used to alter transmission probability. 
```julia
#          Sunday     Monday     Tuesday    Wednesday  Thursday   Friday     Saturday
dowcont = (0.1043502, 0.1402675, 0.1735913, 0.1437642, 0.1596205, 0.1445298, 0.1338766)

```

### Infection severity
The severity of infection is used to determine an individual's care pathway and changes their probability of being sampled, either in the community or in an intensive care unit (ICU).
A WHO report published on 20 February 2020 [WHO, 2020] estimated that *"approximately **80% of laboratory confirmed patients have had mild to moderate disease**, which includes non-pneumonia and pneumonia cases,* **13.8% have severe disease** *(dyspnea, respiratory frequency ≥30/minute, blood oxygen saturation ≤93%, PaO2/FiO2 ratio <300, and/or lung infiltrates >50% of the lung field within 24-48 hours) and **6.1% are critical*** (respiratory failure, septic shock, and/or multiple organ dysfunction/failure)."* At that time asymptomatic infection was judged to be relatively rare.

Knock et al (2021) estimate the infection hospitalisation ratio at 2.55% (95% CrI: 2.17%-3.04%) for the period from 18 March to 31 May 2020. They also disaggregate this by age group and also estimate the probability of admission to ICU as a function of age group (5 year bands).

```julia
propmild = 0.60
propsevere = 0.05
```

### Rate for progressing through care pathway
In the model, infectees progress along a care pathway comprising three stages: visiting a GP, admission to hospital, and admission to ICU. We define the rates at which infected individuals progress along this pathway. How far the infected individual progresses along the pathway depends on the severity of infection, which is described above. 
```julia
gprate = 1/3 
hospadmitrate = 1/3
icurate = 1/5 
```

### Proportion of ICU cases sampled
only a proportion of infected individuals that are admitted to ICU will be sampled for metagenomic sequencing. This is a variable value in the model which can be tuned to investigate the impact on detection time for an epidemic outbreak.
```julia
psampled = 0.05
```

### Transmission reduction
We assume that the opportunity for the virus to be transmitted is reduced once indvidiuals are admitted to hospital and ICU as a result of a reduction in the contacts and the use of protective equipment. We assume this reduction is approximately 75% [$\color{red}{\text{ref}}$], and we therefore multiply by 0.25.

```julia
ρ = 0.25
```

### Time for metagenomic sequencing and analysis
We assume that the time duration between a sample being taken and a new variant or novel pathogen being flagged is between 3 and 7 days, from which we sample uniformly in the model. This includes metagenomic sequencing, bioinformatics and inclusion in the database. 

The duration estimate includes 6.7 hours for metagenomic sampling [Charalampous et al, 2024].
```julia
lagsseqdblb = 3
lagsseqdbub = 7
```

### Phylogenetic tree simulation
Future versions of the model will include 'nodes' in the Infection object for tracking MRCAs. Chronological and genetic distances between sample units will also be computed using the parameters below.
```julia
μ = 0.001 # mean clock rate -- additive relaxed clock
ω = 0.5   # variance inflation
```

### References
Charalampous et al. (2024), Routine Metagenomics Service for ICU Patients with Respiratory Infection, *Am J Respir Crit Care Med*, 209(2), pp 164–174. DOI: 10.1164/rccm.202305-0901OC

Delignette-Muller & Dutang (2015), fitdistrplus: An R Package for Fitting Distributions, *Journal of Statistical Software*, 64(4), 1–34. DOI:10.18637/jss.v064.i04

du Plessis et al. (2020), Establishment and lineage dynamics of the SARS-CoV-2 epidemic in the UK, *Science*, 371(6530), pp 708-712. DOI: 10.1126/science.abf2946

Galmiche et al. (2023), SARS-CoV-2 incubation period across variants of concern, individual factors, and circumstances of infection in France: a case series analysis from the ComCor study, *Lancet Microbe*, 4(6), e409-e417. DOI: 10.1016/S2666-5247(23)00005-8

Knock et al. (2021), Key epidemiological drivers and impact of interventions in the 2020 SARS-CoV-2 epidemic in England, *Sci. Transl. Med.*, 13,eabg4262. DOI: 10.1126/scitranslmed.abg4262

Metcalf et al. (2021), $\color{red}{\text{[ADD DETAILS]}}$

Mossong et al. (2008), Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases, PLoS Medicine, 5(3), e74. DOI: 10.1371/journal.pmed.0050074. Survey data at https://zenodo.org/records/3874557

Verity et al. (2020), $\color{red}{\text{[ADD DETAILS]}}$

Volz et al. (2020), Evaluating the Effects of SARS-CoV-2 Spike Mutation D614G on Transmissibility and Pathogenicity, *Cell*, 184(1), 64-75. DOI: 10.1016/J.CELL.2020.11.020

WHO (2020). Report of the WHO-China Joint Mission on coronavirus disease 2019 (COVID-19). Feb 28, 2020. https://www.who.int/publications-detail/report-of-the-who-china-joint-mission-oncoronavirus-disease-2019-(covid-19) (accessed 21 March 2025).
