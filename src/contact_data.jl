### Load contact distributions
# Generated using socialmixr R package and UK POLYMOD survey data
# Artificial split by single age
#contact_distributions_age_single_yr = load( joinpath( @__DIR__, "..", "data", "polymod_contact_age_setting_distribution.rds" ) )
contact_distributions_age_group = load( joinpath( @__DIR__, "..", "data", "contact_setting_age_group_distributions.rds" ) )
# Remove the 'all' age group as not required
const CONTACT_DISTRIBUTIONS = filter(:age_group => !=("all"), contact_distributions_age_group)

### Load contact matrices (disaggregated by age group)
# Generated using socialmixr R package and UK POLYMOD survey data
const CONTACT_MATRIX_AGE_GROUPS = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_age_groups.rds" ) )
#contact_matrix_home_df = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_home.rds" ) , convert = true )
contact_matrix_home = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_home.rds" ) , convert = true )
const CONTACT_MATRIX_HOME = NamedArray( Matrix( contact_matrix_home ) 
                                      , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS))
contact_matrix_school_work = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_school_work.rds" ) )
const CONTACT_MATRIX_SCHOOL_WORK = NamedArray( Matrix( contact_matrix_school_work ) 
                                              , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS) )
contact_matrix_other = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_other.rds" ) )
const CONTACT_MATRIX_OTHER = NamedArray( Matrix( contact_matrix_other ) 
                                        , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS) )
