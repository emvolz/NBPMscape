#= Miscellaneous functions

- allocate_with_rounding:   allocates a number across a number of categories based on weights
                            ensuring integer values are allocated and the sum of allocations
                            is equal to the original total, e.g. total number of samples
                            allocated across NHS Trusts but the allocations must be integer values
                            and the sum must be equal to the total

=#

"""
Function:       allocate_with_rounding

Description:    Allocates a number across a number of categories/groups based on weights
                ensuring integer values are allocated and the sum of allocations
                is equal to the original total, e.g. total number of samples
                allocated across NHS Trusts but the allocations must be integer values
                and the sum must be equal to the total.

Arguments:      total::Int          Number to be allocated across categories/groups
                weights::Vector     Weightings for each category/group

Returns:        Vector of integer values

Examples:       alloc = allocate_with_rounding( total = 10, weights = [0.15, 0.15, 0.3, 0.4])
                # Checks
                println(alloc)  
                sum(alloc)
                alloc / sum(alloc)
"""
function allocate_with_rounding(;total, weights)
    weightsum = sum(weights)
    # Ideal unrounded allocations
    exact = total .* (weights ./ weightsum)
    # Integer part (floor)
    allocation = floor.(Int, exact)
    # Compute total remainder after integer (floor) allocation
    remainder = total - sum(allocation)
    # Fractional remainders
    fractional = exact .- allocation
    # Find indices of categories/groups with the largest fractional remainders (for distributing leftover units)
    idx = partialsortperm(fractional, rev=true, 1:remainder)
    # Add 1 to the allocations of categories/groups with the largest fractional parts
    allocation[idx] .+= 1
    return allocation
end

