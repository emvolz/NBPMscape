#= 
Applying Block model to correlations between median time to detection (TD) by ITL2 region

Correlation matrices are computed using 'median_TD_by_region.jl' and saved to a file

=#

# Load packages
using CSV
using DataFrames
using Pkg
#Pkg.add("Clustering")
using Clustering
using Plots
#Pkg.add("StatsPlots")
using StatsPlots
using Statistics
using LinearAlgebra



# Load correlation matrix
sim_reps = 15000
psampled = 1 #0.25#1#0.25
corrmat_df = CSV.read("scripts/median_TD_by_region/sim_regionentry/$(sim_reps)_replicates/psampled_$(Int(psampled*100))/regional_TD_correlation_matrix_psampled_$(Int(psampled*100)).csv", DataFrame)
corrmat = Matrix(corrmat_df[:,2:end])

# Two different methods used to assign regions to blocks:
# (1) Hierarchical Clustering
# (2) K-means Clustering

# Different numbers of blocks to assign regions to
num_blocks_list = collect(1:39)

# Create df to hold block allocations for different numbers of blocks
column_names = Symbol.(string.(1:39))
empty_columns = [fill(missing, 39) for _ in 1:39]

hc_block_allocation_df = DataFrame(; 
    ITL2_region = corrmat_df[:,1], 
    [name => col for (name, col) in zip(column_names, empty_columns)]...
)

hc_block_allocation_abs_df = hc_block_allocation_df
kmc_block_allocation_df = hc_block_allocation_df
kmc_block_allocation_abs_df = hc_block_allocation_df

# (1) Hierarchical Clustering
# Convert to distance matrix: d_ij = 1 - corr_ij
#distmat_abs = 1 .- abs.(corrmat)
distmat = 1 .- corrmat 
# Perform hierarchical clustering using average linkage
hclust_result_abs = hclust(distmat_abs, linkage=:average)
hclust_result = hclust(distmat, linkage=:average)
#plot(hclust_result, xticks=false, yticks=true, size=(800, 300))

# (2) K-means Clustering
# Create feature vectors for clustering (e.g. each variable's correlation profile)
#features_abs = abs.(corrmat)
features = corrmat

# Loop through the different block sizes and assign block numbers
# using the two different methods
for i in 1:length(num_blocks_list)
    
    # Number of blocks to split regions into
    num_blocks = num_blocks_list[i]

    # (1) Hierarchical Clustering
    # block_assignments: vector giving the block for each variable
    #hc_block_allocation_abs_df[!,i+1] = cutree(hclust_result_abs, k=num_blocks)
    hc_block_allocation_df[!,i+1] = cutree(hclust_result, k=num_blocks)
    
    # (2) K-means Clustering
    # Cluster using K-means on the rows of features
    #kmeans_result_abs = kmeans(features_abs, num_blocks; maxiter=100, display=:none)
    kmeans_result = kmeans(features, num_blocks; maxiter=100, display=:none)
    # Assignments of each region to blocks
    #kmc_block_allocation_abs_df[!,i+1] = kmeans_result_abs.assignments
    kmc_block_allocation_df[!,i+1] = kmeans_result.assignments

end

# Check whether the Hierarchical Clustering and K-means Clustering methods yield the same results
hc_block_allocation_df == kmc_block_allocation_df
hc_block_allocation_abs_df == kmc_block_allocation_abs_df
hc_block_allocation_abs_df == hc_block_allocation_df
kmc_block_allocation_abs_df == kmc_block_allocation_df
# Whether we use the absolute correlation or not doesn't seem to impact the block allocation
# (at least for 5000 sim replicates and ICU sampling of 25% and 100%)
# Not sure why because abs(cor) loses difference between positive and negative correlations

#### Perform Block correlation
function block_model_matrix(corrmat::Matrix{Float64}, blocks::Vector{Int})
    n = size(corrmat, 1)
    nb = maximum(blocks)  # number of blocks
    blockmat = zeros(nb, nb)

    # Calculate mean correlations within/between blocks
    for i in 1:nb, j in 1:nb
        inds_i = findall(x -> x == i, blocks)
        inds_j = findall(x -> x == j, blocks)
        # For within-block, avoid diagonal double-counting
        if i == j
            vals = [corrmat[k,l] for k in inds_i, l in inds_j if k != l]
        else
            vals = [corrmat[k,l] for k in inds_i, l in inds_j]
        end
        blockmat[i,j] = mean(vals)
    end

    # Reconstruct the matrix with block-constant values
    result = zeros(n, n)
    for i in 1:n, j in 1:n
        result[i, j] = blockmat[blocks[i], blocks[j]]
    end
    for i in 1:n
        result[i, i] = 1.0  # Ensure diagonals remain 1
    end
    result
end

# Create array to store individual matrices for each number of blocks
block_model_matrices = Vector{Matrix{Float64}}(undef, 39) #Vector{Matrix{Float64}}()
#ndims(block_model_array)

# Loop through number of blocks, calculate block matrix based on allocated blocks
# and add resulting matrix to array 
for j in 2:40
    # IF THERE ARE DIFFERENCES BETWEEN METHODS THEN THIS DF MAY NEED TO BE CHANGED
    blocks = hc_block_allocation_df[:,j]   
    block_model_matrices[j-1] = block_model_matrix(corrmat, blocks)
end

#### Plot correlations and list block constituents
# Choose number of blocks to split regions into
n_blocks = 10
# Create region codes/names for labelling heatmap
x_names = names(corrmat_df)[2:end]
y_names = corrmat_df[:,1]

# Create correlation values for annotating the heatmap (if required)
anns = []
n=39
for i in 1:n
    for j in 1:n
        val = round(corrmat[i, j], digits=2)
        # Note: y-axis is reversed in Plots.jl heatmap to match matrix layout
        push!(anns, (j, n - i + 1, string(val)))
    end
end

# Plot heatmap
heatmap(  x_names # x-tick labels (columns)
        , xmirror = true # move x-axis labels to the top
        , xrotation = 90 # Rotate x-axis labels
        #, xticks = (1:n, names(td_by_simrep_by_region))
        , xticks = (collect(1:n) .- 0.5, x_names) #names(td_by_simrep_by_region)) # centers of each cell
        , y_names #names(td_by_simrep_by_region) # y-tick labels (rows)
        #, yticks = (1:n, names(td_by_simrep_by_region))
        , yticks = (collect(1:n) .- 0.5, y_names) #names(td_by_simrep_by_region))# centers of each cell
        , block_model_matrices[ n_blocks ];              # the correlation matrix
          color = :RdBu        # diverging color map, blue to red
        , clim = (-1, 1)       # color limits for correlation coefficients
        , xlabel = " \n"#"ITL2 region"
        , ylabel = "ITL2 region"
        , title = "Correlation Matrix with $(n_blocks) blocks for regional TD \n across $(sim_reps) simulation replicates\n at ICU sampling rate of $(psampled*100)%"
        , yflip = true         # flips y-axis to match matrix layout
        , size = (1100, 800)
        #, annotations = anns
)

# List block constituents
# Create temporary df containing block allocation for n blocks
region_block_allocation_df = DataFrame( ITL2_region = corrmat_df[!,1]
                                        , Block = hc_block_allocation_df[!,n_blocks+1] 
                                      )

region_block_allocation_df = sort(region_block_allocation_df, :Block)
println(region_block_allocation_df)

CSV.write("scripts/median_TD_by_region/sim_regionentry/$(n_replicates)_replicates/psampled_$(Int(psampled*100))/regional_TD_correlation_block_allocation_psampled_$(Int(psampled*100)).csv"
            , region_block_allocation_df)

# Group the regions by allocated block number 
gdf = groupby(region_block_allocation_df, :Block) #+1 because column 1 is the region name

# To see each group:
for subdf in gdf
    println("Regions in block $(subdf.:Block[1]):")
    println(subdf)
end