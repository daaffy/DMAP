using Plots
using CSV
using DataFrames

# -------------------------------------------------------------------------------------------------------
# FUNCTIONS

function random_subset(list,subset_length)
    # subset_length > length(list) ??
    global inds = Vector{Int64}(undef,subset_length)
    i = 1
    while (i <= subset_length)
        next = true
        inds[i] = rand(1:length(list))
        println(inds[i])
        for check_j = 1:(i-1)
            if (inds[i]==inds[check_j])
                next = false
                break
            end
        end
        if (next)
            i = i + 1
        end
    end
    return list[inds]
end

# -------------------------------------------------------------------------------------------------------

# uncomment to load into workspace
# df = CSV.read("/Users/jackh/Documents/Honours/DMAp/simple_pipe/time_0_1.csv", DataFrame)
# df = Matrix(df)

n_points = size(df,1)


# extract discrete vector field
X = Array{Float64}(undef,0,3)
dX = Array{Float64}(undef,0,3)
for i = 1:n_points
    global X = cat(X, df[i,2:4]', dims=1)
    global dX = cat(dX, df[i,5:7]', dims=1)
end



display(scatter(X[:,3],X[:,2]))