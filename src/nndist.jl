# G. Froyland

using LinearAlgebra, Statistics, Distances, Plots

export nndist, sp_nndist

function sp_nndist(X;mode=:none,metric=Euclidean())
    n=size(X,1)
    d=size(X,2)

    avg = 0
    curr_pairwise = Vector{Float64}(undef,n)
    curr_sort = similar(curr_pairwise)
    nn_dist_list = []
    for i = 1:n # step through incrementally rather than compute all at once
        curr_pairwise .= colwise(metric,X[i,:].*ones(d,n),X')
        curr = sort(curr_pairwise)[2d+1] # five-point stencil
        append!(nn_dist_list,curr)
        avg = avg + curr 
    end

    # if (mode==:plot)
    #     display(histogram(nn_dist_list))
    # end

    return avg/n, nn_dist_list
end


function nndist(X;metric=Euclidean(),mode=:none,D=pairwise(metric,X,X,dims=1))

    # Input: array X is nxd;  n data points and d dimensions
    # Output: average distance to the (2d)^th nearest neighbour (averaged over all points)

    n=size(X,1)
    d=size(X,2)

    

    # calculate full pairwise distance array (nxn)
    # D=[metric(X[i,:],X[j,:]) for i=1:n, j=1:n]
    # D = pairwise(metric,X,X,dims=1)

    # find distance from each point to its (2d)^th nearest neighbour
    nn_dist = [sort(D[:,i])[2d+1] for i=1:n]

    # take the mean of these distances
    average_nn_dist=mean(nn_dist)

    if (mode==:plot)
        display(histogram(nn_dist))
    end


    return average_nn_dist

end