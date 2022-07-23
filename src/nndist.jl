# G. Froyland

using LinearAlgebra, Statistics, Distances

export nndist

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

    # if (mode==:plot)
    #     plot(nn_dist)
    # end

    return average_nn_dist

end