using Plots, Distances
include("nndist.jl")

export findϵ

function findϵ(X)

    #X is an Nxd array of N points in ℝ^d
    N = size(X, 1)

    #cheap way to compute distances, you can use a faster way
    # D = [norm(X[i, :] - X[j, :]) for i = 1:N, j = 1:N]
    D = pairwise(Euclidean(),X,X,dims=1)

    baseϵ = nndist(X,D=D)
    #adjust the numbers 10 below to something larger if required (if you don't see flat sections at the start and end)nn
    loggridϵ = -20*log(2)+log(baseϵ):log(2):10*log(2)+log(baseϵ)
    #S=Vector{Array{N,N}}
    logsumvec = []
    S = similar(D)
    for logϵ = loggridϵ
        #S is your matrix P using exp(logϵ) as the ϵ
        S .= exp.(-D .^ 2 ./ exp(logϵ))
        # log_sum = log(sum(exp.(-D .^ 2 ./ exp(logϵ))))
        push!(logsumvec, log(sum(S))) 
        # push!(logsumvec, log_sum)
    end

    print("log(eps_default) = "*string(log(baseϵ/sqrt(2))))



    plot(loggridϵ, logsumvec)
    display(scatter!(loggridϵ, logsumvec, xlabel="log(ϵ)", ylabel="log(sum of P)"))

    
end

