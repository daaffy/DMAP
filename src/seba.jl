# G. Froyland, 2022

using LinearAlgebra

function SEBA(V, Rinit = nothing)

    # Inputs: 
    # V is pxr matrix (r vectors of length p as columns)
    # Rinit is an (optional) initial rotation matrix.

    # Outputs:
    # S is pxr matrix with columns approximately spanning the column space of V
    # R is the optimal rotation that acts on V, which followed by thresholding, produces S

    # Define soft-threshold function:  soft threshold scalar z by threshold μ
    function soft_threshold(z, μ)
        tau = sign(z) * max(abs(z) - μ, 0)
        return tau
    end

    # Begin SEBA algorithm

    maxiter = 5000   # maximum number of iterations allowed
    F = qr(V) # Enforce orthonormality
    V = Matrix(F.Q)
    p, r = size(V)
    μ = 0.99 / sqrt(p)

    S = zeros(size(V))
    # Perturb near-constant vectors
    for j = 1:r
        if maximum(V[:, j]) - minimum(V[:, j]) < 1e-14
            V[:, j] = V[:, j] .+ (rand(p, 1) .- 1 / 2) * 1e-12
        end
    end

    # Initialise rotation
    if Rinit ≡ nothing
        Rnew = I
    else
        # Ensure orthonormality of Rinit
        F = svd(Rinit)
        Rnew = F.U * F.Vt
    end

    #preallocate matrices
    R = zeros(r, r)
    Z = zeros(p, r)
    Si = zeros(p, 1)

    iter = 0
    while norm(Rnew - R) > 1e-14 && iter < maxiter
        iter = iter + 1
        R = Rnew
        Z = V * R'
        # Threshold to solve sparse approximation problem
        for i = 1:r
            Si = soft_threshold.(Z[:, i], μ)
            S[:, i] = Si / norm(Si)
        end
        # Polar decomposition to solve Procrustes problem
        F = svd(S' * V, full = false)
        Rnew = F.U * F.Vt
    end

    # Choose correct parity of vectors and scale so largest value is 1
    for i = 1:r
        S[:, i] = S[:, i] * sign(sum(S[:, i]))
        S[:, i] = S[:, i] / maximum(S[:, i])
    end

    # Sort so that most reliable vectors appear first
    ind = sortperm(vec(minimum(S, dims = 1)), rev = true)
    S = S[:, ind]

    return S, R

end
