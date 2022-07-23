using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis
using SparseArrays
using Arpack
using StaticArrays
using ProfileView
using BenchmarkTools
using NearestNeighbors

"""
    # Rotating double gyre flow field.
For 0 <= t <= 1.
"""

s = t -> (t^2)*(3-2*t)
function vec_field!(du,u,p,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
res = 100
X_0 = create_grid([0 1],[0 1],[res res])
t = (0:0.1:1)

@time traj = solve_ensemble(vec_field!,X_0,t,s)

# -------------------------------------------------------------------------------------------------------
# construct adjacency and diagonal matrix 
@time A, D = npg_k(traj, k=5) # use modified k-nearest neighbours approach
# ϵ = nndist(traj.X[:,1,:])
@time A, D = npg_ϵ(traj, ϵ=2*ϵ)

# -------------------------------------------------------------------------------------------------------
# calculate eigens 
r = 2
@time λ, v = eigs(D-A, D, nev=r+1, which=:SM) # sparse method; transpose converts sparse matrix back into a dense matrix?
# v = sqrt.(D)*v # this is done inside eig method above?

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing

V = real(v[:,1:r])
S, R = SEBA(V)
# S = V

# -------------------------------------------------------------------------------------------------------
# # plot
t_i = 10
# display(scatter(traj.X[:,t_i,1], traj.X[:,t_i,2], zcolor=real(S[:,1]), c=:magma, aspectratio=1))
# display(scatter(traj.X[:,t_i,1], traj.X[:,t_i,2], zcolor=real(S[:,2]), c=:magma, aspectratio=1))

display(heatmap(reshape(real(S[:,1]),res,res)',c=:magma))
display(heatmap(reshape(real(S[:,2]),res,res)',c=:magma))