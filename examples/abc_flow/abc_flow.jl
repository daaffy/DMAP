using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis
using SparseArrays
using Arpack
using StaticArrays
using NearestNeighbors

using Distances


"""
    # ABC Flow Field

"""

A = sqrt(3)
B = sqrt(2)
C = 1
function vec_field!(du,u,p,t)
    x,y,z = u
    du[1] = A*sin(z)+C*cos(y)
    du[2] = B*sin(x)+A*cos(z)
    du[3] = C*sin(y)+B*cos(x)
    # mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

res = 25
xyz = [[x y z] for x in LinRange(0,2*pi,res), y in LinRange(0,2*pi,res), z in LinRange(0,2*pi,res)]
X = vcat(xyz...)

t = (0:0.25:6)

@time traj = solve_ensemble(vec_field!,X,t)

# ------
metric = PeriodicEuclidean([2*pi,2*pi,2*pi])
# @time P = dynamic_laplacian(traj, metric=metric)
# ϵ = nndist(traj.X[:,1,:])
@time A, D = npg_ϵ(traj, metric=metric, ϵ=2*ϵ)

r = 7
# @time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?
@time λ, v = eigs(D-A, D, nev=r+1, which=:SM) 

r = 6
V = real(v[:,1:r])
S, R = SEBA(V)

# ------
# display(scatter(traj.X[1:10,:,1]', traj.X[1:10,:,3]'))

exportf(traj,S,"./examples/abc_flow/output/") # export

