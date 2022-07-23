using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations
using SparseArrays
using Arpack
using StaticArrays
using ProfileView
# using JLD
# using Profile

# -------------------------------------------------------------------------------------------------------
# define velocity field 
# function s(t)
#     if (t < 0)
#         return 0
#     elseif(t > 1)
#         return 1
#     else
#         return (t.^2).*(3-2*t);
#     end
# end

# function flow_field_double_gyre(x,y,t)
#     # flow field for the double gyre; calculated as the curl of the double gyre vector stream function
#     temp = (t.^2).*(3-2*t)
#     return [pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)]*[1-temp; temp];
# end

"""
    # Rotating double gyre flow field.
For 0 <= t <= 1.
"""
# function vec_field(u,t)
#     temp = (t.^2).*(3-2*t)
#     # don't think static matrix/vector has helped much here...
#     return vec(SMatrix{2,2}([pi*sin(2*pi*u[1])*cos(pi*u[2]) 2*pi*sin(pi*u[1])*cos(2*pi*u[2]); -2*pi*cos(2*pi*u[1])*sin(pi*u[2]) -pi*cos(pi*u[1])*sin(2*pi*u[2])])*SVector{2}([1-temp; temp]));
# end

function vec_field!(du,u,p,t)
    temp = (t.^2).*(3-2*t)
    temp = [pi*sin(2*pi*u[1])*cos(pi*u[2]) 2*pi*sin(pi*u[1])*cos(2*pi*u[2]); -2*pi*cos(2*pi*u[1])*sin(pi*u[2]) -pi*cos(pi*u[1])*sin(2*pi*u[2])]*[1-temp; temp]
    du = temp
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

X_0 = create_grid([0 1],[0 1],[60 60])
t = (0:0.1:1)

@time traj = solve_trajectory(vec_field,X_0,t)

# -------------------------------------------------------------------------------------------------------
# construct dynamic Laplacian
# custom metric as an input

# ProfileView.@profview
P = dynamic_laplacian(traj, threshold=0.0)

# -------------------------------------------------------------------------------------------------------
# calculate eigens 
r = 2

# create sparse, https://discourse.julialang.org/t/fast-calculation-of-eigenvalues/30359/5
# @time λ, v = eigen(P') # transpose?
@time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing

# V = real(v[:,end:-1:end-r+1]) 
V = real(v[:,1:2])
S, R = SEBA(V)

# export_eig(traj,S,"./examples/rotating_double_gyre/output/") # export

# -------------------------------------------------------------------------------------------------------
# # plot
display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))