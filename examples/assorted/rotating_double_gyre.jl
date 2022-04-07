# OLD SCRIPT; STORED JUST IN CASE

using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations
using JLD
using Profile
# include("CoherentStructures.jl")
# using .CoherentStructures

# -------------------------------------------------------------------------------------------------------
# define velocity field 
function s(t)
    if (t < 0)
        return 0
    elseif(t > 1)
        return 1
    else
        return (t.^2).*(3-2*t);
    end
end

function flow_field_double_gyre(x,y,t)
    # flow field for the double gyre; calculated as the curl of the double gyre vector stream function
    temp = s(t);
    return [pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)]*[1-temp; temp];
end

function vec_field(u,p,t)
    return flow_field_double_gyre(u[1],u[2],t)
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

# xy = [[x y] for x = 0:.025:1, y = 0:.025:1]
# X = vcat(xy...)
X_0 = create_grid([0 1],[0 1],[.025 .025])

# n_s = size(X,1)
# t = (0:0.01:1)
# tspan = (t[1],t[end])
# T = t[end]-t[1]
# n_t = length(t)
# # traj_data = calculate_trajectories(x,y,t,flow_field_double_gyre) # naive method
# # x_traj = traj_data[1]
# # y_traj = traj_data[2]

# # x_traj = zeros(n_s,n_t)
# # y_traj = zeros(n_s,n_t)
# X_traj = Array{Float64}(undef,n_s,n_t,2)
# for i = 1:n_s
#     local u0 = [X[i,1];X[i,2]]
#     local prob = ODEProblem(vec_field,u0,tspan)
#     local sol = solve(prob,saveat=t)
#     for j = 1:n_t
#         # x_traj[i,j] = sol.u[j][1]
#         # y_traj[i,j] = sol.u[j][2]
#         X_traj[i,j,:] = sol.u[j]
#     end
# end

t = (0:0.5:1)
traj = solve_trajectory(vec_field,X_0,t)

# -------------------------------------------------------------------------------------------------------
# construct dynamic Laplacian; using CoherentStructures.jl
# custom metric as an input
@time P = dynamic_laplacian(X_traj,t)

# Profile.print()
# Profile.clear()
# -------------------------------------------------------------------------------------------------------
# calculate eigens 
λ, v = eigen(P')

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing
r = 2
V = real(v[:,end:-1:end-r+1]) 
S, R = SEBA(V)

# -------------------------------------------------------------------------------------------------------
# plot
# scatter?
display(scatter(X_traj[:,1,1], X_traj[:,1,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
display(scatter(X_traj[:,1,1], X_traj[:,1,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))