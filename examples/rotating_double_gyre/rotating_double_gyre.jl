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

X_0 = create_grid([0 1],[0 1],[.025 .025])
t = (0:0.5:1)
traj = solve_trajectory(vec_field,X_0,t)

# -------------------------------------------------------------------------------------------------------
# construct dynamic Laplacian
# custom metric as an input

@time P = dynamic_laplacian(traj)

# -------------------------------------------------------------------------------------------------------
# calculate eigens 
λ, v = eigen(P') # transpose?

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing
r = 2
V = real(v[:,end:-1:end-r+1]) 
S, R = SEBA(V)

# -------------------------------------------------------------------------------------------------------
# plot
display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))