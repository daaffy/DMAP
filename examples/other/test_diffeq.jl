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

# =tt->(tt.^2).*(3-2*tt)

s = t -> (t^2)*(3-2*t)
function vec_field!(du,u,p,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

X_0 = create_grid([0 1],[0 1],[60 60])
t = (0:0.1:1)

@time traj = solve_ensemble(vec_field!,X_0,t,s)

###### TEST HERE ###########

# vec_field_temp(x,p,t) = vec(vec_field(x,t)) # form for ODE package # note vec() to turn from row vector to column vector (required form)


# function _solve(vec_field!,X_0,t)

    # n_t = length(t)
    # n_s = size(X_0,1)
    # d = size(X_0,2)

    # X = Array{Float64}(undef,n_s,n_t,d) # change this # not sure if I like the representation order; (n_s,d,n_t) more consistent?
    # for i = 1:n_s
    #     local u0 = vec(X_0[i,:])
    #     local prob = ODEProblem(vec_field!,u0,(t[1],t[end]),s)
    #     local sol = solve(prob,reltol=1e-6,saveat=t)
    #     for j = 1:n_t
    #         X[i,j,:] = sol.u[j]
    #     end
    # end

    # function _solve(vec_field!,X_0,t)    
    #     prob = ODEProblem(vec_field!,[0;0],(t[1],t[end]),s)
    #     initial_conditions = X_0'
    #     function prob_func(prob,i,repeat)
    #         remake(prob,u0=initial_conditions[:,i])
    #     end
    #     ensemble_prob = EnsembleProblem(prob; prob_func)
    #     sim = solve(ensemble_prob,Tsit5(),EnsembleSerial(),trajectories=size(initial_conditions,2),saveat=t)

    #     # construct X
    #     n_t = length(t)
    #     n_s = size(X_0,1)
    #     d = size(X_0,2)
    #     X = Array{Float64}(undef,n_s,n_t,d)
    #     for ti = 1:length(t)
    #         X[:,ti,:] = hcat(componentwise_vectors_timestep(sim,ti)...)
    #     end

    #     return Trajectory(X,t)
    # end



# X = [hcat(componentwise_vectors_timestep(sim,i)...) for i = 1:length(t)]
    # return Trajectory(X,t)
    # traj = Trajectory(X,t)
# end

# @time traj = _solve(vec_field!,X_0,t)

# @benchmark traj = _solve(vec_field!,reshape(X_0[100,:],:,2),t)
# print(traj.X)

# @benchmark traj = solve_trajectory(vec_field,reshape(X_0[100,:],:,2),t)
# @time traj = solve_trajectory(vec_field,X_0,t)

# -------------------------------------------------------------------------------------------------------
# construct dynamic Laplacian
# custom metric as an input

# ProfileView.@profview
P = dynamic_laplacian(traj, threshold=0.9)

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
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))