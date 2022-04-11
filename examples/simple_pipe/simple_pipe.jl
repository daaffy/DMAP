using Revise
using DMAP
using LinearAlgebra
using Plots
using CSV
using DataFrames
using DifferentialEquations

# -------------------------------------------------------------------------------------------------------
# extract discrete velocity field to X, dX multi-dimensional arrays.

input_dir = "./examples/simple_pipe/input/"
input_files = readdir(input_dir)

dvf = load_dvf(input_files,input_dir)

# -------------------------------------------------------------------------------------------------------
# build interpolated velocity field

@time cvf = uinterp(dvf,"shep")

# -------------------------------------------------------------------------------------------------------
# integrate trajectories

X_0 = dvf.X
t = LinRange(0,0.01,20)
@time traj = solve_trajectory(cvf,X_0,t)

indicator(x) = x[2] < 0.04 && x[2] > 0.01 ? true : false
traj = clip(traj,indicator) # throw away trajectories that leave the domain of interest; defined by indicator()

# -------------------------------------------------------------------------------------------------------
# process coherent structures

@time P = dynamic_laplacian(traj)
@time λ, v = eigen(P') # transpose?

r = 2
V = real(v[:,end:-1:end-r+1]) 
@time S, R = SEBA(V)

export_eig(traj,S,"./examples/simple_pipe/output/") # export

# -------------------------------------------------------------------------------------------------------
# test plots

# display(scatter(traj.X[:,:,1]',traj.X[:,:,2]'))
display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))

# -------------------------------------------------------------------------------------------------------
# junk

# # test
# i = 1
# j = 7
# eval_points = dvf.X[i,:]
# println(dvf.dX[i,:,j])
# println(cvf(eval_points,dvf.t[j]))
