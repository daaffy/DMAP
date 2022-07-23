using Revise
using DMAP
using LinearAlgebra
using Plots
using Profile
using SparseArrays
using Arpack
using NearestNeighbors

using BenchmarkTools

# -------------------------------------------------------------------------------------------------------
# extract discrete velocity field to X, dX multi-dimensional arrays.


### STENOSIS
input_dir = "./examples/stenosed_pipe/input/"
input_files = readdir(input_dir)
t = LinRange(0,0.7,70)
dvf = load_dvf(input_files,input_dir,t=t) # x-axis: -0.016 to 0.059
dvf = DiscreteVelocityField(dvf.X,dvf.dX/10,dvf.t)
X_0 = dvf.X
X_0 = clip(X_0,x->.01<x[1]<.035)

# ### SIMPLER PIPE
# input_dir = "./examples/simpler_pipe/input/"
# input_files = readdir(input_dir)
# t = [0;0.5;1]
# dvf = load_dvf(input_files,input_dir,t=t) # x-axis: -0.016 to 0.059
# dvf = DiscreteVelocityField(dvf.X,dvf.dX,dvf.t)
# X_0 = dvf.X
# # X_0 = clip(X_0,x->.01<x[1]<.035)
# X_0 = sift(X_0,0.1,0)

# #### 2D circle
# X_0 = vcat([[r*cos(θ) r*sin(θ)] for r = LinRange(0,1,30)[2:end], θ = LinRange(0,2*pi,60)[2:end]]...)
# scatter(X_0[:,1],X_0[:,2])

# -------------------------------------------------------------------------------------------------------


# i = 200 # query point

# lambda = 4
# kd = KDTree(X_0')
# k = Integer(floor(size(X_0,1)*0.025))
# # idxs, _ = knn(kd,X_0'[:,i],k,sortres=true)
# idxs, dists = knn(kd,X_0',k,true)
# Z = hcat(dists...)[2,:] # vector of distance to nearest neighbour from ith point
# # X_0[hcat(idxs...),:]
# C_i = sum(X_0[hcat(idxs...),:],dims=1)[1,:,:]/k
# del = colwise(Euclidean(),C_i',X_0')
# boundary_mask = del.>lambda*Z
# println(sum(boundary_mask))
# # C_i = sum(X_0'[:,idxs],dims=2)# centroid

boundary_idxs, boundary_mask, perc = get_edge(X_0,4,80) # for stenosis
# boundary_idxs, boundary_mask, perc = get_edge(X_0,2,40) # for simple cylinder
# boundary_idxs, boundary_mask, perc = get_edge(X_0,1,10) # for 2d disc
print(perc)
# boundary_idxs, boundary_mask, num = get_edge(X_0,12.5,Integer(floor(size(X_0,1)*0.5))) 

# print(num/size(X_0,1))
# exportf(X_0,boundary_mask,pwd()*"/examples/other/out.csv")

###########trajectory calculation
# cvf = uinterp(dvf,"shep")
# t = collect(LinRange(0,0.7,100))
# @time traj = solve_trajectory(cvf,X_

# exportf(traj,boundary_mask,pwd()*"/examples/other/boundary_advection")


# print(C_i)
# print(X_0'[:,i])
# X_0'[:,idxs]
# exportf()

############
temp = [X_0[i,k] for i = 1:size(X_0,1), j = 1:1, k = 1:size(X_0,2)]
traj = Trajectory(temp,[0])

# P = dynamic_laplacian(traj,boundary_idxs=boundary_idxs,eps_mult=1.)
P = dynamic_laplacian(traj,boundary_idxs=boundary_idxs,epsilon=exp(-12.5)) # for stenosis example, based on linear region of findϵ 
# P = dynamic_laplacian(traj,boundary_idxs=boundary_idxs,epsilon=nndist(X_0)/sqrt(2))

r = 6
λ, v = eigs(sparse(P'), nev=r+1, which=:LM)

exportf(X_0,[real(v) boundary_mask],"./examples/other/boundary_test_2.csv")

# for x in [1 2 3]
#     print(x)
# end


# -----
# @time cvf = uinterp(dvf,"shep") # the interpolation function is extremely heavy
# # cvf = (x,t)->[1;0;0]

# # X_0 = sift(dvf.X,0.2,0)

# # indicator(x) = 0.015 < x[1] < 0.03 ? true : false
# X_0 = sift(dvf.X,0.2,0)
# # X_0 = clip(X_0,indicator)


# t = collect(LinRange(0,0.7,100))
# @time traj = solve_trajectory(cvf,X_0,t)
# # ProfileView.@profview 

# # clip here

# P = dynamic_laplacian(traj, threshold=0.8)

# # # --- DVF integrator debug
# #     g(x,p,t) = [1,0,0]
# #     t = collect(0:0.1:1)
# #     u0 = [0 0 0;1 0 0]'
# #     prob = ODEProblem(g,u0,(0,1))
# #     sol = solve(prob,saveat=t)
# # # ---

# r = 6
# λ, v = eigs(sparse(P'), nev=r+1, which=:LM)

# V = real(v[:,1:r])
# S, R = SEBA(V)

# export_eig(traj,S,"./examples/stenosed_pipe/output/") # export





