using Revise
using DMAP
using LinearAlgebra
using Plots
using Profile
using SparseArrays
using Arpack

# -------------------------------------------------------------------------------------------------------
# extract discrete velocity field to X, dX multi-dimensional arrays.

input_dir = "./examples/stenosed_pipe/input/"
input_files = readdir(input_dir)

# t = collect(0:0.01:0.01)
t = LinRange(0,0.7,70)
dvf = load_dvf(input_files,input_dir,t=t) # x-axis: -0.016 to 0.059
dvf = DiscreteVelocityField(dvf.X,dvf.dX/10,dvf.t)

@time cvf = uinterp(dvf,"shep") # the interpolation function is extremely heavy
# cvf = (x,t)->[1;0;0]

# X_0 = sift(dvf.X,0.2,0)

# indicator(x) = 0.015 < x[1] < 0.03 ? true : false
X_0 = sift(dvf.X,0.2,0)
# X_0 = clip(X_0,indicator)


t = collect(LinRange(0,0.7,100))
@time traj = solve_trajectory(cvf,X_0,t)
# ProfileView.@profview 

# clip here

P = dynamic_laplacian(traj, threshold=0.8)

# # --- DVF integrator debug
#     g(x,p,t) = [1,0,0]
#     t = collect(0:0.1:1)
#     u0 = [0 0 0;1 0 0]'
#     prob = ODEProblem(g,u0,(0,1))
#     sol = solve(prob,saveat=t)
# # ---

r = 6
λ, v = eigs(sparse(P'), nev=r+1, which=:LM)

V = real(v[:,1:r])
S, R = SEBA(V)

export_eig(traj,S,"./examples/stenosed_pipe/output/") # export






# temp = sift(dvf.X,0.1,0)
# X_0 = Array{Float64}(undef, size(temp,1), 1, 3)
# X_0[:,1,:] = temp
# traj = Trajectory(X_0,[0])

# export_traj(traj,"examples/stenosed_pipe/output")




# ---
# temp = Array{Float64}(undef,size(dvf.X,1),1,size(dvf.X,2))
# temp[:,1,:] = dvf.X
# traj = Trajectory(temp,vec([0]))

# x = range(-3.0, 3.0; length=100)
# K = kernelmatrix(SqExponentialKernel(), x)
# plot(
#     heatmap.([K...]; yflip=true, colorbar=false)
# )

# @time nndist(traj.X[:,1,:],mode=:plot)
# @time dist_matrix = pairwise(metric,traj.X[:,1,:],traj.X[:,1,:],dims=1)

# single_dl(traj)