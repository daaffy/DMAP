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
function vec_field!(du,u,s,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y);
             -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],
            [1-s(t); s(t)])
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

res = 200
X_0 = create_grid([0 1],[0 1],[res res])
t = (0:0.1:1)

@time traj = solve_ensemble(vec_field!,X_0,t,s)


boundary_idxs, boundary_mask, perc = get_edge(X_0,0.5,10) # lambda = 0.5, k = 0.5

inds = randperm(size(traj.X,1))[1:5000]
display(scatter(X_0[inds,1],X_0[inds,2],zcolor=boundary_mask[inds]*1.0,aspect_ratio=1))
display(scatter(traj.X[inds,end,1],traj.X[inds,end,2],zcolor=traj.X[inds,1,1],aspect_ratio=1))

# survey_edge_params(X_0,(0:0.1:2),(2:1:20))


# eps calculation
# nn_avg,nn_list = sp_nndist(traj.X[:,1,:],mode=:plot)
nn_avg_2,nn_list_2 = sp_nndist(traj.X[:,end,:],mode=:plot)
plot(xlabel="(2d+1)th Nearest Neighbour",ylabel="Number of points in each value group")
# histogram!(nn_list,bins=10,xlims=(0.0,1.0),color=:blue,label="t=1.00")
histogram!(nn_list_2,bins=2000,xlims=(0.0,0.02),color=:black,label=:none)

epsilon = 0.008/sqrt(2)

# -------------------------------------------------------------------------------------------------------
# construct dynamic Laplacian
# custom metric as an input

# ProfileView.@profview

# @time P = dynamic_laplacian(traj, threshold=0.)
# Revise.retry(); 
@time P = sp_dynamic_laplacian(traj,k=5,epsilon=epsilon,boundary_idxs=boundary_idxs)

# @time A, D = npg_k(traj, k=4)

# tree = KDTree(traj.X[:,1,:]')
# idxs = knn(tree,traj.X[1,1,:],2)

# tree = BallTree(traj.X[:,1,:]')
# idxs = inrrange(tree,traj.X[:,1,:],)

# -------------------------------------------------------------------------------------------------------
# calculate eigens 
r = 30

# create sparse, https://discourse.julialang.org/t/fast-calculation-of-eigenvalues/30359/5
# @time λ, v = eigen(P') # transpose?

@time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?
display(scatter(real(λ),color=:black,xlabel="k",ylabel="Eigenvalue Value",legend=:none))

# λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # NEED TO GET RID OF TRANSPOSE

# for npg...
# @time λ, v = eigs(D-A, D, nev=r+1, which=:SM) # sparse method; transpose converts sparse matrix back into a dense matrix?
# v = sqrt.(D)*v

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing

# V = real(v[:,end:-1:end-r+1]) 
V = real(v[:,1:4])
S, R = SEBA(V[:,1:2])
# S = V

# export_eig(traj,S,"./examples/rotating_double_gyre/output/") # export

# -------------------------------------------------------------------------------------------------------
# # plot
display(scatter(traj.X[:,10,1], traj.X[:,10,2], zcolor=real(S[:,1]), c=:magma, aspectratio=1))
display(scatter(traj.X[:,10,1], traj.X[:,10,2], zcolor=real(S[:,2]), c=:magma, aspectratio=1))

# diffusion map embedding
# display(scatter( real(V[:,2]), real(V[:,3]), real(V[:,4]), zcolor=real(S[:,1]), camera = (45,45) ))

# exportf(traj.X[:,1,:],[V S],"/Users/jackh/out.csv")
# display(heatmap(reshape(real(S[:,1]),res,res)',c=:grays))
# display(heatmap(reshape(real(S[:,2]),res,res)',c=:grays))


display(scatter(traj.X[:,end,1], traj.X[:,end,2], aspectratio=1,label=:none,c=:black,xlims=(-0.1,1.1)))




inds = (1:size(traj.X,1))
# forwards
t1 = 5; t2 = 9
ftle_field_f,avg_f = ftle(traj,[t1;t2],0.02)
ftle_field_f .= (ftle_field_f.>=0).*ftle_field_f
display(scatter(traj.X[inds,t1,1], traj.X[inds,t1,2], zcolor=ftle_field_f[inds], c=:magma, aspectratio=1))

# backwards
t1 = 5; t2 = 1
ftle_field_b,avg_b = ftle(traj,[t1;t2],0.02)
ftle_field_b .= (ftle_field_b.>=0).*ftle_field_b
display(scatter(traj.X[inds,t1,1], traj.X[inds,t1,2], zcolor=ftle_field_b[inds], c=:magma, aspectratio=1))

features = 2
r = 2
mask = [(i <= features ? 1 : 0) for i = 1:r]
partition = [ sum(((1:r).*mask).*(S[i,:].>0.5))
    for i = 1:size(traj.X,1)]
display(scatter(traj.X[inds,end,1], traj.X[inds,end,2], zcolor=partition[inds], aspectratio=1))

# OKUBO WEISS
s_1 = (x,y)->4*pi^2*cos(2*pi*x)*cos(pi*y)
s_2 = (x,y)->3*pi^2*sin(2*pi*x)*sin(pi*y)*(1-s(t_val))-3*pi^2*sin(2*pi*y)*sin(pi*x)*s(t_val)
om = (x,y)->5*pi^2*sin(2*pi*x)*sin(pi*y)*(1-s(t_val))-5*pi^2*sin(2*pi*y)*sin(pi*x)*s(t_val)
f = (x,y)-> s_1(x,y)^2 + s_2(x,y)^2 - om(x,y)^2
# set_1 = x -> (f(x[1],x[2]) < 0) & (x[1] <= 0.5)
# set_2 = x -> (f(x[1],x[2]) < 0) & (x[1] > 0.5)
set = x -> (f(x[1],x[2]) < 0)
Q_0 = [set(traj.X[i,1,:]) for i = 1:size(traj.X,1)]
display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=Q_0, aspectratio=1))


exportf(traj,[boundary_mask Q_0 S ftle_field_f ftle_field_b partition real(v[:,1:4])],"./examples/rotating_double_gyre/output_final/")


# # -------------------------------------------------------------------------------------------------------
# # transport illustration

# # function multiply(A,dim=1)
# #     return []
# # end

# # cs_vecs = S
# # thresh = 0.2
# # # vals = sum((cs_vecs.>thresh).*(1:size(cs_vecs,2))',dims=2)
# # # vals = X_0[:,1].*(cs_vecs[:,1].<thresh)
# # gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=20)
# # # display(scatter(traj.X[:,10,1],traj.X[:,10,2], zcolor=vals,legend=false, markerstrokewidth=0,markershape=:circle, c=palette(:grayC)))


# # inds = Vector{Int64}(undef,0)
# # vals = Vector{Float64}(undef,0)
# # for i = 1:80*80
# #     if (X_0[i,1]<=0.5)
# #        append!(inds,i)
# #        append!(vals,cs_vecs[i,1]>thresh)
# #     end
# # end
# # display(scatter(traj.X[inds,10,1],traj.X[inds,10,2], zcolor=vals,legend=false, markerstrokewidth=0,markershape=:circle, c=palette(:jet)))
# # # inds = [(cs_vecs[i,q]>thresh && X_0[i,1]<=0.5) for i = 1:80]

# # # display(heatmap(reshape(vals,80,80)))