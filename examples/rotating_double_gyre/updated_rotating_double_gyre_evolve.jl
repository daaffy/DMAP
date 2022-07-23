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

"""
    # Rotating double gyre flow field.
For 0 <= t <= 1.
"""
seconds = 5
fps = 5
t_query = LinRange(0.01,4,seconds*fps)

s = t -> (t^2)*(3-2*t)
# g = t -> s(t_query-t)
function vec_field!(du,u,p,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

res = 400
X_0 = create_grid([0 1],[0 1],[res res])
n_s = size(X_0,1)

function def_circle(centre,radius,id)
    return x -> ( ((x[1]-centre[1])^2 + (x[2]-centre[2])^2 <= radius^2) ? id : 0 )
end
set_1 = def_circle([0.25,0.25],0.2,1)
set_2 = def_circle([0.45,0.8],0.15,2)
set_3 = def_circle([0.75,0.5],0.1,4)

# set_1 = x -> (x[1] >= 0.5 ? 1 : 0) # half/half

anim = @animate for i ∈ 1:length(t_query)
    g = t -> s(t_query[i]-t)
    t = [0,t_query[i]]
    @time traj = solve_ensemble(vec_field!,X_0,t,g)
    vals = [set_1(traj.X[i,2,:]) + 
            set_2(traj.X[i,2,:]) + 
            set_3(traj.X[i,2,:]) 
            for i = 1:n_s]

    # vals = [set_1(traj.X[i,2,:])
    #         for i = 1:n_s]

    display(heatmap(reshape(vals,res,res)',c=:bone,aspectratio=1,axis=nothing,foreground_color_subplot=colorant"white",legend=false))
    # :linear_bmy_10_95_c71_n256
end
gif(anim, "anim_fps15.gif", fps = fps)



# t = (0:0.1:1)





# vals = Vector{Int32}(undef,size(traj.X,1))






# vals = [ traj.X[i,2,1]
#     for i = 1:n_s
# ]







# @time P = dynamic_laplacian(traj, threshold=0.)

# # -------------------------------------------------------------------------------------------------------
# # calculate eigens 
# r = 4

# # create sparse, https://discourse.julialang.org/t/fast-calculation-of-eigenvalues/30359/5
# # @time λ, v = eigen(P') # transpose?

# @time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?
# # λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # NEED TO GET RID OF TRANSPOSE

# # for npg...
# # @time λ, v = eigs(D-A, D, nev=r+1, which=:SM) # sparse method; transpose converts sparse matrix back into a dense matrix?
# # v = sqrt.(D)*v

# # -------------------------------------------------------------------------------------------------------
# # SEBA post-processing

# # V = real(v[:,end:-1:end-r+1]) 
# V = real(v[:,1:3])
# S, R = SEBA(V)
# # S = V

# # export_eig(traj,S,"./examples/rotating_double_gyre/output/") # export

# # -------------------------------------------------------------------------------------------------------
# # # plot
# # display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))
# # display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))
# display(heatmap(reshape(real(S[:,1]),res,res)',c=:grays))
# display(heatmap(reshape(real(S[:,2]),res,res)',c=:grays))

# # # -------------------------------------------------------------------------------------------------------
# # # transport illustration

# # # function multiply(A,dim=1)
# # #     return []
# # # end

# # # cs_vecs = S
# # # thresh = 0.2
# # # # vals = sum((cs_vecs.>thresh).*(1:size(cs_vecs,2))',dims=2)
# # # # vals = X_0[:,1].*(cs_vecs[:,1].<thresh)
# # # gr(size=(1000,1000),legend=false,markerstrokewidth=0,markersize=20)
# # # # display(scatter(traj.X[:,10,1],traj.X[:,10,2], zcolor=vals,legend=false, markerstrokewidth=0,markershape=:circle, c=palette(:grayC)))


# # # inds = Vector{Int64}(undef,0)
# # # vals = Vector{Float64}(undef,0)
# # # for i = 1:80*80
# # #     if (X_0[i,1]<=0.5)
# # #        append!(inds,i)
# # #        append!(vals,cs_vecs[i,1]>thresh)
# # #     end
# # # end
# # # display(scatter(traj.X[inds,10,1],traj.X[inds,10,2], zcolor=vals,legend=false, markerstrokewidth=0,markershape=:circle, c=palette(:jet)))
# # # # inds = [(cs_vecs[i,q]>thresh && X_0[i,1]<=0.5) for i = 1:80]

# # # # display(heatmap(reshape(vals,80,80)))