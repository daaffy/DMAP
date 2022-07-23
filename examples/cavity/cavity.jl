using Revise
using DMAP
using Plots
using SparseArrays
using LinearAlgebra
using Arpack

using CSV
using DataFrames
using Distances
using NearestNeighbors


"""
    *** LID-DRIVEN CAVITY FLOW *** 

Trajectories are calculated from openfoam solution files in Paraview.


openfoam case: https://www.openfoam.com/documentation/tutorial-guide/2-incompressible-flow/2.1-lid-driven-cavity-flow#x6-60002.1

"""


# --- import trajectories ---

# - add boolean matrix to store active trajectories
# - pre-scan for number of particle ids

# inputs 
input_dir = "./examples/cavity/trajectories_02/"

dimensions = 2
input_files = readdir(input_dir) # (!) note that time slices may be out of order according to how this list is read

n_steps = length(input_files)
df_vec = [CSV.read(input_dir * input_files[i], DataFrame) for i = 1:n_steps]
t = [df_vec[i][1,"Time"] for i = 1:n_steps]
n_pids = maximum([maximum(df_vec[i][:,"ParticleId"]) for i = 1:n_steps])+1
ids = [df_vec[i][:,"ParticleId"] for i = 1:n_steps]
X = [Matrix(df_vec[i][:,["Points:0","Points:1"]]) for i = 1:n_steps]


# define DynTrajectory
struct DynTrajectory
    X::Vector
    ids::Vector
    t::Vector{Float64}
    n_pids::Int
    n_steps::Int
end

dtraj = DynTrajectory(X,ids,t,n_pids,n_steps)


# --- explore coherence ---


# * npg_k *
# function npg_k(dtraj::DynTrajectory)
    k = 5
    metric = Euclidean()

    n_s = dtraj.n_pids
    n_t = dtraj.n_steps

    # idxs_1 = similar(idxs_2)
    # idxs = Vector{Vector{Int64}}(undef,n_s)

    A = spzeros(n_s,n_s)
    for t_i = 1:n_t
    # t_i = 1
        X_i = dtraj.X[t_i]
        tree = KDTree(X_i',metric)
        idxs_2 = vcat([(id+1)*ones(Int32,k) for id in dtraj.ids[t_i]]...)
        idxs_loc = knn(tree,X_i',k)[1] # .=    MethodError: no method matching copyto!(::Tuple{Vector{Vector{Int64}}...
        idxs_1_loc = vcat(idxs_loc...)
        idxs_1 = [dtraj.ids[t_i][j]+1 for j in idxs_1_loc]
        
        # A += sparse(idxs_1,idxs_2,1,n_s,n_s)
        # sparse()
        global A += sparse(idxs_1,idxs_2,1,n_s,n_s) + sparse(idxs_2,idxs_1,1,n_s,n_s)
    end
    A += A'
    A = 1*(A.>0)
    A -= sparse(I,n_s,n_s)
    D = sparse((1:n_s),(1:n_s),vec(sum(A,dims=1)))

    # return A, D
# end

# # A, D = npg_k(dtraj)

r = 10
@time λ, v = eigs(D-A, D, nev=r, which=:SM)

csets = 5
V = real(v[:,1:10])
S, R = SEBA(V[:,1:2])

X_0 = dtraj.X[1]
display(scatter(real(λ)))
display(scatter(X_0[:,1], X_0[:,2],zcolor=real(S[:,1])))

exportf(dtraj.X[1],[V S],"/Users/jackh/out.csv")

# --- plot trajectories ---

# function inn(vec,val)
#     for i = 1:length(vec)
#         if vec[i] == val
#             return true, i
#         end
#     end
#     return false, 0
# end

# function plot_dtraj(dtraj,id)
#     # assume once particles leave the domain they don't come back...
#     pX = Vector{Vector}(undef,0) # ::Any[]
#     for i = 1:dtraj.n_steps
#         is_in, idx = inn(dtraj.ids[i],id)
#         if (is_in)
#             push!(pX,dtraj.X[i][idx,:])
#         end
#     end
#     pX = hcat(pX...)

#     display(scatter!(pX[1,:], pX[2,:]))
    
#     # return pX
# end

# scatter()
# plot_dtraj(dtraj,1000)
# plot_dtraj(dtraj,260)
# plot_dtraj(dtraj,700)
# plot_dtraj(dtraj,555)
# plot_dtraj(dtraj,611)

