using Revise
using DMAP
using Plots
using SparseArrays
using LinearAlgebra
using Arpack
using Clustering
using Random

using CSV
using DataFrames
using Distances
using NearestNeighbors


"""
    *** VON KARMAN VORTEX STREET *** 

Trajectories are calculated from openfoam solution files in Paraview.


openfoam case: ...

"""


# --- import trajectories ---

# - add boolean matrix to store active trajectories
# - pre-scan for number of particle ids

# inputs 
# input_dir = "./examples/von_karman/trajectories_200res_all/"
# input_dir = "./examples/von_karman/traj_thesis/"
input_dir = "./examples/von_karman/traj_final_02/"

get_every = 1
dimensions = 2
input_files = readdir(input_dir)[(1:get_every:end)] # (!) note that time slices may be out of order according to how this list is read

# reorder
n_steps = length(input_files)
df = DataFrame(t=Float64[], sub=DataFrame[])
for i = 1:n_steps
    df_i = CSV.read(input_dir * input_files[i], DataFrame)
    push!(df,(df_i."Time"[1],df_i))
end
sort!(df)
sub_t = (1:3:length(input_files))
df = df[sub_t,:]
n_steps = size(df,1)



df_vec = [df."sub"[i] for i = 1:n_steps] # ensure .DS_Store is not accidentaly in the trajectory folder
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


function is_in(val,vec)
    # (!) assume val only occurs once in vec, if at all
    ind = sum((val.==vec).*(1:length(vec)))
    if (ind>0)
        return true, ind
    end
    return false, ind
end


# clip DynTrajectory into a normal Trajectory
X = dtraj.X
n_t = length(dtraj.t)
X_out = Array{Float64}(undef,dtraj.n_pids,n_steps,2) # generalise dimensions (!)

# df[1,2]."Q Criterion"[keep_ids]

# determine trajectory ids that exists for the entire time interval
keep_ids = Set(dtraj.ids[1])
for t_i = 2:dtraj.n_steps
    intersect!(keep_ids,Set(dtraj.ids[t_i]))
end
keep_ids = collect(keep_ids)
n_s = length(keep_ids)

# construct traj.X
X = Array{Float64}(undef,n_s,dtraj.n_steps,2)  # generalise dimensions (!)
init_slice = Vector{Int64}(undef,n_s)
for t_i = 1:dtraj.n_steps
    for s_i = 1:n_s
        _,idx = is_in(keep_ids[s_i],dtraj.ids[t_i])
        X[s_i,t_i,:] = dtraj.X[t_i][idx,:]
        if (t_i==1)
            init_slice[s_i] = idx
        end
    end
end

traj = Trajectory(X,dtraj.t)

# sub_t = [(1:2:35);36]
# traj = Trajectory(traj.X[:,sub_t,:],traj.t[sub_t])

# --- see which trajectories stay in the field
# display(scatter(traj.X[:,1,1],traj.X[:,1,2])) 

# --- need histogram of nndists
# nndist(traj.X[:,end,:],mode=:plot) 

nn_avg,nn_list = sp_nndist(traj.X[:,1,:],mode=:plot)
# nn_avg_2,nn_list_2 = sp_nndist(traj.X[:,end,:],mode=:plot)
# histogram(nn_list,bins=20,xlims=(0.0,0.1),color=:black)
plot(xlabel="(2d+1)th Nearest Neighbour",ylabel="Number of points in each value group")
# histogram!(nn_list,bins=10,xlims=(0.0,1.0),color=:blue,label="t=1.00")
histogram!(nn_list_2,bins=100,xlims=(0.0,0.08),color=:black,label=:none)
# epsilon = 0.0350004090185647/sqrt(2) # von karman res 200x200 01


epsilon = 0.04/sqrt(2) # von karman res 200x200 02 (larger grid domain)


# # ------- edge detection!!!!
# # lambda = 4, k = 50
@time boundary_idxs, boundary_mask, perc = get_edge(traj.X[:,1,:],1.5,40); display(scatter(traj.X[boundary_idxs,1,1],traj.X[boundary_idxs,1,2],zcolor=boundary_mask*1.0,aspect_ratio=1,c=:grayC))
# display(scatter(traj.X[inds,1,1],traj.X[inds,1,2],aspect_ratio=1))
# # survey_edge_params(traj.X[:,1,:],(0:0.25:7),(2:10:101))



# # # --- explore coherence ---
# # @time P = dynamic_laplacian(traj, threshold=0.,epsilon=0.002/sqrt(2))
# # @time P = dynamic_laplacian(traj, threshold=0.)
@time P = sp_dynamic_laplacian(traj,epsilon=epsilon,k=30,boundary_idxs=boundary_idxs)


r = 40
@time λ, v = eigs(sparse(P'), nev=r, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?
display(scatter(real(λ),color=:black,xlabel="k",ylabel="Eigenvalue Value",legend=:none))

r=18
V = real(v[:,1:r])
S, R = SEBA(V[:,1:r])
inds = randperm(size(traj.X,1))[1:5000]
for i = 1:r
    display(scatter(traj.X[inds,end,1], traj.X[inds,end,2], zcolor=real(S[inds,i]), c=:jet, aspectratio=1))
end

# turn this into a movie?
# forwards
t1 = 5; t2 = 9
ftle_field_f,avg_f = ftle(traj,[t1;t2],0.05)
ftle_field_f .= (ftle_field_f.>=0).*ftle_field_f
display(scatter(traj.X[inds,t1,1], traj.X[inds,t1,2], zcolor=ftle_field_f[inds], c=:magma, aspectratio=1))

# backwards
t1 = 5; t2 = 1
ftle_field_b,avg_b = ftle(traj,[t1;t2],0.05)
ftle_field_b .= (ftle_field_b.>=0).*ftle_field_b
display(scatter(traj.X[inds,t1,1], traj.X[inds,t1,2], zcolor=ftle_field_b[inds], c=:magma, aspectratio=1))

# # # R = kmeans(V[:,1:20]', 20)
# # # a = assignments(R)

# # display(scatter(real(λ)))
# # # display(scatter(traj.X[:,end,1], traj.X[:,end,2], zcolor=real(S[:,1]), c=:viridis, aspectratio=1))
# # # display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(S[:,2]), c=:viridis, aspectratio=1))

# display(scatter(traj.X[:,end,1], traj.X[:,end,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))


# display(scatter(traj.X[:,end,1], traj.X[:,end,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))
# # display(scatter(traj.X[:,end,1], traj.X[:,end,2], zcolor=a, c=palette(:tab20), aspectratio=1))
# # display(scatter( real(V[:,2]), real(V[:,3]), real(V[:,6]), zcolor=real(S[:,1]), camera = (45,45) ))

# --- q-criterion initial
Q_0 = df[1,2]."Q Criterion"[init_slice]
Q_0 = (Q_0.>10)*1.0
display(scatter(traj.X[inds,1,1], traj.X[inds,1,2], zcolor=Q_0[inds], c=:viridis, aspectratio=1))

# partition id vec
features = 6
mask = [(i <= features ? 1 : 0) for i = 1:r]
partition = [ sum(((1:r).*mask).*(S[i,:].>0.5))
    for i = 1:size(traj.X,1)]
display(scatter(traj.X[inds,1,1], traj.X[inds,1,2], zcolor=partition[inds], aspectratio=1))


# exportf(traj.X[:,end,:],S,"/Users/jackh/out.csv")

# exportf(traj.X[:,1,:],S,"/Users/jackh/out_init.csv")
exportf(traj,[boundary_mask Q_0 S ftle_field_f ftle_field_b partition],"./examples/von_karman/output_02/")

# # r = 10
# # @time λ, v = eigs(D-A, D, nev=r, which=:SM)

# # csets = 5
# # V = real(v[:,1:10])
# # S, R = SEBA(V[:,1:2])

# # X_0 = dtraj.X[1]
# # display(scatter(real(λ)))
# # display(scatter(X_0[:,1], X_0[:,2],zcolor=real(S[:,1])))

# # exportf(dtraj.X[1],[V S],"/Users/jackh/out.csv")


