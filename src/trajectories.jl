using DifferentialEquations, DifferentialEquations.EnsembleAnalysis, Random, ProgressMeter, NearestNeighbors # for boundary detection


using Plots # for survey_edge_params

export create_grid, solve_trajectory, test, clip, clip!, solve_ensemble
export sift, get_edge, survey_edge_params

struct ContinuousVelocityField
    # to be fed into solve_trajectory
end

# to-do:
# - creating grids in more complex geometries e.g., in simple pipe
# - cut-off points that leave the geometry
# ** use an indicator function of some kind to solve the above problems?

# function indicator()

# end

#  AbstractRange might be a more natural input
function create_grid(args...)
    # create a uniform grid in the appropriate representation to be fed into the trajectory solver
    
    # checks 
    # ...
    
    d = length(args)-1 # spatial dimensions
    # lims = [args[d_i][1]:args[d+1][d_i]:args[d_i][2] for d_i = 1:d]
    lims = [LinRange(args[d_i][1],args[d_i][2],args[d+1][d_i]) for d_i = 1:d]
    temp = [collect(i)' for i in Iterators.product(lims...)] # the transpose (') is the only way I could figure out to convert to an appropriate type
    return vcat(temp...)
end


"""
        sift

    Takes a collection of points as an input and returns a random sample of those points.
    Useful for defining initial points of a trajectory in complicated geometries (take dvf.X as input).

    Inspired by: https://discourse.julialang.org/t/sampling-without-replacement/1073/5

"""
function sift(
    X::Array{<:Real,2},
    n::Union{Int64,Float64},
    seed = 1 # supply the rng seed for reproducible samples
)   
    @assert seed >= 0

    if (isa(n,Float64))
        @assert 0 <= n <= 1 # express as a percentage
        n = Int(ceil(n*size(X,1)))
    end

    indices = collect(1:size(X,1))

    rng = MersenneTwister(seed)
    sub_indices = Vector{Int64}(undef,n)
    for i = 1:n
        ind = rand(rng,indices)
        sub_indices[i] = ind
        deleteat!(indices,findall(x->x==ind,indices))
    end

    return X[sub_indices,:]
end

"""

Create a function that sifts while preserving point density??

"""

# ...


"""
    get_edge
Boundary detection based on mean shift methods.

'lambda' - the classifcation threshold; requires calibration.
'k'      - use k nearest neighbours.
"""

function survey_edge_params(X,lambda_range=[],k_range=[])
    kd = KDTree(X')
    perc = Matrix{Float64}(undef,length(lambda_range),length(k_range))
    for i = 1:length(lambda_range)
        for j = 1:length(k_range)
            @assert k_range[j] > 1
            idxs, dists = knn(kd,X',k_range[j],true)
            Z = hcat(dists...)[2,:] # vector of distance to nearest neighbour from ith point
            C_i = sum(X[hcat(idxs...),:],dims=1)[1,:,:]/k_range[j] # centroids
            del = colwise(Euclidean(),C_i',X') # distance between centroid and query point
            boundary_mask = del.>lambda_range[i]*Z # classify
            perc[i,j] = sum(boundary_mask)/size(X,1)
            # _,_,perc[i,j] = get_edge(X,lambda_range[i],k_range[j])
        end
    end
    heatmap(k_range,lambda_range,log10.(perc))
end

function get_edge(
    X::Array{<:Real},
    lambda::Real,
    k::Int64
)   
    @assert k > 1
    kd = KDTree(X')
    idxs, dists = knn(kd,X',k,true)
    Z = hcat(dists...)[2,:] # vector of distance to nearest neighbour from ith point
    C_i = sum(X[hcat(idxs...),:],dims=1)[1,:,:]/k # centroids
    del = colwise(Euclidean(),C_i',X') # distance between centroid and query point
    boundary_mask = del.>lambda*Z # classify

    temp = boundary_mask.*(1:size(X,1))
    boundary_idxs = deleteat!(temp, temp.==0)
    return boundary_idxs, boundary_mask, sum(boundary_mask)/size(X,1)
end

# to-do??
# - vec_field::DiscreteVectorField solver for large datasets (pre-interpolating the whole discrete vector field is too expensive)

# function solve_trajectory_ensemble(
#     vec_field::Function, # make this an abstract type
#     X_0::Matrix{<:Real},
#     t::Vector{<:Real}
# )
#     X_0 = create_grid([0 1],[0 1],[60 60])
#     t = (0:0.1:1)

#     n = size(X_0,1)
#     prob = ODEProblem(vec_field,[0.0;0.0],(t[1],t[end]))
#     function prob_func(prob,i,repeat)
#         remake(prob,u0=vec(X_0[i,:]))
#     end
#     ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
#     @time sim = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=n,solveat=t)
# end

"""
    solve_ensemble
Utilises DifferentialEquations.jl's ensemble solving capabilities to calculate
the trajectories.

"""
function solve_ensemble(vec_field!,X_0,t,s=:no_param)
    # note this doesnt work yet if you only input a single point.. should be an easy fix
    
    if (s==:no_param)
        prob = ODEProblem(vec_field!,X_0[1,:],(t[1],t[end]))
    else
        prob = ODEProblem(vec_field!,X_0[1,:],(t[1],t[end]),s)
    end

    initial_conditions = X_0'
    function prob_func(prob,i,repeat)
        remake(prob,u0=initial_conditions[:,i])
    end
    ensemble_prob = EnsembleProblem(prob; prob_func)
    sim = solve(ensemble_prob,Tsit5(),EnsembleSerial(),trajectories=size(initial_conditions,2),saveat=t)

    # construct X
    n_t = length(t)
    n_s = size(X_0,1)
    d = size(X_0,2)
    X = Array{Float64}(undef,n_s,n_t,d)
    for ti = 1:length(t)
        X[:,ti,:] = hcat(componentwise_vectors_timestep(sim,ti)...)
    end

    return Trajectory(X,t)
end

"""
    solve_ensemble_dvf

"""

"""
    solve_trajectory
Replaced by 'solve_ensemble'.

"""
function solve_trajectory(
    vec_field::Function, # make this an abstract type
    X_0::Matrix{<:Real},
    t::Vector{<:Real}
)   
    vec_field_temp(x,p,t) = vec(vec_field(x,t)) # form for ODE package # note vec() to turn from row vector to column vector (required form)
    
    # function vec_field_temp!(dx,x,p,t)
    #     temp = vec_field(x,p,t)
    #     dx[1] = temp[1]
    #     dx[2] = temp[2]
    # end
    
    # return Trajectory (cut-off values outside indicator function, if provided)
    n_t = length(t)
    n_s = size(X_0,1)
    d = size(X_0,2)

    p = Progress(n_s,desc="Calculating trajectories...",barglyphs=BarGlyphs("[=> ]"),barlen=50,color=:yellow)

    X = Array{Float64}(undef,n_s,n_t,d) # change this # not sure if I like the representation order; (n_s,d,n_t) more consistent?
    for i = 1:n_s
        next!(p)
        local u0 = vec(X_0[i,:])
        local prob = ODEProblem(vec_field_temp,u0,(t[1],t[end]))
        local sol = solve(prob,AutoTsit5(Rosenbrock23()),reltol=1e-6,saveat=t)
        for j = 1:n_t
            X[i,j,:] = sol.u[j]
        end
    end

    return Trajectory(X,t)
end

function solve_trajectory(
    vec_field,
    X::Matrix{Float64},
    t::AbstractRange
)   
    return solve_trajectory(vec_field,X,collect(t))
end

function clip(traj::Trajectory,indicator::Function)
    inds = Vector{Int32}(undef,0) # store trajectory indices that leave the indicated domain
    for i = 1:size(traj.X,1)
        for j = 1:length(traj.t)
            if (!indicator(traj.X[i,j,:])) # "if trajectory leaves the indicated domain"
                push!(inds,i) # add index to inds
                break
            end
        end
    end
    return Trajectory(traj.X[Not(inds),:,:],traj.t) # ignore the trajectory indices that leave the indicated domain
end

function clip(points::Array{<:Real},indicator::Function)
    keep = Vector{Int32}(undef,0)
    for i = 1:size(points,1)
        if (indicator(points[i,:]))
            push!(keep,i)
        end
    end
    return points[keep,:] # can we change this to a mutation? i.e., clip!()
end