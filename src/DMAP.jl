module DMAP

using Distances:Euclidean

# plotting
using Plots

# export:
export dynamic_laplacian, calculate_trajectories, SEBA, animate_scatter, plot_vector_field, velocity_field, iso_kernel
export DiscreteVelocityField

struct Trajectory
    X::Array{Float64,3}
    t::Vector{Float64}
end

struct DiscreteVelocityField
    X::Array{Float64,2} # node coordinates
    dX::Array{Float64} # node velocities
    t::Vector{Float64} # time steps

    # perform checks
    # http://web.mit.edu/julia_v0.6.2/julia/share/doc/julia/html/en/manual/constructors.html
    function DiscreteVelocityField(X,dX,t) 
        # e.g.
        if (size(dX,3)!=length(t)) # ...
            error("incorrect dimensions")
        end
        new(X,dX,t)
    end
end

# assisting with each stage of the workflow
include("./loader.jl")
include("./trajectories.jl")
include("./seba.jl")
include("./nndist.jl")
include("./interpolate.jl")
include("./exporter.jl")

function dynamic_laplacian(
    traj::Trajectory,
    metric = Euclidean(),
    epsilon = nndist(traj.X[:,1,:])/sqrt(2) # according to heuristic; see summary
)   
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)
    alpha = 1 # see Coifman/Lafon
    f(x) = iso_kernel(x,0,epsilon)

    temp = traj.t[2:end]-traj.t[1:end-1]
    coeff = 0.5*[temp[1]; temp[1:end-1]+temp[2:end]; temp[end]]
    global P_final = zeros(n_s,n_s)
    for ti = 1:n_t
        dist_matrix = pairwise(metric,traj.X[:,ti,:],traj.X[:,ti,:],dims=1)
        ker_matrix = broadcast(f,dist_matrix)
        q_bar = sum(ker_matrix, dims=2)/n_s
        d_bar = [ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha) for i = 1:n_s, j = 1:n_s]
        d_barsum = sum(d_bar, dims=1)/n_s;
        p_bar = [d_bar[i,j]/d_barsum[j] for i = 1:n_s, j = 1:n_s]
        global P_final += coeff[ti]*p_bar/n_s
    end

    P_final /= traj.t[end]-traj.t[1]

    return P_final
end

function iso_kernel(
    x_1,
    x_2,
    epsilon,
    metric = Euclidean()
)   
    # thresh = 0
    # temp = exp(-metric(x_1,x_2)^2/epsilon)
    # return temp >= thresh ? temp : 0;
    return exp(-metric(x_1,x_2)^2/epsilon)
end

# function iso_kernel(
#     x,
#     epsilon
# )   
#     # thresh = 0
#     # temp = exp(-metric(x_1,x_2)^2/epsilon)
#     # return temp >= thresh ? temp : 0;
#     return exp(-x^2/epsilon)
# end

function animate_scatter(x_traj,y_traj,vec)
    n_t = size(x_traj,2)
    anim = @animate for i ∈ 1:n_t
            display(scatter(x_traj[:,i], y_traj[:,i], zcolor=vec, c=:jet, aspectratio=1))
    end
    gif(anim, "anim.gif", fps = 15)
end

function plot_vector_field(x,y,vector_field)
    display(quiver(x,y, quiver=(vector_field[:,1],vector_field[:,2])))
end

end