module DMAP

using Distances:Euclidean

# plotting
using Plots

# export:
export dynamic_laplacian, dynamic_laplacian0, calculate_trajectories, SEBA, animate_scatter, plot_vector_field, velocity_field, iso_kernel
export DiscreteVelocityField

struct DiscreteVelocityField
    X::Array{Float64,2} # node coordinates
    dX::Array{Float64} # node velocities
    t::Vector{Float64} # time steps

    # perform checks
    # http://web.mit.edu/julia_v0.6.2/julia/share/doc/julia/html/en/manual/constructors.html
    function DiscreteVelocityField(X,dX,t) # ...
        # e.g.
        if (size(dX,3)!=length(t))
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

# to-do;; vectorise /  try to increase the computational speed
# ;; make sparse (zero threshold)
# function dynamic_laplacian(
#     # to-do: trajectory struct (can extract spatial dimensions from this)
#     traj::Trajectory,
#     metric = Euclidean(),
#     epsilon = nndist(traj.X[:,1,:])/sqrt(2) # according to heuristic; see summary
#     # epsilon = 0.25/sqrt(2)
# )   
#     X_traj = traj.X
#     t = traj.t
#     n_s = size(X_traj, 1)
#     n_t = size(X_traj, 2)
#     alpha = 1

#     global P_final = zeros(n_s,n_s)
#     for ti = 1:n_t
#         # ker_matrix = [iso_kernel(X_traj[i,ti,:],X_traj[j,ti,:],epsilon,metric) for i = 1:n_s, j = 1:n_s]
#         ker_matrix = Array{Float64}(undef, n_s, n_s);
#         for i = 1:n_s
#             for j = 1:n_s
#                 # global ker_matrix[i,j] = iso_kernel([x_traj[i,ti];y_traj[i,ti]],[x_traj[j,ti];y_traj[j,ti]],epsilon)
#                 global ker_matrix[i,j] = iso_kernel(X_traj[i,ti,:],X_traj[j,ti,:],epsilon,metric)
#             end
#         end
        

#         q_bar = sum(ker_matrix, dims=2)/n_s;

#         d_bar = zeros(n_s,n_s);
#         for i = 1:n_s
#             for j = 1:n_s
#                 global d_bar[i,j] = ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha); 
#             end
#         end
#         d_barsum = sum(d_bar, dims=1)/n_s;

#         global p_bar = zeros(n_s,n_s);
#         for i = 1:n_s
#             for j = 1:n_s
#                 global p_bar[i,j] = d_bar[i,j]/d_barsum[j];
#             end
#         end

#         global P_final += p_bar/n_s
#     end


#     P_final /= n_t # THIS NEEDS TO BE GENERALISED

#     return P_final
# end

function dynamic_laplacian(
    traj::Trajectory,
    metric = Euclidean(),
    epsilon = nndist(traj.X[:,1,:])/sqrt(2) # according to heuristic; see summary
)   
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)
    alpha = 1 # see Coifman/Lafon
    f(x) = iso_kernel(x,epsilon)

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

function iso_kernel(
    x,
    epsilon
)   
    # thresh = 0
    # temp = exp(-metric(x_1,x_2)^2/epsilon)
    # return temp >= thresh ? temp : 0;
    return exp(-x^2/epsilon)
end

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

# function velocity_field(
#     a::Int
# )
#     print("one")
# end

# function velocity_field(
#     analytic_field::Function,
#     b::Int
# )
#     print("two")
# end


end