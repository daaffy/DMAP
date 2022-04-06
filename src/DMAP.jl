module DMAP

using Distances:Euclidean

# plotting
using Plots

# export:
export dynamic_laplacian, calculate_trajectories, SEBA, animate_scatter, plot_vector_field, Trajectory, velocity_field, iso_kernel

# assisting with each stage of the workflow
include("loader.jl")
include("seba.jl")
include("nndist.jl")
include("interpolate.jl")

# DATA STRUCTURES

struct DiscreteVelocityField
    X::Array{Float64} # node coordinates
    dX::Array{Float64} # node velocities
    t::Vector{Float64} # time steps

    # perform checks
    # http://web.mit.edu/julia_v0.6.2/julia/share/doc/julia/html/en/manual/constructors.html
    function DiscreteVelocityField(X,dX,t)
        # e.g.
        if (size(dX,3)!=length(t))
            error("incorrect dimensions")
        end
        new(X,dX,t)
    end
end

# create trajectory structure
# 
# struct Trajectory
#     x::Array
#     y::Array
#     z::Array
#     t::Vector
# end

# to-do;; vectorise /  try to increase the computational speed
# ;; make sparse (zero threshold)
function dynamic_laplacian(
    # to-do: trajectory struct (can extract spatial dimensions from this)
    X_traj::Array{Float64},
    t,
    metric = Euclidean(),
    epsilon = nndist(X_traj[:,1,:])/sqrt(2) # according to heuristic; see summary
    # epsilon = 0.25/sqrt(2)
)   
    n_s = size(X_traj, 1)
    n_t = size(X_traj, 2)
    alpha = 1

    global P_final = zeros(n_s,n_s)
    for ti = 1:n_t
        # ker_matrix = [iso_kernel(X_traj[i,ti,:],X_traj[j,ti,:],epsilon,metric) for i = 1:n_s, j = 1:n_s]
        ker_matrix = Array{Float64}(undef, n_s, n_s);
        for i = 1:n_s
            for j = 1:n_s
                # global ker_matrix[i,j] = iso_kernel([x_traj[i,ti];y_traj[i,ti]],[x_traj[j,ti];y_traj[j,ti]],epsilon)
                global ker_matrix[i,j] = iso_kernel(X_traj[i,ti,:],X_traj[j,ti,:],epsilon,metric)
            end
        end
        # print(sum(sum(ker_matrix>0 ? 1 : 0))/(n_s*n_s))
        

        q_bar = sum(ker_matrix, dims=2)/n_s;

        d_bar = zeros(n_s,n_s);
        for i = 1:n_s
            for j = 1:n_s
                global d_bar[i,j] = ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha); 
            end
        end
        d_barsum = sum(d_bar, dims=1)/n_s;

        global p_bar = zeros(n_s,n_s);
        for i = 1:n_s
            for j = 1:n_s
                global p_bar[i,j] = d_bar[i,j]/d_barsum[j];
            end
        end

        global P_final += p_bar/n_s
    end


    P_final /= n_t # THIS NEEDS TO BE GENERALISED

    return P_final
end

function iso_kernel(
    x_1,
    x_2,
    epsilon,
    metric = Euclidean()
)   
    thresh = 0
    temp = exp(-metric(x_1,x_2)^2/epsilon)
    return temp >= thresh ? temp : 0;
end

function calculate_trajectories(
    x,
    y,
    t,
    flow_field
)
    n_p = length(x);
    max_del_t = (t[2]-t[1])*10^(-2); # UPDATE THIS TO MORE GENERAL
    if (length(t) == 2)
        time_steps = LinRange(t[1],t[2],ceil(Int,(t[2] - t[1])/max_del_t)+1);
        n_t = length(time_steps);
        x_int = Array{Float64}(undef, n_p, n_t);
        x_int[:,1] = x;
        y_int = Array{Float64}(undef, n_p, n_t);
        y_int[:,1] = y;
        for i_t = 1:n_t-1
            for i_p = 1:n_p
                
                next_xy = [x_int[i_p,i_t];y_int[i_p,i_t]] + flow_field(x_int[i_p,i_t],y_int[i_p,i_t],time_steps[i_t])*(time_steps[i_t+1] - time_steps[i_t])
                x_int[i_p,i_t+1] = mod(next_xy[1],20); # note mod due specificially for BJ periodic boundary conditions; need to make this more flexible!
                y_int[i_p,i_t+1] = next_xy[2];
            end
        end
        x_out = [x_int[:,1] x_int[:,end]];
        y_out = [y_int[:,1] y_int[:,end]];

        return x_out, y_out
    elseif(length(t) > 2)
        n_t = length(t); # NOTE: different definition than above
        x_out = Array{Float64}(undef, n_p, n_t);
        y_out = Array{Float64}(undef, n_p, n_t);
        x_out[:,1] = x;
        y_out[:,1] = y;
        for i_t = 1:n_t-1
            traj_out = calculate_trajectories(x_out[:,i_t], y_out[:,i_t], [t[i_t];t[i_t+1]],flow_field);
            x_out[:,i_t+1] = traj_out[1][:,2];
            y_out[:,i_t+1] = traj_out[2][:,2];
        end 
        return x_out, y_out
    end
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