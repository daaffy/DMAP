module DMAP

using Distances:Euclidean

# plotting
using Plots

# export:
export DiscreteVelocityField, Trajectory # need to export?

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
        if (size(dX,3)!=length(t) || size(dX,1)!=size(X,1) || size(dX,2)!=size(X,2))
            error("Dimensions mismatch.")
        end
        new(X,dX,t)
    end
end

# assisting with each stage of the workflow
include("./loader.jl")
include("./trajectories.jl")
include("./seba.jl")
include("./nndist.jl")
include("./methods.jl")
include("./interpolate.jl")
include("./exporter.jl")
include("./findepsilon.jl")

# function iso_kernel(
#     x,
#     epsilon
# )   
#     # thresh = 0
#     # temp = exp(-metric(x_1,x_2)^2/epsilon)
#     # return temp >= thresh ? temp : 0;
#     return exp(-x^2/epsilon)
# end

# function animate_scatter(x_traj,y_traj,vec)
#     n_t = size(x_traj,2)
#     anim = @animate for i ∈ 1:n_t
#             display(scatter(x_traj[:,i], y_traj[:,i], zcolor=vec, c=:jet, aspectratio=1))
#     end
#     gif(anim, "anim.gif", fps = 15)
# end

# function plot_vector_field(x,y,vector_field)
#     display(quiver(x,y, quiver=(vector_field[:,1],vector_field[:,2])))
# end

end