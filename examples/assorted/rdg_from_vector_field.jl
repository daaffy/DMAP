# using Revise
using CoherentStructures
using LinearAlgebra
using Plots
# using Dierckx # unstructured interpolation

# function interp_vel_field!(du,u,p,t,interpf)
#     function vel_field!(du,u,p,t)
#         du = interpf(u[1],u[2])
#     end
# end

# -------------------------------------------------------------------------------------------------------
# define velocity field 
    function s(t)
        if (t < 0)
            return 0
        elseif(t > 1)
            return 1
        else
            return (t.^2).*(3-2*t);
        end
    end

    function flow_field_double_gyre(x,y,t)
        # flow field for the double gyre; calculated as the curl of the double gyre vector stream function
        temp = s(t);
        return [pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)]*[1-temp; temp];
    end

# -------------------------------------------------------------------------------------------------------
# (this part is invisible in practice)

    flow_field = flow_field_double_gyre
    xy = [[x y] for x = 0:.2:1, y = 0:.2:1]
    X = vcat(xy...)
    n = size(X,1)
    # z = real([flow_field(X[i,1],X[i,2],0) for i = 1:n])
    global z = Matrix{Float64}(undef,0,2);
    for i = 1:n
        global z = cat(dims=1,z,flow_field(X[i,1],X[i,2],0)')
    end

# -------------------------------------------------------------------------------------------------------

# interpolate on sampled velocity field; e.g., ansys output
interp = uinterp_2D(X[:,1],X[:,2],z,0.2)

# evaluate interp on a grid of points
xy = [[x y] for x = 0:.05:1, y = 0:.05:1]
X = vcat(xy...)
n = size(X,1)
global z = Matrix{Float64}(undef,0,2);
for i = 1:n
    global z = cat(dims=1,z,interp(X[i,1],X[i,2]))
    # global z = cat(dims=1,z,flow_field(X[i,1],X[i,2],0)')
    # global z = cat(dims=1,z,interp(X[i,1],X[i,2])-flow_field(X[i,1],X[i,2],0)')
end

plot_vector_field(X[:,1],X[:,2],z/50)
# display(scatter(X[:,1],X[:,2],zcolor=z[:,1].^2 + z[:,2].^2))

# temp = interp_vel_field!(0,0,0,0,interp)

# -------------------------------------------------------------------------------------------------------



