using Revise
using DMAP
using Plots
using LinearAlgebra
# using ScatteredInterpolation

# we want to approximate the following function (vector valued function, in general)
function g(x,y)
    return [sin(x*y);2*sin(x*y)]
end

# let's sample the function g on a subdomain:
# q = 1 # = 1 for uniformly spaced grid
# xy_0 = [[x^q y^q] for x = 0:.1:1, y = 0:.1:1] # uneven grid example
# X_0 = vcat(xy_0...)
# n_0 = size(X_0,1)
# z_0 = vcat([g(X_0[i,1],X_0[i,2]) for i = 1:n_0]...) # Vector{Vector{Float64}}

X_0 = [0.0 0; 0 1; 1 0]
z_0 = Array{Float64}(undef,size(X_0,1),size(X_0,2),3)
z_0[:,:,1] = [0.0 0; 0 0; 0 0]
z_0[:,:,2] = [1 1; 0.5 0.5; 0.4 0.4]
z_0[:,:,3] = [10 10; 10 10; 10 10]
t = [0;1;2]

# # test uinterp_0 of a time-slice
# itp = uinterp_0(X_0,z_0[:,:,1],"ubk")
# eval_points = X_0[1:2,:]
# println(evaluate_itp(itp,eval_points))
# println(z_0[1:2,:,1])

dvf = DiscreteVelocityField(X_0,z_0,t)

# build function here ------
# itps = [uinterp_0(dvf.X,dvf.dX[:,:,t_i],"ubk") for t_i = 1:length(dvf.t)]
# t = 1.2 # input to eval function
# points = [0 0]
# if (t <= dvf.t[1])
#     return evaluate_itp(itps[1],points)
# elseif((t >= dvf.t[end]))
#     return evaluate_itp(itps[end],points)
# else
#     for ti = 1:length(dvf.t)
#         if (dvf.t[ti] > t)
#             del1 = (dvf.t[ti]-t)/((dvf.t[ti]-dvf.t[ti-1]))
#             del2 = (t-dvf.t[ti-1])/((dvf.t[ti]-dvf.t[ti-1]))
#             return del1*evaluate_itp(itps[ti-1],points)+del2*evaluate_itp(itps[ti],points)
#             break
#         end
#     end
# end
# determine which 
# --------------------------

f = uinterp(dvf,"nn")

t = 1.5 # input to eval function
points = [0 0; 0 1; 1 0; 1 0.1]
println(f(points,t))

# looks like it works well!








# create interpolated function
# f = uinterp_0(X_0,z_0,"ubk")

# evaluate g on a fine grid for comparison with f
# xy = [[x y] for x = 0:.01:1, y = 0:.01:1]
# X = vcat(xy...)
# n = size(X,1)
# z = [g(X[i,1],X[i,2]) for i = 1:n]



# err1 = test_uinterp_0(f,X_0,z_0) # get max error on the original sample grid
# err2 = test_uinterp_0(f,X,z) # get max error on the fine grid

# println("max error on sampled grid: ",err1)
# println("max error on fine grid: ", err2)

# using ScatteredInterpolation.jl:
# itp = interpolate(Gaussian(1),X_0',z_0)

# 


# https://eljungsk.github.io/ScatteredInterpolation.jl/stable/

