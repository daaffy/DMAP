using Revise
using DMAP
using Plots
using LinearAlgebra

# we want to approximate the following function (vector valued function, in general)
function g(x,y)
    return [sin(x*y);2*sin(x*y)]
end

# let's sample the function g on a subdomain:
q = 2 # = 1 for uniformly spaced grid
xy_0 = [[x^q y^q] for x = 0:.1:1, y = 0:.1:1] # uneven grid example
X_0 = vcat(xy_0...)
n_0 = size(X_0,1)
z_0 = [g(X_0[i,1],X_0[i,2]) for i = 1:n_0]

# create interpolated function
f = uinterp_0(X_0,z_0,"uniform")

# evaluate g on a fine grid for comparison with f
xy = [[x y] for x = 0:.01:1, y = 0:.01:1]
X = vcat(xy...)
n = size(X,1)
z = [g(X[i,1],X[i,2]) for i = 1:n]

err1 = test_uinterp_0(f,X_0,z_0) # get max error on the original sample grid
err2 = test_uinterp_0(f,X,z) # get max error on the fine grid

println("max error on sampled grid: ",err1)
println("max error on fine grid: ", err2)
