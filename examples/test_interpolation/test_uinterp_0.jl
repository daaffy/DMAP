# using Revise
using DMAP
using Plots
using LinearAlgebra

# want to approximate the following function
function g(x,y)
    return [sin(x*y);2*sin(x*y)]
end

# let's sample the function g on a subdomain:
# xy_0 = [[x y] for x = 0:.1:1, y = 0:.1:1] # even grid
xy_0 = [[x^1.5 y^1.5] for x = 0:.1:1, y = 0:.1:1] # uneven grid example
X_0 = vcat(xy_0...)
n_0 = size(X_0,1)
z_0 = real([g(X_0[i,1],X_0[i,2]) for i = 1:n_0])

# display(scatter(X_0[:,1], X_0[:,2]))

eps = calculate_eps(X_0) # rbf "width", play around with this

# coord_vec = X_0
# zval = z_0

# n = size(coord_vec,1)
# ker_matrix = [rbf(coord_vec[i,:],coord_vec[j,:],eps) for i = 1:n, j = 1:n]
# # coeff_vec = ker_matrix\zval
# coeff_vec = inv(ker_matrix)*zval
# U = [0;0]
# temp = sum(coeff_vec.*[rbf(U,coord_vec[i,:],eps) for i = 1:n],dims=1)
# length(temp) == 1 ? sum(temp) : temp

# f = uinterp_2D(X[:,1],X[:,2],z,eps) # interpolated function defined on unstructured points
f = uinterp_0(X_0,z_0,eps)

# let's evaluate the interpolant f on the domain:
xy = [[x y] for x = 0:.01:1, y = 0:.01:1]
X = vcat(xy...)
# X = X_0
n = size(X,1)
z = [g(X[i,1],X[i,2]) for i = 1:n]
# z = [real(f(X[i,:]))[1]-g(X[i,1],X[i,2]) for i = 1:n]
# z = [real(f(X[i,:]))-g(X[i,1],X[i,2]) for i = 1:n] # calculate error between interpolated function and actual function

# calculate error
# (!) needs to weighted for accurate integral
# err = sum(z.^2) 

# test_uinterp_0(f,X_0,z_0)

# plot
# display(scatter(X[:,1], X[:,2], zcolor=z, c=:jet, aspectratio=1))
