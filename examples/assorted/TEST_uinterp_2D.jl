# using Revise
using CoherentStructures
using Plots



# want to approximate the following function
function g(x,y)
    return sin(x*y)
end

# let's sample the function g on a subdomain:
xy = [[x y] for x = -2:.3:2, y = -2:.3:2]
X = vcat(xy...)
n = size(X,1)
z = real([g(X[i,1],X[i,2]) for i = 1:n])



# function calculate_eps(X)
#     n = size(X,1)
#     dist_matrix = [norm(X[i,:]-X[j,:]) for i = 1:n, j = 1:n]
#     # d_vec = [min(filter(x -> x == 0,dist_matrix[i,1])) for i = 1:n]
#     # return 0.815*sum(d_vec)/n # Hardy's, obtained from S. Rippa paper.
#     return 0.815*0.3
# end

eps = calculate_eps(X) # rbf "width", play around with this

# f = uinterp_2D(X[:,1],X[:,2],z,eps) # interpolated function defined on unstructured points
f = uinterp_2D(X,z,eps)

# let's evaluate the interpolant f on the domain:
xy = [[x y] for x = -4:.1:4, y = -4:.1:4]
X = vcat(xy...)
n = size(X,1)
# z = [real(f(X[i,:])) for i = 1:n] # calculate error between interpolated function and actual function
z = [real(f(X[i,:])) for i = 1:n]

# plot
display(scatter(X[:,1], X[:,2], zcolor=z, c=:jet, aspectratio=1))
