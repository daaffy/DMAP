using Revise
using DMAP
using LinearAlgebra
using Plots
using CSV
using DataFrames
using DifferentialEquations

# -------------------------------------------------------------------------------------------------------]
# extract discrete velocity field to X, dX multi-dimensional arrays.

input_dir = "/Users/jackh/Documents/Honours/DMAp/DMAP/examples/simple_pipe/input/"
input_files = readdir(input_dir)

dvf = load_dvf(input_files,input_dir)

# -------------------------------------------------------------------------------------------------------
# build interpolated velocity field
# eps = calculate_eps(dvf.X)

# f = uinterp_0(dvf.X,dvf.dX[:,1,1],eps)

# eps = 0.0002798035466992911
# f = uinterp_2D(X,dX[:,:,1],0.0002798035466992911*1E-2)
# f = uinterp(X,t,dX,eps)

coord_vec = dvf.X
eps = 0.01
# times = t
zval = dvf.dX[:,:,1]

n = size(coord_vec,1)
ker_matrix = [rbf(coord_vec[i,:],coord_vec[j,:],eps) for i = 1:n, j = 1:n]
# coeff_vec = ker_matrix\zval
coeff_vec = inv(ker_matrix)*zval
U = dvf.X[1,:]
temp = sum(coeff_vec.*[rbf(U,coord_vec[i,:],eps) for i = 1:n],dims=1)
length(temp) == 1 ? sum(temp) : temp


# -------------------------------------------------------------------------------------------------------
# integrate trajectories


# -------------------------------------------------------------------------------------------------------
# test plots
