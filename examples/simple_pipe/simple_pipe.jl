using Revise
using CoherentStructures
using LinearAlgebra
using Plots
using CSV
using DataFrames
using DifferentialEquations

# -------------------------------------------------------------------------------------------------------
# extract discrete velocity field to X, dX multi-dimensional arrays.
# ASSUMPTION: velocity field is evaluated on a constant grid

# input_dir = "/Users/jackh/Documents/Honours/DMAp/DMAP_0/examples/simple_pipe/input/"
# input_files = readdir(input_dir)
# discrete_steps = length(input_files)

# # this can surely be cleaned up...
# i = 1
# df = CSV.read(input_dir * input_files[i], DataFrame)
# df = Array(df)
# n_points = size(df,1)
# #

# X = Array{Float64}(undef,n_points,3)
# for j = 1:n_points
#     global X[j,:] = df[j,2:4]'
# end

# dX = Array{Float64}(undef,n_points,3,discrete_steps)
# for i = 1:discrete_steps
#     global df = CSV.read(input_dir * input_files[i], DataFrame)
#     global df = Matrix(df)

#     for j = 1:n_points
#         # check that grid doesn't change through time (within margin; use initial grid anyway)
#         # if (X[j,:] != df[j,2:4]')
#         #     print(-df[j,2:4]')
#         #     return
#         # end
#         global dX[j,:,i] = df[j,5:7]'
#     end
# end

# for this example we need to manually construct the time vector since the data is not supplied
t = LinRange(0,1,discrete_steps)

# -------------------------------------------------------------------------------------------------------
# build interpolated velocity field

eps=calculate_eps(X)
# eps = 0.0002798035466992911
# f = uinterp_2D(X,dX[:,:,1],0.0002798035466992911*1E-2)
# f = uinterp(X,t,dX,eps)

# coords = X
# times = t
# zvals = dX

# n = size(coords,1)
# ker_matrix = [iso_kernel(coords[i,:],coords[j,:],eps) for i = 1:n, j = 1:n]
# coeff_vec = similar(zvals)
# inv_ker = inv(ker_matrix)
# for i = 1:length(times)
#     coeff_vec[:,:,i] = inv_ker*zvals[:,:,i]
# end
# function eval(X,t)
#     # interpolate through (simple linear interpolation)
#     if (t <= times[1])
#         coeff_vec = coeff_vec[:,:,1]
#     elseif(t >= times[length(times)])
#         coeff_vec = coeff_vec[:,:,length(times)]
#     else
#         for ti = 1:length(times)
#             if (times[ti] > t)
#                 del1 = (times[ti]-t)/((times[ti]-times[ti-1]))
#                 del2 = (t-times[ti-1])/((times[ti]-times[ti-1]))
#                 coeff_vec = del1*coeff_vec[:,:,ti-1]+del2*coeff_vec[:,:,ti]
#                 break
#             end
#         end
#     end

#     return sum(coeff_vec.*[iso_kernel(X',coord_vec[i,:],eps) for i = 1:n],dims=1)[1]
# end


# -------------------------------------------------------------------------------------------------------
# integrate trajectories


# -------------------------------------------------------------------------------------------------------
# test plots
