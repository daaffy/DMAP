using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis
using SparseArrays
using Arpack
using StaticArrays
using ProfileView
using BenchmarkTools
using NearestNeighbors
# using JLD
# using Profile

# -------------------------------------------------------------------------------------------------------
# define velocity field 
# function s(t)
#     if (t < 0)
#         return 0
#     elseif(t > 1)
#         return 1
#     else
#         return (t.^2).*(3-2*t);
#     end
# end

# function flow_field_double_gyre(x,y,t)
#     # flow field for the double gyre; calculated as the curl of the double gyre vector stream function
#     temp = (t.^2).*(3-2*t)
#     return [pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)]*[1-temp; temp];
# end

"""
    # Rotating double gyre flow field.
For 0 <= t <= 1.
"""
# function vec_field(u,t)
#     temp = (t.^2).*(3-2*t)
#     # don't think static matrix/vector has helped much here...
#     return vec(SMatrix{2,2}([pi*sin(2*pi*u[1])*cos(pi*u[2]) 2*pi*sin(pi*u[1])*cos(2*pi*u[2]); -2*pi*cos(2*pi*u[1])*sin(pi*u[2]) -pi*cos(pi*u[1])*sin(2*pi*u[2])])*SVector{2}([1-temp; temp]));
# end

# =tt->(tt.^2).*(3-2*tt)


# s = t -> sin(0.5*pi*t)

# forward = 1 # -1 for backwards
s = t -> (t^2)*(3-2*t) # forward
# s = t -> 0 # backward
function vec_field!(du,u,p,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

# test on basic dynamics
# function vec_field!(du,u,p,t)
#     x,y = u
#     du[1] = sqrt(x)
#     du[2] = y
# end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

res = 200
X_0 = create_grid([0 1],[0 1],[res res])
t = (0:0.01:1)

@time traj = solve_ensemble(vec_field!,X_0,t,s)

# -------------------------------------------------------------------------------------------------------
# let's watch what happens to δ-ball of points over a short period of time...







# # # --- VISUALISE ADVECTION IN HERE SOMEWHERE
# point_ind = 6320
# point_ind = 7321
# t1_ind = 12
# t2_ind = 15
# balltree = BallTree(traj.X[:,t1_ind,:]') 
# # idxs,_ = get_δball(balltree,point,δ=0.03)
# print("point coordinates: "*string(traj.X[point_ind,t1_ind,:]))
# # method
# centre_0 = traj.X[point_ind,t1_ind,:]
# centre_1 = traj.X[point_ind,t2_ind,:]

# # idxs,_ = get_δball(balltree,centre_0,δ=0.010)
# # δ_ball = traj.X[idxs,t1_ind,:]
# # adv_δ_ball = traj.X[idxs,t2_ind,:]
# # J,perts,advected_perts = calculate_jacobian(δ_ball,adv_δ_ball,centre_0,centre_1)
# # test = (J*perts')'

# idxs,_ = get_δball(balltree,centre_0,δ=0.025)
# # centre_0 = traj.X[point_ind,t1_ind,:]
# # centre_1 = traj.X[point_ind,t2_ind,:]
# δ_ball = traj.X[idxs,t1_ind,:]
# adv_δ_ball = traj.X[idxs,t2_ind,:]
# _,perts,advected_perts = calculate_jacobian(δ_ball,adv_δ_ball,centre_0,centre_1)

# scatter()
# pyplot(markershape = :circle)
# scatter!(centre_0[1].+perts[:,1],centre_0[2].+perts[:,2],aspectratio=1,label="t_0",color="black",)# perterbations
# # pyplot(markershape = :diamond)
# scatter!(centre_1[1].+advected_perts[:,1],centre_1[2].+advected_perts[:,2],aspectratio=1,label="t_0+T",color="black",markershape = :diamond)# Jacobian*perterbations
# # display(scatter!(advected_perts[:,1],advected_perts[:,2],aspectratio=1)) # perterbations post advection



# centre = centre_0
# advected_centre = centre_1



# --- build func. here

# inputs:
# t_inds = [1;50]




# t1_ind = 1
# t2_ind = 50
# points = traj.X[:,t1_ind,:]
# advected_points = traj.X[:,t2_ind,:]
# balltree = BallTree(points') 

# λ_vec = vec(zeros(size(points,1),1))
# avg_neighbours = 0
# for i = axes(points,1)
#     centre = points[i,:]
#     advected_centre = advected_points[i,:]
#     idxs,_ = get_δball(balltree,centre,δ=0.0075)
#     global avg_neighbours += length(idxs)
#     δ_ball = points[idxs,:]
#     advected_δ_ball = advected_points[idxs,:]
#     J,_,_ = calculate_jacobian(δ_ball,advected_δ_ball,centre,advected_centre)
#     λs,_ = eigen(J'*J)
#     λ_vec[i] = max(λs...)
# end
# avg_neighbours /= size(points,1)
# ftle_field = 1/(traj.t[t2_ind]-traj.t[t1_ind])*log.(sqrt.(λ_vec))

ftle_field,_ = ftle(traj,[1;50],0.015)

# display(scatter(X_0[:,1],X_0[:,2],zcolor=ftle_field, c=:jet, aspectratio=1))
# print("avg_neighbours: "*string(avg_neighbours))
display(heatmap(LinRange(0,1,res),LinRange(0,1,res),reshape(ftle_field,res,res)',interpolate=true,aspectratio=1,c =:magma))











# colour candidates
# :ice, :flag is fun, 

# J = [1 0; 0 1]
# ball_1 = J*δ_ball'

# # the initial δ-ball looks like:
# display(scatter(δ_ball[:,1],δ_ball[:,2],aspectratio=1))

# # what does the advected δ-ball look like?
# display(scatter!(traj.X[idxs,t2_ind,1],traj.X[idxs,t2_ind,2],aspectratio=1))

# ftle calculation


# -------------------------------------------------------------------------------------------------------
# # plot
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))
# display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))