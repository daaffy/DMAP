
using Revise
using DMAP
using LinearAlgebra
using Plots
using DifferentialEquations
using SparseArrays
using Arpack
using Distances
# using JLD
# using Profile

function iso_kernel(
    x_1,
    x_2,
    epsilon,
    metric = Euclidean()
)   
    # thresh = 0
    # temp = exp(-metric(x_1,x_2)^2/epsilon)
    # return temp >= thresh ? temp : 0;
    return exp(-metric(x_1,x_2)^2/epsilon)
end


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

function vec_field(u,t)
    return flow_field_double_gyre(u[1],u[2],t)
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

X_0 = create_grid([0 1],[0 1],[40 40])
t = (0:0.05:1)
traj = solve_trajectory(vec_field,X_0,t)

metric = Euclidean()
epsilon = nndist(traj.X[:,1,:])/sqrt(2) # according to heuristic; see summary
threshold = 0.9

n_s = size(traj.X, 1)
n_t = size(traj.X, 2)
alpha = 1 # see Coifman/Lafon
f(x) = (iso_kernel(x,0,epsilon) >= threshold ? iso_kernel(x,0,epsilon) : 0) # can we do it without writing iso_kernel twice?

k = 10

temp = traj.t[2:end]-traj.t[1:end-1]
coeff = 0.5*[temp[1]; temp[1:end-1]+temp[2:end]; temp[end]]
global P_final = spzeros(n_s,n_s)
for ti = 1:k
    dist_matrix = pairwise(metric,traj.X[:,ti,:],traj.X[:,ti,:],dims=1)
    ker_matrix = broadcast(f,dist_matrix)
    ker_matrix = sparse(ker_matrix)
    q_bar = sum(ker_matrix, dims=2)/n_s
    d_bar = [ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha) for i = 1:n_s, j = 1:n_s]
    d_barsum = sum(d_bar, dims=1)/n_s;
    p_bar = [d_bar[i,j]/d_barsum[j] for i = 1:n_s, j = 1:n_s]
    global P_final += coeff[ti]*p_bar/n_s
end

P_final /= traj.t[end]-traj.t[1]

# sparse(P_final)
A = sparse(P_final)
print(length(A.nzval)/(n_s*n_s))
spy(A,ms=3)