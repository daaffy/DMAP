using Revise
using CoherentStructures
using LinearAlgebra
using Plots
using DifferentialEquations
using Distances


# DEFINE FLOW FIELD
function vec_field(u,p,t)
    A = sqrt(3)
    B = sqrt(2)
    C = 1
    return [A*sin(u[3])+C*cos(u[2]);
            B*sin(u[1])+A*cos(u[3]);
            C*sin(u[2])+B*cos(u[1])]
end

del = 0.1
x_ticks = 0:del:1; nx = length(x_ticks)
y_ticks = 0:del:1; ny = length(y_ticks)
z_ticks = 0:del:1; nz = length(z_ticks)
# xyz = [[x y z] for x = 0:del:1, y = 0:del:1, z = 0:del:1]
c = 1
X = Array{Float64}(undef,0,3)
map2 = Array{Float64}(undef,nx,ny,nz)
for i = 1:nx, j = 1:ny, k = 1:nz
    global X = cat(X,[x_ticks[i] y_ticks[j] z_ticks[k]],dims=1)
    map2[i,j,k] = c
    global c += 1
end
# X = vcat(xyz...)
n_s = size(X,1)

t = (0:0.01:1)
tspan = (t[1],t[end])
T = t[end]-t[1]
n_t = length(t)

X_traj = Array{Float64}(undef,n_s,n_t,3)
for i = 1:n_s
    local u0 = [X[i,1];X[i,2];X[i,3]]
    local prob = ODEProblem(vec_field,u0,tspan)
    local sol = solve(prob,saveat=t)
    for j = 1:n_t
        X_traj[i,j,:] = sol.u[j]
    end
end

function peuc(x,y)
    return peuclidean(x,y,[1,1,1])
end


@time P = dynamic_laplacian(X_traj,t,peuc)

# calculate eigens 
λ, v = eigen(P')

# -------------------------------------------------------------------------------------------------------
# SEBA post-processing
r = 6
V = real(v[:,end:-1:end-r+1]) 
S, R = SEBA(V)

# mask to obtain a z-plane slice:
layer = 1
mask = zeros(n_s,nx*ny)
inds = vec(map2[:,:,layer])
for j = 1:length(inds)
    mask[Int(inds[j]),j] = 1
end
mask = mask'

display(scatter(mask*X_traj[:,1,1], mask*X_traj[:,1,2], zcolor=mask*real(v[:,1]), c=:jet, aspectratio=1))
