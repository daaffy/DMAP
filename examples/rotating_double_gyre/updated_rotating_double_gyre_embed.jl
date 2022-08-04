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

s = t -> (t^2)*(3-2*t)
function vec_field!(du,u,p,t)
    x,y = u
    mul!(du,[pi*sin(2*pi*x)*cos(pi*y) 2*pi*sin(pi*x)*cos(2*pi*y); -2*pi*cos(2*pi*x)*sin(pi*y) -pi*cos(pi*x)*sin(2*pi*y)],[1-p(t); p(t)])
end

# -------------------------------------------------------------------------------------------------------
# integrate to determine trajectories 
# note: combining integration with dynamic Laplacian construction will save time if velocity vector field is given as input

res = 80
X_0 = create_grid([0 1],[0 1],[res res])
t = (0:0.0125:1)

@time traj = solve_ensemble(vec_field!,X_0,t,s)

id = 1 # 900?
function sgn(x)
    if (x >= 0)
        return 1
    else
        return -1
    end
end

for i = 1:length(t)
    traj_curr = Trajectory( traj.X[:,1:i,:],traj.t[1:i])
    @time P = dynamic_laplacian(traj_curr, threshold=0.)
    
    
    
    #
    # -------------------------------------------------------------------------------------------------------
    # calculate eigens 
    r = 10
    @time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?

    # -------------------------------------------------------------------------------------------------------
    # SEBA post-processing

    V = real(v[:,1:8])
    S, R = SEBA(V[:,1:2])

    if (i == 1)
        global sgn_2 = sgn(V[id,2])
        global sgn_3 = sgn(V[id,3])
        global sgn_4 = sgn(V[id,4])
    end

    if (sgn(V[id,2]) != sgn_2)
        V[:,2] *= -1
    end

    if (sgn(V[id,3]) != sgn_3)
        V[:,3] *= -1
    end

    if (sgn(V[id,4]) != sgn_4)
        V[:,4] *= -1
    end
    # S = V

    # export_eig(traj,S,"./examples/rotating_double_gyre/output/") # export

    # -------------------------------------------------------------------------------------------------------
    # # plot
    # display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(V[:,1]), c=:jet, aspectratio=1))
    # display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))

    # diffusion map embedding
    display(scatter( real(V[:,2]), real(V[:,3]), real(V[:,4]), zcolor=real(S[:,1]), camera = (45,45) ))
    
    # display(scatter( real(V[:,2]), real(V[:,3]), zcolor=real(V[:,1]), camera = (45,45) ))

    exportf(traj.X[:,1,:],[V S traj.X[:,end,:]],"./examples/rotating_double_gyre/embedding_output_02/out_"*string(i)*".csv")
end