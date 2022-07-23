using Revise
using DMAP
using LinearAlgebra
using Plots
using Arpack
# using CSV
# using DataFrames
using DifferentialEquations, DifferentialEquations.EnsembleAnalysis

# -------------------------------------------------------------------------------------------------------
# extract discrete velocity field to X, dX multi-dimensional arrays.

# input_dir = "./examples/simple_pipe/input/"
# input_files = readdir(input_dir)

# dvf = load_dvf(input_files,input_dir)

# -------------------------------------------------------------------------------------------------------
# build interpolated velocity field

# @time cvf = uinterp(dvf,"shep")

# A = rand(3000,3000)
# @time eigen(A)

# function f(x,t)
#     return [0;1;0]
# end
# cvf = f

# -------------------------------------------------------------------------------------------------------
# integrate trajectories

"""
    let's sample an well known analytic example to test out the interpolation method
"""

# rotating double gyre velocity field
# function vec_field(u,t)
#     temp = (t.^2).*(3-2*t)
#     return [pi*sin(2*pi*u[1])*cos(pi*u[2]) 2*pi*sin(pi*u[1])*cos(2*pi*u[2]); -2*pi*cos(2*pi*u[1])*sin(pi*u[2]) -pi*cos(pi*u[1])*sin(2*pi*u[2])][1-temp; temp];
# end





input_dir = "./examples/simpler_pipe/input/"
input_files = readdir(input_dir)

dvf = load_dvf(input_files,input_dir)

# cvf = uinterp(dvf,"shep")

X_0 = dvf.X[1:10,:]
t = LinRange(0,1,20)


# set up vec_field!
using ScatteredInterpolation

interval = 1
itp1 = ScatteredInterpolation.interpolate(Shepard(),dvf.X',dvf.dX[:,:,interval])
itp2 = ScatteredInterpolation.interpolate(Shepard(),dvf.X',dvf.dX[:,:,interval+1])
# ScatteredInterpolation.evaluate(itp,[0;0;0])

function itp_obj(t,t_interval,dvf_interval,X)
    t_prime = (t-t_interval[1])/(t_interval[2]-t_interval[1]) # reparameterise t; 0 <= t_prime <= 1
    dvf_ti = t_prime*dvf_interval[2]+(1-t_prime)*dvf_interval[1]
    return ScatteredInterpolation.interpolate(Shepard(),X',dvf_ti)
end

p = t -> itp_obj(dvf.t[2]*0.5,dvf.t[1:2],[dvf.dX[:,:,1], dvf.dX[:,:,2]],dvf.X)
function vec_field!(du,u,p,t)
    # itp = itp_obj(t,dvf.t[i:i+1],dXs,dvf.X)
    temp = vec(ScatteredInterpolation.evaluate(p(t),u))
    du[1] = temp[1]
    du[2] = temp[2]
    du[3] = temp[3]
end

# itp = itp_obj(dvf.t[2]*0.5,dvf.t[1:2],[dvf.dX[:,:,1], dvf.dX[:,:,2]],dvf.X)

# i = 1
# p = t -> itp_obj(t,dvf.t[i:i+1],[dvf.dX[:,:,i], dvf.dX[:,:,i+1]],dvf.X)

# p(t) = ScatteredInterpolation.interpolate(Shepard(),dvf.X',dvf.dX[:,:,1])

# function _affect!(integrator,i)
#     # return t -> itp_obj(t,dvf.t[i:i+1],[dvf.dX[:,:,i], dvf.dX[:,:,i+1]],dvf.X)
#     integrator.p = t -> itp_obj(t,dvf.t[i:i+1],[dvf.dX[:,:,i], dvf.dX[:,:,i+1]],dvf.X)
# end

# function affect!(integrator)
#     integrator.p = t -> itp_obj(t,dvf.t[2:2+1],[dvf.dX[:,:,2], dvf.dX[:,:,2+1]],dvf.X)
# end


interval = [1]
dXs = [dvf.dX[:,:,1], dvf.dX[:,:,2]] # initialise
function temp_affect!(integrator, interval::Array{Int,1}, dXs)
    # i = integrator.p += 1

    # interval[1] += 1
    # dXs .= [dvf.dX[:,:,i], dvf.dX[:,:,i+1]]
    # print("hello")
end

# differential equation stuff
prob = ODEProblem(vec_field!,X_0[1,:],(t[1],t[end]),p)
initial_conditions = X_0'
function prob_func(prob,i,repeat)
    remake(prob,u0=initial_conditions[:,i])
end
ensemble_prob = EnsembleProblem(prob; prob_func)

# (!) if there are enough dvf slices to warrant callbacks
cbs = CallbackSet([DiscreteCallback(
    (u,t,integrator)->t in [dvf.t[i]], # condition
    # integrator -> _affect!(integrator,i) # affect!
    integrator -> temp_affect!(integrator,interval,dXs)
    )
    for i = 2:2]...)

tstop = dvf.t[2:2]

sim = solve(ensemble_prob,Tsit5(),EnsembleSerial(),reltol=1e-6,trajectories=size(initial_conditions,2),saveat=t,callback = cbs, tstops=tstop) # ,callback = cbs, tstops=tstop

# construct X
n_t = length(t)
n_s = size(X_0,1)
d = size(X_0,2)
X = Array{Float64}(undef,n_s,n_t,d)
for ti = 1:length(t)
    X[:,ti,:] = hcat(componentwise_vectors_timestep(sim,ti)...)
end

traj = Trajectory(X,t)

# traj.X

###### let's get here:
# @time traj = solve_ensemble(vec_field!,X_0,t,p)



# @time traj = solve_trajectory(cvf,X_0,t)

# indicator(x) = x[2] < 0.04 && x[2] > 0.01 ? true : false
# traj = clip(traj,indicator) # throw away trajectories that leave the domain of interest; defined by indicator()

plot(traj.t,traj.X[2,:,3])

















# # -------------------------------------------------------------------------------------------------------
# # process coherent structures

# @time P = dynamic_laplacian(traj, threshold=0.90)
# r = 2
# @time λ, v = eigs(sparse(P'), nev=r+1, which=:LM) # sparse method; transpose converts sparse matrix back into a dense matrix?
# # @time λ, v = eigen(P') # transpose?


# # V = real(v[:,end:-1:end-r+1]) 
# V = real(v[:,1:r])
# @time S, R = SEBA(V)

# # export_eig(traj,S,"./examples/simple_pipe/output/") # export

# # -------------------------------------------------------------------------------------------------------
# # test plots

# # display(scatter(traj.X[:,:,1]',traj.X[:,:,2]'))

# display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(S[:,1]), c=:jet, aspectratio=1))

# # display(scatter(traj.X[:,1,1], traj.X[:,1,2], zcolor=real(v[:,end-1]), c=:jet, aspectratio=1))
# # display(scatter(X_0[:,1], X_0[:,2], zcolor=real(S[:,2]), c=:jet, aspectratio=1))

# -------------------------------------------------------------------------------------------------------
# junk

# # test
# i = 1
# j = 7
# eval_points = dvf.X[i,:]
# println(dvf.dX[i,:,j])
# println(cvf(eval_points,dvf.t[j]))
