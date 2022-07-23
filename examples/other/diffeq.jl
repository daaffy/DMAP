# following https://tutorials.sciml.ai/html/introduction/03-optimizing_diffeq_code.html

using DifferentialEquations, BenchmarkTools, Plots

# function lorenz(u,p,t)
#     dx = 10.0*(u[2]-u[1])
#     dy = u[1]*(28.0-u[3]) - u[2]
#     dz = u[1]*u[2] - (8/3)*u[3]
#     return [dx,dy,dz]
# end

# u0 = [1.0;0.0;0.0]
# tspan = (0.0,100.0)
# prob = ODEProblem(lorenz,u0,tspan)
# @benchmark solve(prob,Tsit5())

# ---
function lorenz!(du,u,p,t)
    du[1] = 10*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

function vec_field!(du,u,p,t)
    du[1] = 1
    du[2] = p
    du[3] = 0
end

# testing callbacks
function _affect!(integrator,val)
    integrator.p = val
end

# function affect!(integrator)
#     integrator.p = -2
# end

vals=[-0.0;2.0]
tstop = [0.00;0.75]
cbs = CallbackSet([
    DiscreteCallback(
        (u,t,integrator) -> t in [tstop[i]], # condition
        # affect!
        integrator -> _affect!(integrator,vals[i])
    )   
for i = 1:2]...)

u0 = [0.0;0.0;0.0]
tspan = [0.0;1.0]
prob = ODEProblem(vec_field!,u0,tspan,10)
sol = solve(prob,Tsit5(),reltol=1e-6,callback = cbs, tstops=tstop)

u = hcat(sol.u...)

plot(u[1,:],u[2,:])
