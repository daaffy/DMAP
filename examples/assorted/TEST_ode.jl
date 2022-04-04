using DifferentialEquations
using Plots

# One dimensional example
# f(u,p,t) = 1.01*u # p = ?
# u0 = 1/2
# tspan = (0.0,1.0) 
# prob = ODEProblem(f,u0,tspan) # how do i choose the timesteps?
# sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8) 


# plot(sol)

function lorenz!(du,u,p,t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
prob = ODEProblem(lorenz!,u0,tspan)
sol = solve(prob)

plot(sol,vars=(1,2,3))