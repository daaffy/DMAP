using DifferentialEquations

export create_grid, solve_trajectory, test

struct Trajectory
    X::Array{Float64,3}
    t::Vector{Float64}
end

struct ContinuousVelocityField
    # to be fed into solve_trajectory
end

# to-do:
# - creating grids in more complex geometries e.g., in simple pipe
# - cut-off points that leave the geometry
# ** use an indicator function of some kind to solve the above problems?

function indicator()

end

function create_grid(args...)
    # create a uniform grid in the appropriate representation to be fed into the trajectory solver
    
    # checks 
    # ...
    
    d = length(args)-1 # spatial dimensions
    lims = [args[d_i][1]:args[d+1][d_i]:args[d_i][2] for d_i = 1:d]
    temp = [collect(i)' for i in Iterators.product(lims...)] # the transpose (') is the only way I could figure out to convert to an appropriate type
    return vcat(temp...)
end

function solve_trajectory(
    vec_field, # make this an abstract type
    X_0::Matrix{Float64},
    t::Vector{Float64}
)
    # return Trajectory (cut-off values outside indicator function, if provided)
    n_t = length(t)
    n_s = size(X_0,1)
    d = size(X_0,2)

    X = Array{Float64}(undef,n_s,n_t,d) # not sure if I like the representation order; (n_s,d,n_t) more consistent?
    for i = 1:n_s
        local u0 = [X_0[i,1];X_0[i,2]]
        local prob = ODEProblem(vec_field,u0,(t[1],t[end]))
        local sol = solve(prob,saveat=t)
        for j = 1:n_t
            X[i,j,:] = sol.u[j]
        end
    end

    return Trajectory(X,t)
end

function solve_trajectory(
    vec_field,
    X::Matrix{Float64},
    t::AbstractRange
)   
    return solve_trajectory(vec_field,X,collect(t))
end