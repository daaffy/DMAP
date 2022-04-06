using CSV
using DataFrames

export load_dvf

function load_dvf(input_files,input_dir)
    # assists in loading data from a ".csv time-series" to a DiscreteVelocityField in Julia
    # ASSUMPTION: velocity field is evaluated on a constant grid
    discrete_steps = length(input_files)

    # this can surely be cleaned up...
    i = 1
    df = CSV.read(input_dir * input_files[i], DataFrame)
    df = Array(df)
    n_points = size(df,1)
    #

    X = Array{Float64}(undef,n_points,3)
    for j = 1:n_points
        global X[j,:] = df[j,2:4]'
    end

    dX = Array{Float64}(undef,n_points,3,discrete_steps)
    for i = 1:discrete_steps
        global df = CSV.read(input_dir * input_files[i], DataFrame)
        global df = Matrix(df)

        for j = 1:n_points
            # check that grid doesn't change through time (within margin; use initial grid anyway)
            # if (X[j,:] != df[j,2:4]')
            #     print(-df[j,2:4]')
            #     return
            # end
            global dX[j,:,i] = df[j,5:7]'
        end
    end

    # for this example we need to manually construct the time vector since the data is not supplied
    t = LinRange(0,1,discrete_steps)

    return DiscreteVelocityField(X,dX,t)
end