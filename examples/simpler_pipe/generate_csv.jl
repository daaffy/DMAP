using DataFrames, CSV

# X = Array{Float64}(undef,0,3)
# dX = Array{Float64}(undef,0,3)

# for r = 0:0.1:1, phi = LinRange(0,2*pi,10)[1:end-1], z = 0:0.5:10
    
# end

temp = [[r*cos(phi) r*sin(phi) z] for r = 0:0.1:1, phi = LinRange(0,2*pi,20)[1:end-1], z = 0:0.2:10]
X = vcat(temp...)
dX = vcat([[0. 0. 0.] for i = 1:size(X,1)]...)
nodes = collect(1:size(X,1))
df = DataFrame(
                nodenum=nodes,
                xcoord=X[:,1],
                ycoord=X[:,2],
                zcoord=X[:,3],
                dx=dX[:,1],
                dy=dX[:,2],
                dz=dX[:,3],
                )

CSV.write("./examples/simpler_pipe/input/time_0_2.csv",df)