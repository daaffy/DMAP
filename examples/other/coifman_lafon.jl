using Distributions
using Plots

function generate_cluster(centre; r = 1, n = 10)
    C = [r 0; 0 r]
    d = MvNormal(centre, C)
    return rand(d, n)
end

# generate 3 clusters
points = [generate_cluster([1,0], r=0.3, n=300) generate_cluster([-3,0], r=0.3, n=300) generate_cluster([-2,2], r=0.3, n=300)]
# scatter(points[1,:],points[2,:])
