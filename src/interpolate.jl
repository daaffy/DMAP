# ------------------------------------------------------------------------------------------------------------------------
# interpolate.jl
# Interpolate a discrete vector field in both space and time
# ------------------------------------------------------------------------------------------------------------------------


export rbf, uinterp_0, uinterp, test_uinterp_0, calculate_eps

function rbf(
    x,
    y,
    eps,
    metric = Euclidean()
)   
    # gaussian rbf centred at x evaluated at y
    return exp(-metric(x,y)^2/eps)
end

function uinterp(coords,times,zvals,eps=1/sqrt(2))
    # second edition unstructed interpolator
    # zval is an array of zvalues indexed by time
    n = size(coords,1)
    ker_matrix = [iso_kernel(coords[i,:],coords[j,:],eps) for i = 1:n, j = 1:n]
    coeff_vec = similar(zvals)
    inv_ker = inv(ker_matrix)
    for i = 1:length(times)
        coeff_vec[:,:,i] = inv_ker*zvals[:,:,i]
    end
    function eval(X,t)
        # interpolate through (simple linear interpolation)
        if (t <= times[1])
            coeff_vec = coeff_vec[:,:,1]
        elseif(t >= times[length(times)])
            coeff_vec = coeff_vec[:,:,length(times)]
        else
            for ti = 1:length(times)
                if (times[ti] > t)
                    del1 = (times[ti]-t)/((times[ti]-times[ti-1]))
                    del2 = (t-times[ti-1])/((times[ti]-times[ti-1]))
                    coeff_vec = del1*coeff_vec[:,:,ti-1]+del2*coeff_vec[:,:,ti]
                    break
                end
            end
        end

        return sum(coeff_vec.*[iso_kernel(X',coords[i,:],eps) for i = 1:n],dims=1)[1]
    end
end


# to-do: add variable bandwidth kernels / choice over kernel mode
# can we use splat (...) here to fix some of the earlier issues???
function uinterp_0(
    coord_vec,
    zval,
    eps=1/sqrt(2)
)
    # unstructured interpolation using radial basis functions with shape parameter = eps
    # returns interpolation function eval

    # uninterp_0 is a single time slice; can be used to debug the more general uinterp
    n = size(coord_vec,1)
    ker_matrix = [rbf(coord_vec[i,:],coord_vec[j,:],eps) for i = 1:n, j = 1:n]
    # coeff_vec = ker_matrix\zval
    coeff_vec = inv(ker_matrix)*zval
    function eval(U)
        temp = sum(coeff_vec.*[rbf(U,coord_vec[i,:],eps) for i = 1:n],dims=1)
        return length(temp) == 1 ? sum(temp) : temp
        # returns datatype: array of vectors?
    end
end

function test_uinterp_0(f_interp,X,z)
    diff = [f_interp(X[i,:])-z[i] for i = 1:length(z)]
    err = sum(diff.^2) # error here
    print(err)
end

function calculate_eps(X,factor=0.815)
    n = size(X,1)
    dist_matrix = [norm(X[i,:]-X[j,:]) for i = 1:n, j = 1:n]
    d_vec = sort(dist_matrix,dims=2)[:,2]
    return factor*sum(d_vec)/n # Hardy's, obtained from S. Rippa paper.
end