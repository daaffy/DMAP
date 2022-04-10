# ------------------------------------------------------------------------------------------------------------------------
# interpolate.jl
# Interpolate a discrete vector field in both space and time
# ------------------------------------------------------------------------------------------------------------------------

using Distances, ScatteredInterpolation

export 
    rbf, 
    uinterp_0, 
    uinterp, # spatial interpolation specified by uinterp_0 (Gaussian), followed by linear interpolation in time
    test_uinterp_0, 
    calculate_uniform_bandwidth,
    evaluate_itp,
    example

# how do we include from the main DMAP.jl file:
# struct DiscreteVelocityField

function rbf(
    x,
    y,
    eps,
    metric = Euclidean()
)   
    # gaussian rbf centred at x evaluated at y
    return exp(-metric(x,y)^2/eps)
end

function uinterp(dvf::DiscreteVelocityField,mode::String)
    itps = [uinterp_0(dvf.X,dvf.dX[:,:,t_i],mode) for t_i = 1:length(dvf.t)]
    function eval(points,t)
        if (t <= dvf.t[1])
            return evaluate_itp(itps[1],points)
        elseif((t >= dvf.t[end]))
            return evaluate_itp(itps[end],points)
        else
            for ti = 1:length(dvf.t)
                if (dvf.t[ti] > t)
                    del1 = (dvf.t[ti]-t)/((dvf.t[ti]-dvf.t[ti-1]))
                    del2 = (t-dvf.t[ti-1])/((dvf.t[ti]-dvf.t[ti-1]))
                    return del1*evaluate_itp(itps[ti-1],points)+del2*evaluate_itp(itps[ti],points)
                    break
                end
            end
        end
    end
end

# function uinterp(coords,times,zvals,eps=1/sqrt(2))
#     # second edition unstructed interpolator
#     # zval is an array of zvalues indexed by time
#     n = size(coords,1)
#     ker_matrix = [iso_kernel(coords[i,:],coords[j,:],eps) for i = 1:n, j = 1:n]
#     coeff_vec = similar(zvals)
#     inv_ker = inv(ker_matrix)
#     for i = 1:length(times)
#         coeff_vec[:,:,i] = inv_ker*zvals[:,:,i]
#     end
#     function eval(X,t)
#         # interpolate through (simple linear interpolation)
#         if (t <= times[1])
#             coeff_vec = coeff_vec[:,:,1]
#         elseif(t >= times[length(times)])
#             coeff_vec = coeff_vec[:,:,length(times)]
#         else
#             for ti = 1:length(times)
#                 if (times[ti] > t)
#                     del1 = (times[ti]-t)/((times[ti]-times[ti-1]))
#                     del2 = (t-times[ti-1])/((times[ti]-times[ti-1]))
#                     coeff_vec = del1*coeff_vec[:,:,ti-1]+del2*coeff_vec[:,:,ti]
#                     break
#                 end
#             end
#         end

#         return sum(coeff_vec.*[iso_kernel(X',coords[i,:],eps) for i = 1:n],dims=1)[1]
#     end
# end


# to-do: add variable bandwidth kernels / choice over kernel mode
# can we use splat (...) here to fix some of the earlier issues???
function uinterp_0(points::AbstractArray{<:Real,2},samples::AbstractArray{<:Number,2},mode::String="ubk")
    if (cmp(mode,"ubk")==0) 
        # uniform bandwidth across all kernels
        eps = calculate_uniform_bandwidth(points)
        return uinterp_0(points,samples,Gaussian(1/sqrt(eps)))
    # elseif (cmp(mode,"vbk")==0)
    #     # variable bandwidth kernel
    #     eps = calculate_variable_bandwidth(coord_vec)
    elseif(cmp(mode,"nn")==0)
        return uinterp_0(points,samples,NearestNeighbor())
    else
        error("mode not recognised")
    end
    
end

function uinterp_0(points::AbstractArray{<:Real,2},samples::AbstractArray{<:Number,2},mode::ScatteredInterpolation.InterpolationMethod) # supplying eps as a scalar converts to a uniform kernel vector  
    # wrap ScatteredInterpolation.jl
    return itp = ScatteredInterpolation.interpolate(mode,points',samples) # ScatteredInterpolation accepts points in a different form to which I like
    # function eval(eval_points::AbstractArray{<:Real,2})
    #     return ScatteredInterpolation.evaluate(itp,eval_points)
    # end
end

function evaluate_itp(itp::ScatteredInterpolation.ScatteredInterpolant,points::AbstractVector{<:Real}) # single point input usually given as a vector; unsure if this is the most elegant solution
    return evaluate_itp(itp,vcat(points'))
end

function evaluate_itp(itp::ScatteredInterpolation.ScatteredInterpolant,points::AbstractArray{<:Real,2}) # type itp???
    # wrap ScatteredInterpolation.jl
    # println(ScatteredInterpolation.evaluate(itp,points'))
    return ScatteredInterpolation.evaluate(itp,points')
end

# function uinterp_0(coord_vec,zval,eps::Float64) # supplying eps as a scalar converts to a uniform kernel vector  
#     eps_vec = [eps for i = 1:size(coord_vec,1)]
#     return uinterp_0(coord_vec,zval,eps_vec)
# end

# function uinterp_0(
#     coord_vec,
#     zval, # restrict to Matrix or Vector{Vector}?
#     eps::Vector{Float64}
# )
#     # unstructured interpolation using radial basis functions with shape parameter = eps
#     # returns interpolation function eval

#     # checks
#     # ...

#     # uninterp_0 is a single time slice; can be used to debug the more general uinterp
#     n = size(coord_vec,1)
#     # eps[i] or eps[j] below?
#     ker_matrix = [rbf(coord_vec[i,:],coord_vec[j,:],eps[j]) for i = 1:n, j = 1:n]
#     # coeff_vec = ker_matrix\zval
#     coeff_vec = inv(ker_matrix)*zval
#     function eval(U)
#         temp = sum(coeff_vec.*[rbf(U,coord_vec[i,:],eps[i]) for i = 1:n],dims=1)
#         return length(temp) == 1 ? sum(temp) : temp
#         # returns datatype: array of vectors?
#     end
# end

# function uninterp_0(
#     coord_vec,
#     zval, # restrict to Matrix or Vector{Vector}?
#     eps::Vector{Float64},
#     dist_matrix::Matrix{Float64}
# )
#     # what if dist_matrix has already been pre-calculated, e.g. in vbf bandwidth calculation?
#     # this function allows dist_matrix to be supplied

# end

function test_uinterp_0(f_interp,X::Matrix{Float64},z::Matrix{Float64})
    diff = [(f_interp(X[i,:])-z[i,:]) for i = 1:length(z)]
    # err = sum(diff.^2) # error here
    return err = maximum(broadcast(abs,diff))
end

function test_uinterp_0(f_interp,X::Matrix{Float64},z::Vector{Vector{Float64}})
    diff = vcat([(f_interp(X[i,:])-z[i])' for i = 1:length(z)]...)
    # err = sum(diff.^2) # error here
    return err = maximum(broadcast(abs,diff))
end

function calculate_uniform_bandwidth(X,factor=0.815)
    n = size(X,1)
    dist_matrix = [norm(X[i,:]-X[j,:]) for i = 1:n, j = 1:n]
    d_vec = sort(dist_matrix,dims=2)[:,2]
    return factor*sum(d_vec)/n # Hardy's, obtained from S. Rippa paper.
end

function calculate_variable_bandwidth(X,factor=0.815)
    n = size(X,1)
    dist_matrix = [norm(X[i,:]-X[j,:]) for i = 1:n, j = 1:n]
    d_vec = sort(dist_matrix,dims=2)[:,2]
    return factor*d_vec
end