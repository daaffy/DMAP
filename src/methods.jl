using SparseArrays, ProgressMeter, NearestNeighbors, PyCall
# https://docs.julialang.org/en/v1/stdlib/SparseArrays/
# https://www.caam.rice.edu/software/ARPACK/

export dynamic_laplacian, single_dl, _kernel, sp_dynamic_laplacian
export npg_k, npg_ϵ

"""
    dynamic_laplacian

Dense calculation of the dynamic laplacian. Need to implement sparse version.
"""

function dynamic_laplacian(
    traj::Trajectory;
    metric = Euclidean(),
    eps_mult = 1.0,
    epsilon = eps_mult*nndist(traj.X[:,1,:])/sqrt(2), # according to heuristic; see summary
    threshold = 0.0,
    boundary_idxs = []
    )   
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)
    alpha = 1 # see Coifman/Lafon


    # p = Progress(n_t,desc="Calculating dynamic Laplacian...",barglyphs=BarGlyphs("[=> ]"),barlen=50,color=:yellow)

    f(x) = iso_kernel(x,0,epsilon)
    if (length(traj.t)>1)
        temp = traj.t[2:end]-traj.t[1:end-1]
        coeff = 0.5*[temp[1]; temp[1:end-1]+temp[2:end]; temp[end]]/ (traj.t[end]-traj.t[1])
    else
        # temp = 1.
        coeff = 1.
    end
    P_final = zeros(n_s,n_s)
    temp = similar(P_final)
    for ti = 1:n_t
        # next!(p)
        temp .= pairwise(metric,traj.X[:,ti,:],traj.X[:,ti,:],dims=1)
        # temp .= broadcast(f,temp)
        temp .= _kernel.(temp,epsilon)
        temp .*= broadcast(>=,temp,threshold) # threshold the kernel matrix using broadcasting
        
        # # diagnose
        # println(sum(P_final.>=1e-8)/(n_s^2))
        # display(heatmap(temp))

        q_bar = (sum(temp, dims=2)./n_s).^alpha
        temp ./= (q_bar .* q_bar')
        d_barsum = sum(temp, dims=1)./n_s
        temp ./= d_barsum
        P_final .+= coeff[ti].*temp./n_s
    end
    # P_final /= traj.t[end]-traj.t[1] # make sure this hasn't changed anything
    P_final

    if (boundary_idxs != [])
        for idx in boundary_idxs
            P_final[:,idx] = zeros(n_s)
            P_final[idx,:] = zeros(n_s)
        end
    end

    return P_final
end

function sp_dynamic_laplacian(
    traj::Trajectory;
    metric = Euclidean(),
    k = 20,
    # eps_mult = 1.0,
    epsilon = nndist(traj.X[:,1,:])/sqrt(2), # according to heuristic; see summary
    # threshold = 0.0,
    boundary_idxs = []
)
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)
    alpha = 1 # see Coifman/Lafon
    f(x) = iso_kernel(x,0,epsilon)
    T = traj.t[end] - traj.t[1]

    # tree = BallTree(traj.X[:,1,:]')
    # idxs = inrange(tree,traj.X[:,1,:]',0.1)
    # idxs_1 = vcat(idxs...)
    # idxs_2 = vcat([i*ones(Int32,length(idxs[i])) for i = 1:n_s]...)
    # # vals = similar(idxs_1)
    # vals = colwise(Euclidean(),traj.X[idxs_1,1,:]',traj.X[idxs_2,1,:]')
    # A = sparse(idxs_1,idxs_2,vals)
    
    idxs = Vector{Vector{Int64}}(undef,n_s) # can we further initialize the nested vector?
    idxs_2 = vcat([i*ones(Int32,k) for i = 1:n_s]...)
    idxs_1 = similar(idxs_2)
    vals_vec = zeros(Float64,length(idxs_1))
    vals = Vector{Vector{Int64}}(undef,n_s) # can we further initialize the nested vector?
    
    # print("here")
    P = spzeros(n_s,n_s)
    for t_i = 1:n_t
        temp = spzeros(n_s,n_s)
        
        tree = KDTree(traj.X[:,t_i,:]',metric,reorder=true)
        idxs,vals = knn(tree,traj.X[:,t_i,:]',k) # broadcast?
        idxs_1 .= vcat(idxs...)
        vals_vec .= vec(vcat(vals...))
        vals_vec .= _kernel.(vals_vec,epsilon)
        
        temp = sparse(idxs_1,idxs_2,vals_vec)
        # q_bar = (sum(temp, dims=2)./n_s).^alpha
        q_bar_inv = ((sum(temp, dims=2)./n_s).^alpha).^-1
        temp .*= q_bar_inv .* q_bar_inv'
        # d_barsum = sum(temp, dims=1)./n_s
        d_barsum_inv = (sum(temp, dims=1)./n_s).^-1
        temp .*= d_barsum_inv
        P .+= temp./n_s

        # temp = sparse(idxs_1,idxs_2,vals_vec)
        # q_bar = (sum(temp, dims=2)./n_s).^alpha
        # temp .*= (q_bar .* q_bar').^-1
        # d_barsum = sum(temp, dims=1)./n_s
        # temp ./= d_barsum
        # P .+= temp./n_s

    end

    # P ./= T
    return P


end

"""
may need this later for step evaluation?
"""
# function single_dl(traj::Trajectory;
#     metric = Euclidean(),
#     epsilon = nndist(traj.X[:,1,:])/sqrt(2), # according to heuristic; see summary
#     threshold = 0.0
#     )
#     @assert length(traj.t)==1

#     n_s = size(traj.X, 1)
#     alpha = 1

#     if (threshold != 0.0)
#         f(x) = (iso_kernel(x,0,epsilon) >= threshold ? iso_kernel(x,0,epsilon) : 0) # can we do it without writing iso_kernel twice?
#     else
#         f(x) = iso_kernel(x,0,epsilon)
#     end

#     dist_matrix = pairwise(metric,traj.X[:,1,:],traj.X[:,1,:],dims=1)
#     ker_matrix = broadcast(f,dist_matrix)
#     q_bar = sum(ker_matrix, dims=2)/n_s
#     # d_bar = 
#     d_bar = [ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha) for i = 1:n_s, j = 1:n_s]
#     d_barsum = sum(d_bar, dims=1)/n_s;
#     p_bar = [d_bar[i,j]/d_barsum[j] for i = 1:n_s, j = 1:n_s]
#     P_final += p_bar/n_s

#     return P_final
# end


"""
    npg
Padberg-Gehle and Schneide, "Network-based study of Lagrangian transport and mixing"; 2017.

"""
function npg_ϵ(
    traj::Trajectory;
    ϵ = nndist(traj.X[:,1,:]),
    metric = Euclidean()
)
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)

    # BallTree version has poor allocation strategy...
    A = spzeros(n_s,n_s)
    for t_i = 1:n_t
        tree = BallTree(traj.X[:,t_i,:]',metric)
        idxs = inrange(tree,traj.X[:,t_i,:]',ϵ)
        idxs_1 = vcat(idxs...)
        idxs_2 = vcat([i*ones(Int32,length(idxs[i])) for i = 1:n_s]...)
        
        A += sparse(idxs_1,idxs_2,1.0)
    end
    A = 1*(A.>0)
    A -= sparse(I,n_s,n_s)

    D = sparse((1:n_s),(1:n_s),vec(sum(A,dims=1)))

    return A, D
        
end

function npg_k(
    traj::Trajectory;
    k = 5,
    metric = Euclidean()
)   
    n_s = size(traj.X, 1)
    n_t = size(traj.X, 2)
    idxs_2 = vcat([i*ones(Int32,k) for i = 1:n_s]...)
    idxs_1 = similar(idxs_2)
    idxs = Vector{Vector{Int64}}(undef,n_s) # can we further initialize the nested vector?

    A = spzeros(n_s,n_s)
    for t_i = 1:n_t
        tree = KDTree(traj.X[:,t_i,:]',metric)
        idxs .= knn(tree,traj.X[:,t_i,:]',k)[1]
        idxs_1 .= vcat(idxs...)

        A += sparse(idxs_1,idxs_2,1) + sparse(idxs_2,idxs_1,1)
    end
    A += A'
    A = 1*(A.>0)
    A -= sparse(I,n_s,n_s)

    D = sparse((1:n_s),(1:n_s),vec(sum(A,dims=1)))
    
    return A, D

end








# function sparse_dynamic_laplacian(
#     traj::Trajectory;
#     metric = Euclidean(),
#     epsilon = nndist(traj.X[:,1,:])/sqrt(2), # according to heuristic; see summary
#     threshold = 0.0
#     )   
#     n_s = size(traj.X, 1)
#     n_t = size(traj.X, 2)
#     alpha = 1 # see Coifman/Lafon
#     f(x) = (iso_kernel(x,0,epsilon) >= threshold ? iso_kernel(x,0,epsilon) : 0) # can we do it without writing iso_kernel twice?

#     temp = traj.t[2:end]-traj.t[1:end-1]
#     coeff = 0.5*[temp[1]; temp[1:end-1]+temp[2:end]; temp[end]]
#     global P_final = spzeros(n_s,n_s)
#     for ti = 1:n_t
#         dist_matrix = pairwise(metric,traj.X[:,ti,:],traj.X[:,ti,:],dims=1)
#         ker_matrix = broadcast(f,dist_matrix)
#         ker_matrix = sparse(ker_matrix)
#         # ker_matrix = sparse(I(n_s))
#         # print(ker_matrix)
#         q_bar = sum(ker_matrix, dims=2)/n_s
#         d_bar = 
#         d_bar = [ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha) for i = 1:n_s, j = 1:n_s]
#         d_barsum = sum(d_bar, dims=1)/n_s;
#         p_bar = [d_bar[i,j]/d_barsum[j] for i = 1:n_s, j = 1:n_s]
#         global P_final += coeff[ti]*p_bar/n_s
#     end

#     P_final /= traj.t[end]-traj.t[1]

#     return P_final
# end

# function sparse_dynamic_laplacian(
#     traj::Trajectory,
#     metric = Euclidean(),
#     epsilon = nndist(traj.X[:,1,:])/sqrt(2) # according to heuristic; see summary
# )   
#     n_s = size(traj.X, 1)
#     n_t = size(traj.X, 2)
#     alpha = 1 # see Coifman/Lafon
#     f(x) = iso_kernel(x,0,epsilon)

#     temp = traj.t[2:end]-traj.t[1:end-1]
#     coeff = 0.5*[temp[1]; temp[1:end-1]+temp[2:end]; temp[end]]
#     global P_final = zeros(n_s,n_s)
#     for ti = 1:n_t
#         dist_matrix = pairwise(metric,traj.X[:,ti,:],traj.X[:,ti,:],dims=1)
#         ker_matrix = broadcast(f,dist_matrix)
#         q_bar = sum(ker_matrix, dims=2)/n_s
#         d_bar = [ker_matrix[i,j]/(q_bar[i]^alpha*q_bar[j]^alpha) for i = 1:n_s, j = 1:n_s]
#         d_barsum = sum(d_bar, dims=1)/n_s;
#         p_bar = [d_bar[i,j]/d_barsum[j] for i = 1:n_s, j = 1:n_s]
#         global P_final += coeff[ti]*p_bar/n_s
#     end

#     P_final /= traj.t[end]-traj.t[1]

#     return P_final
# end

function iso_kernel(
    x_1,
    x_2,
    epsilon,
    metric = Euclidean()
)   
    # thresh = 0
    # temp = exp(-metric(x_1,x_2)^2/epsilon)
    # return temp >= thresh ? temp : 0;
    return exp(-metric(x_1,x_2)^2/epsilon)
end

function _kernel(val,epsilon;)
    val = exp(-(val^2)/epsilon)
    # @.
end


# FTLE field function