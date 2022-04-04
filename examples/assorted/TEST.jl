# using Revise
# using CoherentStructures

# r = 6
# V = real(v[:,end:-1:end-r+1]) 
# S, R = SEBA(V)

layer = 4
mask = zeros(n_s,nx*ny)
inds = vec(map2[:,:,layer])
for j = 1:length(inds)
    mask[Int(inds[j]),j] = 1
end
mask = mask'

# display(scatter(mask*X_traj[:,1,1], mask*X_traj[:,1,2], zcolor=mask*real(v[:,end-2]), c=:jet, aspectratio=1))
display(scatter(mask*X_traj[:,1,1], mask*X_traj[:,1,2], zcolor=mask*real(S[:,6]), c=:jet, aspectratio=1))