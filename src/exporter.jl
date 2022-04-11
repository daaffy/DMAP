using DataFrames, CSV

export export_eig

function export_eig(traj::Trajectory,vecs::AbstractArray{<:Real,2},directory_name::String)
    
    nt = size(traj.X,2)
    for ti = 1:nt
        if (size(traj.X,3)==2)
            df = DataFrame(
                xcoord=traj.X[:,ti,1],
                ycoord=traj.X[:,ti,2],
                eig1=vecs[:,1],
                eig2=vecs[:,2]
                )
        elseif (size(traj.X,3)==3)
            df = DataFrame(
                xcoord=traj.X[:,ti,1],
                ycoord=traj.X[:,ti,2],
                zcoord=traj.X[:,ti,3],
                eig1=vecs[:,1],
                eig2=vecs[:,2]
                )
        else
            error("not supported")
        end
        
        CSV.write(directory_name*"out"*string(ti)*".csv",df)
    end
    
end