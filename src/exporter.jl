using DataFrames, CSV

export export_eig, export_traj, exportf

function exportf(traj::Trajectory,vecs::AbstractArray{},directory_name::String)
    
    nt = size(traj.X,2)
    for ti = 1:nt
        if (size(traj.X,3)==2)
            df = DataFrame(
                xcoord=traj.X[:,ti,1],
                ycoord=traj.X[:,ti,2],
                # eig1=vecs[:,1],
                # eig2=vecs[:,2]
                )
        elseif (size(traj.X,3)==3)
            df = DataFrame(
                xcoord=traj.X[:,ti,1],
                ycoord=traj.X[:,ti,2],
                zcoord=traj.X[:,ti,3],
                )
        else
            error("not supported")
        end
        for eig_index = 1:size(vecs,2)
            df[!,"eig"*string(eig_index)] = vecs[:,eig_index]
        end
        CSV.write(directory_name*"out"*string(ti)*".csv",df)
    end
    
end


"""
    export_traj

Somewhat thrown together at this stage in order to debug trajectory data for large data sets.

"""
function export_traj(traj::Trajectory,directory_name::String)
    nt = size(traj.X,2)
    for ti = 1:nt
        df = DataFrame(
                xcoord=traj.X[:,ti,1],
                ycoord=traj.X[:,ti,2],
                zcoord=traj.X[:,ti,3]
                )
        CSV.write(directory_name*"/traj_step_"*string(ti)*".csv",df)
    end
end


"""
    exportf

! Would be better to exportf a variable collection of arrays,
  and have optional coordinate/function format. 

"""
function exportf(
    X::Array{<:Real},
    f::AbstractArray{},
    export_path::String
)
    @assert size(X,1) == size(f,1)

    df = DataFrame()
    for coord_index = 1:size(X,2)
        df[!,"coord"*string(coord_index)] = X[:,coord_index]
    end     

    for f_index = 1:size(f,2)
        df[!,"f"*string(f_index)] = f[:,f_index]
    end
    
    CSV.write(export_path,df)
end
