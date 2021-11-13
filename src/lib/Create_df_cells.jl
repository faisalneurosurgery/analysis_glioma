module Create_df_cells

using Revise
using PyPlot
using FileIO, DataFrames, ProgressMeter, Statistics, LinearAlgebra
# perso
using Winding_number

export create_df_velocity, create_df_velocity_zone, create_df_t_xy_uv_zone,
    create_df_correlation, create_df_correlation_zone

""" create the dataframe with all the speed, angle, u,v """
function create_df_velocity(all_path,xMin,xMax,yMin,yMax)
    # init
    N = length(all_path)            # number of tracks
    # create dataframe
    #-----------------
    df = DataFrame(u = Float64[], v = Float64[], speed = Float64[], θ = Float64[])
    @showprogress 1 "Computing..."  for i=1:N
        path_i = all_path["path"*string(i)]
        n_i = nrow(path_i)
        for k = 1:n_i
            # init
            x,y = path_i[k,:x], path_i[k,:y]
            # save if in box
            if (x>xMin) & (x<xMax) & (y>yMin) & (y<yMax)
                u,v = path_i[k,:u], path_i[k,:v]
                θ = atan(v,u)
                speed = sqrt(u^2+v^2)
                push!(df, [u v speed atan(v,u)])
            end
        end
    end
    return df
end

""" Create a dataframe with all the velocity inside a zone. """
function create_df_velocity_zone(all_path,P,xMin,xMax,yMin,yMax)
    # init
    N = length(all_path)            # number of tracks
    # create dataframe
    #-----------------
    df = DataFrame(u = Float64[], v = Float64[], speed = Float64[], θ = Float64[])
    @showprogress 1 "Computing..."  for i=1:N
        path_i = all_path["path"*string(i)]
        n_i = nrow(path_i)
        for k = 1:n_i
            # init
            x,y = path_i[k,:x], path_i[k,:y]
            ## temporal
            if (path_i[k,:time] < 160*598.542)
            # save if in box
            if (x>xMin) & (x<xMax) & (y>yMin) & (y<yMax)
                # in the box
                if (abs(wind_number([x,y],P))>.9)
                    # i in the zone
                    u,v = path_i[k,:u], path_i[k,:v]
                    θ = atan(v,u)
                    speed = sqrt(u^2+v^2)
                    push!(df, [u v speed atan(v,u)])
                end
            end
            end #temporal
        end
    end
    return df
end


""" Create a dataframe with time/position/velocity inside a zone. """
function create_df_t_xy_uv_zone(all_path,P,xMin,xMax,yMin,yMax)
    # init
    N = length(all_path)            # number of tracks
    # create dataframe
    #-----------------
    df = DataFrame(idx_cell = Int64[], t = Float64[], x = Float64[], y = Float64[], u = Float64[], v = Float64[], speed = Float64[], θ = Float64[])
    @showprogress 1 "Computing..."  for i=1:N
        path_i = all_path["path"*string(i)]
        n_i = nrow(path_i)
        for k = 1:n_i
            # init
            t,x,y = path_i[k,:time], path_i[k,:x], path_i[k,:y]
            # save if in box
            if (x>xMin) & (x<xMax) & (y>yMin) & (y<yMax)
                # in the box
                if (abs(wind_number([x,y],P))>.9)
                    # i in the zone
                    u,v = path_i[k,:u], path_i[k,:v]
                    θ = atan(v,u)
                    speed = sqrt(u^2+v^2)
                    push!(df, [i t x y u v speed atan(v,u)])
                end
            end
        end
    end
    return df
end


""" Create a dataframe with all the relative position and correlation. """
function create_df_correlation(all_path,xMin,xMax,yMin,yMax,rMax)
    # init
    N = length(all_path)         # number of tracks
    # create dataframe
    #-----------------
    df = DataFrame(r = Float64[], theCor = Float64[], x_ij = Float64[], y_ij = Float64[], x_ji = Float64[], y_ji = Float64[])
    @showprogress 1 "Computing..."  for i=1:(N-1)
        path_i = all_path["path"*string(i)]
        for j=(i+1):N
            # init
            path_j = all_path["path"*string(j)]
            path_ij = innerjoin(path_i, path_j, on = :time, makeunique=true)
            #path_ij = join(path_i, path_j, on = :time)
            n_ij = nrow(path_ij)
            for k = 1:n_ij
                # init
                X_i = [path_ij[k,:x]; path_ij[k,:y]]
                ω_i = [path_ij[k,:u]; path_ij[k,:v]]
                X_j = [path_ij[k,:x_1]; path_ij[k,:y_1]]
                ω_j = [path_ij[k,:u_1]; path_ij[k,:v_1]]
                # normalize
                ω_i ./= norm(ω_i)
                ω_j ./= norm(ω_j)
                if (min(X_i[1],X_j[1])>xMin) & (max(X_i[1],X_j[1])<xMax) & (min(X_i[2],X_j[2])>yMin) & (max(X_i[2],X_j[2])<yMax)
                    # stat
                    r = norm(X_j-X_i)
                    if (r<rMax)
                        theCor = ω_i⋅ω_j
                        x_ij = (X_j-X_i)⋅ω_i
                        y_ij = (X_j-X_i)⋅[-ω_i[2];ω_i[1]]
                        x_ji = (X_i-X_j)⋅ω_j
                        y_ji = (X_i-X_j)⋅[-ω_j[2];ω_j[1]]
                        # saving
                        push!(df,  [r theCor x_ij y_ij x_ji y_ji])
                    end
                end
            end
        end
    end
    return df
end

""" Create a dataframe with all the relative position and correlation inside a zone """
function create_df_correlation_zone(all_path,P,xMin,xMax,yMin,yMax,rMax)
    # init
    N = length(all_path)         # number of tracks
    # create dataframe
    #-----------------
    df = DataFrame(r = Float64[], theCor = Float64[], x_ij = Float64[], y_ij = Float64[], x_ji = Float64[], y_ji = Float64[])
    @showprogress 1 "Computing..."  for i=1:(N-1)
        path_i = all_path["path"*string(i)]
        for j=(i+1):N
            # init
            path_j = all_path["path"*string(j)]
            path_ij = innerjoin(path_i, path_j, on = :time, makeunique=true)
            #path_ij = join(path_i, path_j, on = :time)
            n_ij = nrow(path_ij)
            for k = 1:n_ij
                # init
                X_i = [path_ij[k,:x]; path_ij[k,:y]]
                ω_i = [path_ij[k,:u]; path_ij[k,:v]]
                X_j = [path_ij[k,:x_1]; path_ij[k,:y_1]]
                ω_j = [path_ij[k,:u_1]; path_ij[k,:v_1]]
                # normalize
                ω_i ./= norm(ω_i)
                ω_j ./= norm(ω_j)
                if (min(X_i[1],X_j[1])>xMin) & (max(X_i[1],X_j[1])<xMax) & (min(X_i[2],X_j[2])>yMin) & (max(X_i[2],X_j[2])<yMax)
                    # stat
                    r = norm(X_j-X_i)
                    if (r<rMax)
                        if (abs(wind_number(X_i,P))>.9)
                            if (abs(wind_number(X_j,P))>.9)
                                # both in box and in zone
                                theCor = ω_i⋅ω_j
                                x_ij = (X_j-X_i)⋅ω_i
                                y_ij = (X_j-X_i)⋅[-ω_i[2];ω_i[1]]
                                x_ji = (X_i-X_j)⋅ω_j
                                y_ji = (X_i-X_j)⋅[-ω_j[2];ω_j[1]]
                                # saving
                                push!(df,  [r theCor x_ij y_ij x_ji y_ji])
                            end
                        end
                    end
                end
            end
        end
    end
    return df
end

end                             # module
