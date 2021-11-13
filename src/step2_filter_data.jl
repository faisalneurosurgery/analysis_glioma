#=
Smooth the cell's trajectories
=#

using ProgressMeter
using CSV, JSON, DataFrames, FileIO
using DSP, Statistics

# parameter
list_nameExp = ["NPA_stich_3"]
folderData = "../data_tracking/step1_data_trackMate/"
folderSave = "../data_tracking/step2_data_filtered/"
l = 4                           # filter
σ = 2

function smooth_path(t,x,y,theFilter)
    """ smooth the data (t,x,y) using a filter and return the velocity
    assuming the interval time  is uniform
    """
    # init
    bigL = length(theFilter)
    l = floor(Int64,(bigL-1)/2)    # l needs to be odd
    Δt = t[2]-t[1]
    # smooth position
    x_smooth = DSP.conv(x,theFilter)[bigL:(end-bigL+1)] # need to 'cut' beginning and ending
    y_smooth = DSP.conv(y,theFilter)[bigL:(end-bigL+1)]
    # compute velocity
    u = (x_smooth[3:end] - x_smooth[1:(end-2)])/(2*Δt)*3600    # in μm/h
    v = (y_smooth[3:end] - y_smooth[1:(end-2)])/(2*Δt)*3600

    return [t[(l+2):(end-(l+1))],
            x_smooth[2:(end-1)],
            y_smooth[2:(end-1)],
            u,
            v]
end


#---------------------------------------------------#
#                     Mega loop                     #
#---------------------------------------------------#
for nameExp in list_nameExp
    println(nameExp)
    #-----------------------#
    # A) Initialization     #
    #-----------------------#
    # A.1) load data
    spots = CSV.read(folderData*nameExp*"/"*nameExp*"_spots.csv",DataFrame)
    tracks = CSV.read(folderData*nameExp*"/"*nameExp*"_tracks.csv",DataFrame)
    N = length(tracks[:,:Label])         # number of tracks
    global infoFile = Dict()
    open(folderData*nameExp*"/"*nameExp*"_info.json", "r") do f
        global infoFile
        infoFile = JSON.parse(f)
    end
    # A.2) filter
    theFilter = exp.(-(-l:l).^2/(2*σ))
    theFilter /= sum(theFilter)
    # A.3) saving
    Δt = infoFile["Frame interval (s)"]
    dico_path = Dict()
    global count_path = 0

    #--------------------------------------------------------------#
    #----------          B) the loop              -----------------#
    #--------------------------------------------------------------#
    @showprogress 1 "Computing..." for k=1:N
        if (tracks[k,:NUMBER_SPOTS]>10)  # needs at least 10 points
            # load
            index_track = findall( (spots[:,:TRACK_ID] .== tracks[k,:TRACK_ID]))
            x = convert(Vector{Float64}, spots[index_track,:POSITION_X])        # position (μm)
            y = convert(Vector{Float64}, spots[index_track,:POSITION_Y])
            t = convert(Vector{Float64}, spots[index_track,:FRAME])*Δt # or spots[:POSITION_T] for NPA data
            # smooth
            tS,xS,yS,u,v=smooth_path(t,x,y,theFilter)
            # check if path is not at the boundary
            if (mean(abs.(u))>1e-10 && mean(abs.(v))>1e-10)
                # save
                global count_path += 1
                df_path = hcat( DataFrame([round.(Int,tS/Δt)],[:frame]), # int
                                DataFrame([tS xS yS u v],[:time,:x,:y,:u,:v]) ) # all float
                #df_path = convert(DataFrame, [round.(Int,tS/Δt) tS xS yS u v]) # wrong type for the 'frame'
                #names!(df_path,map(parse,["frame","time","x","y","u","v"]))
                dico_path["path"*string(count_path)] = df_path
            end
        end
    end
    #--------------------------------------------------------------#
    #--------------------------------------------------------------#
    #-----------------------#
    # C) Saving             #
    #-----------------------#
    folderName = folderSave*nameExp
    if (!isdir(folderName))
        # check if the folder does not exist
        mkpath(folderName)          # pour un chemin
    end
    save(folderName*"/dico_path_"*nameExp*".jld2",dico_path)
end
