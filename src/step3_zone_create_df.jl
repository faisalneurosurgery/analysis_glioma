#=
Create/Plot individual velocities
=#

using FileIO, JSON, DataFrames
# perso
push!(LOAD_PATH,"lib")
using Create_df_cells

# parameters
#-----------
list_nameExp = ["NPA_stich_3"]
velocity_or_correlation = "correlation" # "velocity" or "correlation"
dist_boundary = 20.0              # in microns
rMax = 80.0                       # max distance pair for correlation
# loading
folderDataTrackMate = "../data_tracking/step1_data_trackMate/"
folderDataFilter = "../data_tracking/step2_data_filtered/"
folderDataZone = "../data_tracking/step3_data_zone/"
name_zones = ["zoneA","zoneB","zoneC","zoneD"]

#---------------------------------------------------#
#                     Mega loop                     #
#---------------------------------------------------#
for nameExp in list_nameExp
    println(nameExp)
    # find zones
    dictZone = Dict()
    open(folderDataZone*nameExp*"/coordinates.json", "r") do f
        dictZone = JSON.parse(f)
    end
    nbr_zones = length(dictZone)-1
    xMin_pixel = dictZone["bounding_box"]["xMin"]
    xMax_pixel = dictZone["bounding_box"]["xMax"]
    yMin_pixel = dictZone["bounding_box"]["yMin"]
    yMax_pixel = dictZone["bounding_box"]["yMax"]
    # size image
    dictP = Dict()
    open(folderDataTrackMate*nameExp*"/"*nameExp*"_info.json", "r") do f
        dictP = JSON.parse(f)
    end
    Lx = dictP["Width_microns"]
    Ly = dictP["Height_microns"]
    #------------------------#
    #    create dataframe    #
    #------------------------#
    # load path
    all_path = load(folderDataFilter*nameExp*"/dico_path_"*nameExp*".jld2")
    # create df
    #Tmin,Tmax = 0,1e10
    for k in 1:nbr_zones
        x_pixel = dictZone[name_zones[k]]["x"]
        y_pixel = dictZone[name_zones[k]]["y"]
        x = (x_pixel .- xMin_pixel)/(xMax_pixel-xMin_pixel)*Lx
        y = (y_pixel .- yMin_pixel)/(yMax_pixel-yMin_pixel)*Ly
        P = [x y]
        if (velocity_or_correlation == "velocity")
            df_zone = create_df_velocity_zone(all_path,P,dist_boundary,Lx-dist_boundary,dist_boundary,Ly-dist_boundary)
            folderName = folderDataZone*nameExp*"/df_velocity/"
            if (!isdir(folderName))
                mkpath(folderName)
            end
            save(folderName*nameExp*"_"*name_zones[k]*"_df_velocity.jld2","df",df_zone)
        else
            df_zone = create_df_correlation_zone(all_path,P,dist_boundary,Lx-dist_boundary,dist_boundary,Ly-dist_boundary,rMax)
            folderName = folderDataZone*nameExp*"/df_correlation/"
            if (!isdir(folderName))
                mkpath(folderName)
            end
            save(folderName*nameExp*"_"*name_zones[k]*"_df_correlation.jld2","df",df_zone)
        end
        println(size(df_zone))
    end
end
