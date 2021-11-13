#=
Create/Plot individual velocities
=#

using FileIO, JSON, DataFrames
# perso
push!(LOAD_PATH,"lib")
using Plot_df_cells


# parameters
#-----------
list_nameExp = ["NPA_stich_3"]


#---------------------------------------------------#
#                     Mega loop                     #
#---------------------------------------------------#
for nameExp in list_nameExp
    println(nameExp)
    # loading
    folderDataTrackMate = "../data_tracking/step1_data_trackMate/"
    folderDataFilter = "../data_tracking/step2_data_filtered/"
    folderDataZone = "../data_tracking/step3_data_zone/"
    name_zones = ["zoneA","zoneB","zoneC","zoneD"]
    # find zones
    dictZone = Dict()
    open(folderDataZone*nameExp*"/coordinates.json", "r") do f
        dictZone = JSON.parse(f)
    end
    nbr_zones = length(dictZone)-1
    #------------------------#
    # load path
    dict_df_zone = Dict()
    str_folder = folderDataZone*nameExp*"/df_correlation/"
    for k in 1:nbr_zones
        df_zone = load(str_folder*nameExp*"_"*name_zones[k]*"_df_correlation.jld2")["df"]
        plot_corr2D(df_zone,42,7,
                    nameExp*" "*name_zones[k],
                    str_folder*nameExp*name_zones[k])
    end
end
