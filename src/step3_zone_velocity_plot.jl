#=
Create/Plot individual velocities
=#

using FileIO, JSON, DataFrames, Statistics, CSV
# perso
push!(LOAD_PATH,"lib")
using Plot_df_cells

# loading
folderDataTrackMate = "../data_tracking/step1_data_trackMate/"
folderDataFilter = "../data_tracking/step2_data_filtered/"
folderDataZone = "../data_tracking/step3_data_zone/"
name_zones = ["zoneA","zoneB","zoneC","zoneD"]

# parameters
#-----------
list_nameExp = ["NPA_stich_3"]

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
    #------------------------#
    # load path
    dict_df_zone = Dict()
    str_folder = folderDataZone*nameExp*"/df_velocity/"
    for k in 1:nbr_zones
        df_zone = load(str_folder*nameExp*"_"*name_zones[k]*"_df_velocity.jld2")["df"]
        plot_distri_velocity(df_zone,5,.5,.04,
                             nameExp*" "*name_zones[k],
                             str_folder*nameExp*name_zones[k]*"_velocity")
        dict_df_zone[name_zones[k]] = df_zone
        # print stat
        pola = sqrt(mean(cos.(df_zone[:,:θ])).^2 + mean(sin.(df_zone[:,:θ])).^2)
        pola_nem = sqrt(mean(cos.(2*df_zone[:,:θ])).^2 + mean(sin.(2*df_zone[:,:θ])).^2)
        println(" stat "*name_zones[k])
        println("   polarization  : ",pola)
        println("   pola. nematic : ",pola_nem)
        # save
        CSV.write(str_folder*nameExp*"_"*name_zones[k]*".csv",select(df_zone,:θ))
    end
    # plot
    plot_distri_angle_zones(dict_df_zone,2π/20,3,[1,2,3],nameExp,str_folder*nameExp*"_angle_velocity",true)
    plot_distri_speed_zones(dict_df_zone,.5,30,nameExp,str_folder*nameExp*"_speed",true)
end
