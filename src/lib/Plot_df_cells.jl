module Plot_df_cells

using Revise
using FileIO, DataFrames, Statistics, LinearAlgebra
using PyPlot
using LaTeXStrings
# perso
using PIC_method

export plot_distri_angle, plot_distri_angle_zones, plot_distri_angle_BeforeAfter
export plot_distri_speed, plot_distri_speed_zones, plot_distri_velocity, plot_corr2D


""" Plot the distribution of velocities """
function plot_distri_velocity(df,maxUV,dUV,vmax,str_title,str_folder)
    # compute histogram
    xInt2,yInt2,ρ2 = PIC_method.PIC_2D(-maxUV,maxUV,-maxUV,maxUV,dUV,dUV,df[:,:u],-df[:,:v],true)
    # plot
    figure(4);clf()
    imshow(adjoint(ρ2),interpolation="bicubic",origin="lower",
           extent=(-maxUV,maxUV,-maxUV,maxUV),aspect="equal",vmin=0,vmax=vmax,cmap="jet")
    colorbar(extend="max",format="%.0e")
    scatter(0,0,s=20,color="black")
    xlabel(L"velocity $x \,(\mu m/h)$")
    ylabel(L"velocity $y \,(\mu m/h)$")
    title(str_title)
    savefig(str_folder*".pdf",bbox_inches="tight")
end


""" Plot the distribution of angle velocity """
function plot_distri_angle_zones(dict_df_zone,Δθ,zMax,array_yticks,str_title,str_folder,saving=true)
    # init
    name_zones = sort(collect(keys(dict_df_zone)))
    name_color = ["blue","gold","red","teal"]
    # plot
    figure(2);clf()
    for k in 1:length(dict_df_zone)
        ρ_zone,intTheta = distri_angle_zone(dict_df_zone[name_zones[k]],Δθ)
        polar(intTheta,2*π*ρ_zone,linewidth=3,label=name_zones[k],color=name_color[k])
    end
    intTheta = 0:.1:2π
    polar(intTheta,1*ones(length(intTheta)),linewidth=1,"--",color="black")
    # deco
    legend()
    axis((0,2π,0,zMax))
    yticks(array_yticks)
    xlabel(str_title*L", velocity angle $\theta$")
    if (saving)
        savefig(str_folder*".pdf",bbox_inches="tight")
    end
end


""" Plot the distribution of speed """
function plot_distri_speed_zones(dict_df_zone,dSpeed,maxSpeed,str_title,str_folder,saving=true)
    # init
    name_zones = sort(collect(keys(dict_df_zone)))
    name_color = ["blue","gold","red","teal"]
    # stat
    figure(3);clf()
    for k in 1:length(dict_df_zone)
        df_zone = dict_df_zone[name_zones[k]]
        println("-- stat "*name_zones[k]*" --")
        println("  mean = ",mean(df_zone[:,:speed]))
        println("  std  = ", std(df_zone[:,:speed]))
        intS,ρ = PIC_method.PIC_1D(0,maxSpeed,dSpeed,df_zone[:,:speed],true)
        plot(intS,ρ,label=name_zones[k],color=name_color[k],linewidth=3)
    end
    grid()
    legend()
    xlabel(L"speed distribution $\mu m/h$")
    ylabel(L"histogram")
    title(str_title)
    if (saving)
        savefig(str_folder*".pdf",bbox_inches="tight")
    end
end


""" distribution angle velocity """
function distri_angle_zone(df_zone,Δθ)
    ave_velocity_zone = [mean(cos.(-df_zone[:,:θ]));mean(-sin.(df_zone[:,:θ]));]
    ave_nem_velocity_zone = [mean(cos.(-2*df_zone[:,:θ]));mean(sin.(-2*df_zone[:,:θ]));]
    polarization_zone = norm(ave_velocity_zone)
    nematic_polarization_zone = norm(ave_nem_velocity_zone)
    intTheta,ρ_zone = PIC_method.PIC_1D(-π,π,Δθ,-df_zone[:,:θ],true)
    ρ_zone[1] = (ρ_zone[1]+ρ_zone[end])/2
    ρ_zone[end] = ρ_zone[1]
    return ρ_zone,intTheta
end

function plot_distri_angle(df,Δθ,zMax,array_yticks,str_title,str_name)
    """ Plot the distribution of angle velocity
    """
    # stat
    ave_velocity = [mean(cos.(-df[:,:θ]));mean(-sin.(df[:,:θ]));]
    ave_nem_velocity = [mean(cos.(-2*df[:,:θ]));mean(sin.(-2*df[:,:θ]));]
    polarization = norm(ave_velocity)
    nematic_polarization = norm(ave_nem_velocity)
    println("-- stat --")
    println("  polarization = ",polarization)
    println("  nem. polari. = ",nematic_polarization)
    # plot
    #-----
    figure(2)
    clf()
    intTheta,ρ = PIC_method.PIC_1D(-π,π,Δθ,-df[:,:θ],true)
    # periodicity
    ρ[1] = (ρ[1]+ρ[end])/2
    ρ[end] = ρ[1]
    polar(intTheta,2*π*ρ,linewidth=3)
    polar(intTheta,1*ones(length(intTheta)),linewidth=1,"--",color="black")
    #quiver(0,0,ave_velocity[1],-ave_velocity[2],scale_units="xy",scale=.3,color="red") # # θ,r,dx,dy
    quiver(0,0,ave_velocity[1],ave_velocity[2],scale=4,color="red") # # θ,r,dx,dy
    # deco
    axis((0,2π,0,zMax))
    yticks(array_yticks)
    xlabel(str_title*L", velocity angle $\theta$")
    savefig(str_name*".pdf",bbox_inches="tight")
end


function plot_distri_angle_BeforeAfter(df_before,df_after,Δθ,zMax,array_yticks,str_name)
    """ Plot the distribution of angle velocity
    """
    # stat before
    ave_velocity_before = [mean(cos.(-df_before[:,:θ]));mean(-sin.(df_before[:,:θ]));]
    ave_nem_velocity_before = [mean(cos.(-2*df_before[:,:θ]));mean(sin.(-2*df_before[:,:θ]));]
    polarization_before = norm(ave_velocity_before)
    nematic_polarization_before = norm(ave_nem_velocity_before)
    println("-- stat before--")
    println("  polarization = ",polarization_before)
    println("  nem. polari. = ",nematic_polarization_before)
    intTheta,ρ_before = PIC_method.PIC_1D(-π,π,Δθ,-df_before[:,:θ],true)
    ρ_before[1] = (ρ_before[1]+ρ_before[end])/2
    ρ_before[end] = ρ_before[1]
    # stat after
    ave_velocity_after = [mean(cos.(-df_after[:,:θ]));mean(-sin.(df_after[:,:θ]));]
    ave_nem_velocity_after = [mean(cos.(-2*df_after[:,:θ]));mean(sin.(-2*df_after[:,:θ]));]
    polarization_after = norm(ave_velocity_after)
    nematic_polarization_after = norm(ave_nem_velocity_after)
    println("-- stat after --")
    println("  polarization = ",polarization_after)
    println("  nem. polari. = ",nematic_polarization_after)
    intTheta,ρ_after = PIC_method.PIC_1D(-π,π,Δθ,-df_after[:,:θ],true)
    ρ_after[1] = (ρ_after[1]+ρ_after[end])/2
    ρ_after[end] = ρ_after[1]
    # plot
    #-----
    figure(2)
    clf()
    # periodicity
    polar(intTheta,2*π*ρ_before,linewidth=3,label="before",color="blue")
    polar(intTheta,2*π*ρ_after,linewidth=3,label="after",color="orange")
    polar(intTheta,1*ones(length(intTheta)),linewidth=1,"--",color="black")
    #quiver(0,0,ave_velocity[1],-ave_velocity[2],scale_units="xy",scale=.3,color="red") # # θ,r,dx,dy
    # deco
    legend()
    axis((0,2π,0,zMax))
    yticks(array_yticks)
    xlabel(str_name*L", velocity angle $\theta$")
    savefig(str_name*"_angle_velocity_BeforeAfter.pdf",bbox_inches="tight")
end


""" Plot the distribution of speed """
function plot_distri_speed(df,dSpeed,maxSpeed,str_title,str_name,df_after=nothing)
    # stat
    println("-- stat --")
    println("  mean = ",mean(df[:,:speed]))
    println("  std  = ",std(df[:,:speed]))
    intS,ρ = PIC_method.PIC_1D(0,maxSpeed,dSpeed,df[:,:speed],true)
    if (df_after!=nothing)
        println("-- stat (after) --")
        println("  mean (after) = ",mean(df_after[:,:speed]))
        println("  std (after)  = ",std(df_after[:,:speed]))
        intS,ρ_after = PIC_method.PIC_1D(0,maxSpeed,dSpeed,df_after[:,:speed],true)
    end
    # plot
    figure(3);clf()
    if (df_after==nothing)
        plot(intS,ρ,linewidth=3.0)
    else
        plot(intS,ρ,label="before")
        plot(intS,ρ_after,label="after")
    end
    xlabel(L"speed distribution $\mu m/h$")
    grid()
    xlabel(L"speed distribution $\mu m/h$")
    ylabel(L"histogram")
    title(str_title)
    savefig(str_name*".pdf",bbox_inches="tight")
end


function plot_corr2D(df,L,dx,str_title,str_folder)
    # init
    xInt2,yInt2,ρ2,w2 = PIC_method.PIC_2D(-L,L,-L,L,dx,dx,[df[:,:x_ij];df[:,:x_ji]],[df[:,:y_ij];df[:,:y_ji]],true,[df[:,:theCor];df[:,:theCor]])
    C = 2*L*2*L
    # relative position
    #------------------
    figure(5);clf()
    imshow(C*adjoint(ρ2),interpolation="bicubic",origin="lower",extent=(-L,L,-L,L),vmin=0,vmax=1.4,cmap="jet") # 3000, 1800, 2200
    colorbar(extend="max",format="%1.1f")#format="%1.5f")
    cs = contour(C*adjoint(ρ2),origin="lower",extent=(-L,L,-L,L),colors="black",levels=[.2,.5,1,1.2]) # 3000, 1800, 2200
    clabel(cs)
    xlabel("left -------- right")
    ylabel("back -------- front")
    title("Histogram neighbors, "*replace(str_title,"_"=>" "))
    s = range(0,stop=2π,length=100)
    #plot(20*cos.(s),20*sin.(s),linestyle="dashed",linewidth=3,color="white",alpha=.5)
    quiver(0,0,0,40, scale_units="xy",scale=1,color="white",alpha=.7,linewidth=5)
    axis("square")
    savefig(str_folder*"_relative_position.pdf",bbox_inches="tight")
    # correlation
    #------------
    figure(6);clf()
    imshow(adjoint(w2),interpolation="bicubic",origin="lower",extent=(-L,L,-L,L),vmin=-.6,vmax=.6,cmap="seismic")
    colorbar(extend="both",format="%1.1f")
    cs = contour(adjoint(w2),origin="lower",extent=(-L,L,-L,L),colors="black",levels=[-.2,-.1,.1,.2,.3,.4,.5])
    clabel(cs)
    #colorbar(extend="both",format="%.0e")
    title("Correlation, "*replace(str_title,"_"=>" "))
    xlabel("left -------- right")
    ylabel("back -------- front")
    plot(20*cos.(s),20*sin.(s),linestyle="dashed",linewidth=3,color="yellow",alpha=.8)
    quiver(0,0,0,40, scale_units="xy",scale=1,color="yellow",alpha=.8,linewidth=3)
    axis("square")
    savefig(str_folder*"_correlation_uv.pdf",bbox_inches="tight")
    # nematic correlation
    #--------------------
    cos_2θ = 2*df[:,:theCor].^2 .- 1
    xInt2,yInt2,ρ2,wNematic = PIC_method.PIC_2D(-L,L,-L,L,dx,dx,[df[:,:x_ij];df[:,:x_ji]],[df[:,:y_ij];df[:,:y_ji]],true,[cos_2θ;cos_2θ])
    figure(7);clf()
    imshow(adjoint(wNematic),interpolation="bicubic",origin="lower",extent=(-L,L,-L,L),vmin=-.25,vmax=.25,cmap="seismic")
    colorbar(extend="both",format="%1.1f")
    cs = contour(adjoint(wNematic),origin="lower",extent=(-L,L,-L,L),colors="black",levels=[-.2,-.15,-.1,-.05,0,.05,.1,.15,.2])
    clabel(cs)
    #colorbar(extend="both",format="%.0e")
    title("Nematic Correlation, "*replace(str_title,"_"=>" "))
    xlabel("left -------- right")
    ylabel("back -------- front")
    plot(20*cos.(s),20*sin.(s),linestyle="dashed",linewidth=3,color="yellow",alpha=.8)
    quiver(0,0,0,40, scale_units="xy",scale=1,color="yellow",alpha=.8,linewidth=3)
    axis("square")
    savefig(str_folder*"_nematic_correlation_uv.pdf",bbox_inches="tight")
end


end                             # module
