module PIC_method

export PIC_2D, PIC_1D, Ngrid_2D, Ngrid_1D

#=
    Particle-In-Cell method to estimate density in 1D/2D:
       . nearest-grid: classical 'histogram' (order 1)
       . PIC: method of order 2
    Input:
       . xMin,xMax,yMin,yMax: bounding box
       . Δx,Δy: grid size
       . x,y: positions particles
       . vx,vy: velocity
    Output: xInt,yInt, density, u,v 

   Example:
   xMin,xMax,yMin,yMax = 0,5,0,8
   Δx,Δy = 1,1
   x = 5*rand(500000)
   y = 8*rand(500000)
   h = cos.(x.*y) + randn(500000)
   xInt,yInt,ρ = PIC_method.PIC_2D(xMin,xMax,yMin,yMax,Δx,Δy,x,y)
   xInt,yInt,ρ,u = PIC_method.PIC_2D(xMin,xMax,yMin,yMax,Δx,Δy,x,y,h)
   xInt,ρ,u= PIC_method.Ngrid_1D(xMin,xMax,Δx,x,h)

                                                             (nX,nY+1)   (nX+1,nY+1)
          +----------+----------+----------+----------+----------o----------o
          |          |          |          |          |          |          |
          |  (1,nY)  |          |          |          |          | (nX,nY)  |
          |          |          |          |          |          |          |
          +----------+----------+----------+----------+----------+----------+
          |          |          |          |          |          |          |
          |          |          |          |          |          |          |
          |          |          |          |          |          |          |
          +----------+----------+----------+----------+----------+----------+
          |          |          |          |          |          |          |
          |          |          |          |          |          |          |
          |          |          |          |          |          |          |
          +----------+----------+----------+----------+----------+----------+
          |          |          |          |          |          |          |
      Δy  |  (1,2)   |          |          |          |          |          |
          |          |          |          |          |          |          |
          +----------+----------+----------+----------+----------+----------+
          |          |          |          |          |          |          |
          |  (1,1)   |   (2,1)  |          |          |          |  (nX,1)  |
          |          |          |          |          |          |          |
          o----------+----------+----------+----------+----------+----------+
     (xMin,yMin)          Δx

=#



function PIC_2D(xMin,xMax,yMin,yMax,Δx,Δy,x,y,normalize=false,vx=nothing,vy=nothing)
    """
    Compute the average density and velocity of particles using PIC method of order 2
      . x,y   : lists of x and y components of the vector positions of the particles
      . xMin,...,yMax: edges of the grid
      . Δx,Δy : meshgrid
      . vx,vy : x and y components of the velocity variable (or other things)
    """
    # init
    N = length(x)
    nX = floor(Int64,(xMax-xMin)/Δx)
    nY = floor(Int64,(yMax-yMin)/Δy)
    ρ  = zeros(nX+1,nY+1)
    if (vx != nothing)
        fu = zeros(nX+1,nY+1)
        u  = zeros(nX+1,nY+1)
    end
    if (vy != nothing)
        fv = zeros(nX+1,nY+1)
        v  = zeros(nX+1,nY+1)
    end
    #------#
    # Grid #  'at the boundary'
    #------#
    for k in 1:N
        i = floor(Int64,(x[k]-xMin)/Δx) + 1
        j = floor(Int64,(y[k]-yMin)/Δy) + 1
        if (1<=i) & (i<=nX) & (1<=j) & (j<=nY)
            #---     (i-1)Δx < x < iΔx  ,   (j-1)Δy < y < jΔy
            x_p = x[k] - (xMin + (i-1)*Δx)
            y_p = y[k] - (yMin + (j-1)*Δy)
            #- Add mass ΔxΔy on the grid
            ρ[i,j]      += (Δx - x_p)*(Δy - y_p)
            ρ[i+1,j+1]  += x_p*y_p
            ρ[i+1,j]    += x_p*(Δy-y_p)
            ρ[i,j+1]    += (Δx-x_p)*y_p
            if (vx != nothing)
                fu[i,j]     += (Δx - x_p)*(Δy - y_p)*vx[k]
                fu[i+1,j+1] += x_p*y_p*vx[k]
                fu[i+1,j]   += x_p*(Δy-y_p)*vx[k]
                fu[i,j+1]   += (Δx-x_p)*y_p*vx[k]
            end
            if (vy != nothing)
                fv[i,j]     += (Δx - x_p)*(Δy - y_p)*vy[k]
                fv[i+1,j+1] += x_p*y_p*vy[k]
                fv[i+1,j]   += x_p*(Δy-y_p)*vy[k]
                fv[i,j+1]   += (Δx-x_p)*y_p*vy[k]
            end
        end
    end
    # deduce u and v
    if (vx != nothing)
        for i in 1:(nX+1), j in 1:(nY+1)
            if (ρ[i,j]>0)
                u[i,j] = fu[i,j]/ρ[i,j]
            end
        end
    end
    if (vy != nothing)
        for i in 1:(nX+1), j in 1:(nY+1)
            if (ρ[i,j]>0)
                v[i,j] = fv[i,j]/ρ[i,j]
            end
        end
    end
    # normalize ρ
    if (normalize)
        ρ ./= sum(ρ)*Δx*Δy
    else
        ρ ./= Δx*Δy
    end
    # boundary
    for i in 1:(nX+1)
        ρ[i,1]  *= 2
        ρ[i,nY+1] *= 2
    end
    for j in 1:(nY+1)
        ρ[1,j]  *= 2
        ρ[nX+1,j] *= 2
    end
    # finally
    xInt = xMin .+ (0:nX)*Δx
    yInt = yMin .+ (0:nY)*Δy
    if (vx == nothing)
        return xInt,yInt,ρ
    elseif (vy == nothing)
        return xInt,yInt,ρ,u
    else
        return xInt,yInt,ρ,u,v
    end
end

function PIC_1D(xMin,xMax,Δx,x,normalize=false,v=nothing)
    """
     Same as PIC_2D with only one variable.
    """
    # init
    N = length(x)
    nX = floor(Int64,(xMax-xMin)/Δx)
    ρ  = zeros(nX+1)
    if (v != nothing)
        f = zeros(nX+1)
        u = zeros(nX+1)
    end
    #------#
    # Grid #  'at the boundary'
    #------#
    for k in 1:N
        i = floor(Int64,(x[k]-xMin)/Δx) + 1
        if (1<=i) & (i<=nX)
            #---     (i-1)Δx < x < iΔx
            x_p = x[k] - (xMin + (i-1)*Δx)
            #- Add mass ΔxΔy on the grid
            ρ[i]    += (Δx - x_p)
            ρ[i+1]  += x_p
            if (v != nothing)
                f[i]   += (Δx - x_p)*v[k]
                f[i+1] += x_p*v[k]
            end
        end
    end
    # deduce u
    if (v != nothing)
        for i in 1:(nX+1)
            if (ρ[i]>0)
                u[i] = f[i]/ρ[i]
            end
        end
    end
    # normalize ρ
    if (normalize)
        ρ ./= sum(ρ)*Δx
    else
        ρ ./= Δx
    end
    # boundary
    ρ[1]    *= 2
    ρ[nX+1] *= 2
    # finally
    xInt = xMin .+ (0:nX)*Δx
    if (v == nothing)
        return xInt,ρ
    else
        return xInt,ρ,u
    end
end


function Ngrid_2D(xMin,xMax,yMin,yMax,Δx,Δy,x,y,normalize=false,vx=nothing,vy=nothing)
    """
     Compute the density and velocity of particle using PIC method of order 1 (nearest grid point)
    """
    # init
    N = length(x)
    nX = floor(Int64,(xMax-xMin)/Δx)
    nY = floor(Int64,(yMax-yMin)/Δy)
    ρ  = zeros(nX,nY)
    if (vx != nothing)
        fu = zeros(nX,nY)
        u  = zeros(nX,nY)
    end
    if (vy != nothing)
        fv = zeros(nX,nY)
        v  = zeros(nX,nY)
    end
    #------#
    # Grid #  'at the boundary'
    #------#
    for k in 1:N
        i = floor(Int64,(x[k]-xMin)/Δx) + 1
        j = floor(Int64,(y[k]-yMin)/Δy) + 1
        if (1<=i) & (i<=nX) & (1<=j) & (j<=nY)
            #---     (i-1)Δx < x < iΔx  ,   (j-1)Δy < y < jΔy
            x_p = x[k] - (xMin + (i-1)*Δx)
            y_p = y[k] - (yMin + (j-1)*Δy)
            #- Add mass ΔxΔy on the grid
            ρ[i,j]  += 1
            if (vx != nothing)
                fu[i,j] += vx[k]
            end
            if (vy != nothing)
                fv[i,j] += vy[k]
            end
        end
    end
    # deduce u and v
    if (vx != nothing)
        for i in 1:(nX), j in 1:(nY)
            if (ρ[i,j]>0)
                u[i,j] = fu[i,j]/ρ[i,j]
            end
        end
    end
    if (vy != nothing)
        for i in 1:(nX), j in 1:(nY)
            if (ρ[i,j]>0)
                v[i,j] = fv[i,j]/ρ[i,j]
            end
        end
    end
    # normalize ρ
    if (normalize)
        ρ ./= sum(ρ)*Δx*Δy
    else
        ρ ./= Δx*Δy
    end
    # finally
    xInt = xMin .+ (1:nX)*Δx .- Δx/2
    yInt = yMin .+ (1:nY)*Δy .- Δy/2
    if (vx == nothing)
        return xInt,yInt,ρ
    elseif (vy == nothing)
        return xInt,yInt,ρ,u
    else
        return xInt,yInt,ρ,u,v
    end
end

function Ngrid_1D(xMin,xMax,Δx,x,normalize=false,v=nothing)
    """
     Same as PIC_2D with only one variable.
    """
    # init
    N = length(x)
    nX = floor(Int64,(xMax-xMin)/Δx)
    ρ  = zeros(nX)
    if (v != nothing)
        f = zeros(nX)
        u = zeros(nX)
    end
    #------#
    # Grid #  'at the boundary'
    #------#
    for k in 1:N
        i = floor(Int64,(x[k]-xMin)/Δx) + 1
        if (1<=i) & (i<=nX)
            #---     (i-1)Δx < x < iΔx
            x_p = x[k] - (xMin + (i-1)*Δx)
            #- Add mass ΔxΔy on the grid
            ρ[i] += 1
            if (v != nothing)
                f[i] += v[k]
            end
        end
    end
    # deduce u
    if (v != nothing)
        for i in 1:nX
            if (ρ[i]>0)
                u[i] = f[i]/ρ[i]
            end
        end
    end
    # normalize ρ
    if (normalize)
        ρ ./= sum(ρ)*Δx
    else
        ρ ./= Δx
    end
    # finally
    xInt = xMin .+ (1:nX)*Δx .- Δx/2
    if (v == nothing)
        return xInt,ρ
    else
        return xInt,ρ,u
    end
end

end                             # end module
