module Winding_number

using LinearAlgebra

export wind_number

function wind_number(X,P)
    """
       Return the winding number for the point X with respect to the (closed) path P.
    """
    # step 1: normalize the vector X-P_i
    N = length(P[:,1])
    u = zeros(N,2)
    for i=1:N
        u[i,:] = (P[i,:] - X)/norm(P[i,:] - X)
    end
    # step 2: angle between two successive vectors
    dθ = zeros(N-1)
    for i=1:(N-1)
        x = u[i,1]*u[i+1,1] + u[i,2]*u[i+1,2]           # dot product
        y = -u[i,2]*u[i+1,1] + u[i,1]*u[i+1,2]          # cross product
        dθ[i] = atan(y,x)
    end
    # step 3:
    return sum(dθ)/(2π)
end

# # testing
# P = [0 0;
#      1 0;
#      1 1;
#      0 1;
#      0 0]
# X = [.2; .5]
# wind_number(X,P)

end                             # module
