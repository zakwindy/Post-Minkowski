function lyaExponent(rb, ra, pb, pa, t)
    # ra and pa refer to the position and momentum of body a, similar thing for b
    n = size(ra)[1] # Number of snapshots

    rsum = 0
    psum = 0
    startTime = 10 # The time step corresponding to t = 1 N body units; chosen to be 1 to ensure virial equillibrium
    endTime = 0 # If loop stops prematurely, it is recorded which step it stops on 
    phaseSpaceDistances = [] # Measured in natural log space; slope of a fitted line for these is Lyapunov exponent

    for i in startTime:n
        rsum += (rb[i] - ra[i])^2 # Might need to replace running sums with just one term, need to check
        psum += (pb[i] - pa[i])^2 # Sum over all particles
        logdelta = .5 * log(rsum + psum)
        phaseSpace = exp(logdelta)

        push!(phaseSpaceDistances, phaseSpace)

        if phaseSpace >= .1
            endTime = i
            break
        end
    end
    
    
    



end