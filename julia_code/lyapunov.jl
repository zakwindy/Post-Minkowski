function lyaExponent(rb, ra, pb, pa, numbodies)
    # ra and pa refer to the position and momentum of body a, similar thing for b
    n = size(ra)[1] # Number of snapshots
    nb = numbodies # Number of bodies in simulation
    pbd = pertBod # Body index being perturbed

    rsum = 0
    psum = 0
   
    endTime = 0 # If loop stops prematurely, it is recorded which step it stops on 
    phaseSpaceDistances = [] # Measured in natural log space; slope of a fitted line for these is Lyapunov exponent

    for i in 1:n # Check this later
        
        for j in 1:nb
            rsum += (rb[i][j] - ra[i][j])^2 # Might need to replace running sums with just one term, need to check
            psum += (pb[i][j] - pa[i][j])^2 # Sum over all particles
        end

        logdelta = .5 * log(rsum + psum)
        phaseSpace = exp(logdelta)

        push!(phaseSpaceDistances, phaseSpace)

        if phaseSpace >= .1
            endTime = i
            break
        end
        
    end
        
    return phaseSpaceDistances

end