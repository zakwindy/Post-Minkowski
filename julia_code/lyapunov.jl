using Plots
using CSV
using DataFrames

##############################
#Functions to be called later#
##############################

function deltaYnorm(dataNP, dataP) # NP stands for not perturbed; input one timestep of data
    dataSize = size(dataNP)[1]
    sum = 0;

    #realDataNP = dataNP[:, 1:end .!=1] # Getting rid of time steps from phase space
    #realDataP = dataP[:, 1:end .!=1]   # Thats not what actually happened, keeping this in here just in case

    for i in 1:dataSize
        sum += (dataNP[i] - dataP[i])^2
    end
    return sqrt(sum)
end

function piProduct(numberList)
    listSize = size(numberList)[1]
    product = 1;

    for i in 1:listSize
        product *= numberList[i]
        
    end
    return product
end

# how to deal with different timesteps, sol to csv, preserving dense solutions
# Double check the final time and time steps

function getRvalues(dataNP, dataP, f) # Put in full data sets; f should be something like 10^6
    dataSize = size(dataNP)[1]
    Rs = [1.0]
    dy0 = deltaYnorm(dataNP[1,:], dataP[1,:]) # Get rows in call to data sets

    for i in 1:dataSize
        dy = deltaYnorm(dataNP[i,:], dataP[i,:])
        potentialR = dy/dy0
        
        if dy >= piProduct(Rs)*f*dy0
            potentialR = dy/dy0
            push!(Rs, potentialR)
        end
    end

    return Rs
end

function readData(dataName::String)
    df = DataFrame(CSV.File(dataName))
    data = Matrix(df)
    return data
end


#########################################
#Main function to get Lyapunov exponents#
#########################################

function lyapunov(dataNP, dataP, f)
    Rs = getRvalues(dataNP, dataP, f)
    tFinal = dataNP[end, 1]
    dyfdy0 = deltaYnorm(dataNP[end,:],dataP[end,:])/deltaYnorm(dataNP[begin,:],dataP[begin,:])
    tSize = size(dataNP)[1]

    lambdaMax = (1/tFinal)*(log(dyfdy0)+ log(piProduct(Rs)))
    print(string("The lyapunov exponent is ",lambdaMax))

    plotValues = []
    tValues = dataNP[:, 1]

    #for i in 1:tSize
    #    dydy0 = deltaYnorm(dataNP[i,:],dataP[i,:])/deltaYnorm(dataNP[begin,:],dataP[begin,:])
    #    plotValue = log(dydy0)
    #    push!(plotValues, plotValue)
    #end

    #plot(tValues, plotValues)
end

workingDirectory = "/Users/justin_tackett/Documents/Research/Dr. Neilsen/Code/localGit/Post-Minkowski/julia_code/solver/src"

cd(workingDirectory)

dataNP_1 = readData("data_pn_np_1.csv") # Getting data from files
dataP_1 = readData("data_pn_p_1.csv")

dataNP_temp = Matrix(dataNP_1)
dataP_temp = Matrix(dataP_1)

dataNP_2 = dataNP_temp[:,1:end .!=1]
dataP_2 = dataP_temp[:,1:end .!=1]

lyapunov(dataNP_2, dataP_2, 10^6)

# What to do next:
#     go through and figure out why its a 27x0 array access, check line 142 to see whats going on with [1]
#     figure out what the heck is going on with sizeof(), going to 216 instead of 27