##################
### PM VERSION ###
##################

using Plots
using CSV
using DataFrames

#####################################
# Variables to adjust for given run #
#####################################

# Location of data files
workingDirectory = "/Users/justin_tackett/Documents/Research/Dr. Neilsen/Code/localGit/Post-Minkowski/julia_code/solver/src/Lyapunov Data"

pertDataName = "data_pn_p_1.csv" # Name of data file, perturbed
notPertDataName = "data_pn_np_1.csv" # Name of data file, not perturbed

f = 10^6 # Ratio of rescaling threshold to perturbation amount

################################
# Functions to be called later #
################################

function deltaYnorm(dataNP, dataP) # NP stands for not perturbed; input one timestep of data
    dataSize = size(dataNP)[1]
    sum = 0;

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


###########################################
# Main function to get Lyapunov exponents #
###########################################

function lyapunov(dataNP, dataP, f, ts)
    Rs = getRvalues(dataNP, dataP, f)
    print(string("These are the R values ", Rs, "\n")) # FIXME
    tFinal = dataNP[end, 1]
    dyfdy0 = deltaYnorm(dataNP[end,:],dataP[end,:])/deltaYnorm(dataNP[begin,:],dataP[begin,:])
    tSize = size(ts)[1]

    lambdaMax = (1/tFinal)*(log(dyfdy0)+ log(piProduct(Rs)))
    print(string("The lyapunov exponent is ",lambdaMax, "\n"))

    plotValues = [];
    tValues = ts

    for i in 1:tSize
        dydy0 = deltaYnorm(dataNP[i,:],dataP[i,:])/deltaYnorm(dataNP[begin,:],dataP[begin,:])
        plotValue = log(dydy0)
        push!(plotValues, plotValue)
    end

    plot(tValues, plotValues)
end

##############################################
# Calling main equations and reading in data #
##############################################

cd(workingDirectory)

dataNP_1 = readData(notPertDataName) # Getting data from files
dataP_1 = readData(pertDataName)

# Checking to make sure time steps and starting and end points for both datasets are the same
if dataNP_1[:,1] != dataP_1[:,1]
    print("Times for perturbed and unperturbed data sets are not equal; check data runs")
    return 1
else
    times = dataNP_1[:,1]
end

dataNP_2 = dataNP_1[:,1:end .!=1] # Data sets without the timesteps
dataP_2 = dataP_1[:,1:end .!=1]

lyapunov(dataNP_2, dataP_2, f, times)