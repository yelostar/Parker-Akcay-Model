#!/usr/bin/env julia

using JLD2
using StatsBase
using FileIO
using Graphs
#using Plots


mutable struct NetworkParameters

    #measurement data
    meanCoopFreq::Float64
    meanProbNeighborCoop::Float64
    meanProbNeighborDef::Float64
    meanProbNeighborAcc::Float64
    meanProbRandomCoop::Float64
    meanProbRandomDef::Float64
    meanProbRandomAcc::Float64
    meanDegree::Float64
    meanDistInclusion::Float64
    meanFitness::Float64
    meanShortestPaths::Float64
    meanConnComponents::Float64
    meanConnCompSize::Float64
    meanLargestConnComp::Float64
    meanDistConnection::Float64

    #Node characteristics
    popPNC::Array{Float64, 1}
    popPND::Array{Float64, 1}
    popPRC::Array{Float64, 1}
    popPRD::Array{Float64, 1}
    popPNA::Array{Float64, 1}
    popPRA::Array{Float64, 1}
    popStrategies::Array{Int64, 1}
    popPayoff::Array{Float64, 1}
    popFitness::Array{Float64, 1}

    #structural variables
    numGens::Int64
    popSize::Int64
    edgeMatrix::Array{Int64, 2}

    #social variables
    cost::Float64
    benefit::Float64
    synergism::Float64
    linkCost::Float64
    muS::Float64
    muP::Float64    
    delta::Float64

    #simulation determining
    pnd::Bool
    prd::Bool
    allowReject::Bool
    distInherit::Bool

    #evolving links variables
    sigmapn::Float64
    sigmapr::Float64

    #distance inherit variables
    distFactor::Float64 

    function NetworkParameters(b::Float64, c::Float64, d::Float64, cL::Float64, gen::Int, distInherit::Bool, distFactor::Float64, pn::Float64, pnd::Bool, pr::Float64, prd::Bool, allowReject::Bool, pa::Float64, muP::Float64, delta::Float64, sigmapn::Float64, sigmapr::Float64)

        popSize = 100
        popPNC = zeros(Float64, popSize)
        popPNC[:] .= pn
        popPND = zeros(Float64, popSize)
        popPND[:] .= pn
        popPRC = zeros(Float64, popSize)
        popPRC[:] .= pr
        popPRD = zeros(Float64, popSize)
        popPRD[:] .= pr

        #probability to accept connections: PNA (paternal contacts) and PRA (random contacts)
        popPNA = zeros(Float64, popSize)
        popPRA = zeros(Float64, popSize)
        if(allowReject) #if allowReject is false, PNA and PRA are set to 1.0 and are not allowed to separately evolve
            popPNA[:] .= pa
            popPRA[:] .= pa
        else
            popPNA[:] .= 1.0
            popPRA[:] .= 1.0
        end

        popStrategies = zeros(Int64, popSize)
        popStrategies[2:2:popSize] .= 1
        popFitness = zeros(Float64, popSize)
        popFitness[:] .= 1.0
        
        if(distInherit==true) #alternate set up of edgeMatrix
            edgeMatrix = zeros(popSize, popSize)
            for(i) in 1:popSize
                linkWeights = distFactor.^locSelect(popSize, i) #creates weighted list of locations based on distFactor
                links = zeros(Float64, Int(popSize/2)) 
                sample!((1:popSize), Weights(linkWeights), links, replace=false) #picks 50 links for each individual with probabilities weighted by locations
                for(j) in links
                    edgeMatrix[i, Int(j)] = 1    
                end
            end       
        else
            edgeMatrix=rand([0,1],(popSize,popSize))
        end
        edgeMatrix = edgeMatrix .* transpose(edgeMatrix) #ensures that connections are reciprocated
        for(i) in 1:popSize
            edgeMatrix[i, i] = 0 #cannot connect with themself
        end

        cost = c
        synergism = d
        benefit = b
        linkCost = cL
        muS = muP #changing strategies
        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, popPNC, popPND, popPRC, popPRD, popPNA, popPRA, popStrategies, zeros(Float64, popSize), popFitness, gen, popSize, edgeMatrix, cost, benefit, synergism, linkCost, muS, muP, delta, pnd, prd, allowReject, distInherit, sigmapn, sigmapr, distFactor)
    end
end

function locSelect(popSize::Int64, i::Int64)
    locs = 1:popSize
    distances = min.(abs.(locs.-locs[i]), abs.(locs.-locs[i].-popSize),abs.(locs.-locs[i].+popSize))
    return distances
end   

#measurement functions
function coopRatio(network::NetworkParameters)
    coopCount = 0
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            coopCount += 1
        end
    end
    coopCount /= network.popSize
    network.meanCoopFreq += coopCount
end

function probInherit(network::NetworkParameters)
    network.meanProbRandomCoop += sum(network.popPRC)/network.popSize
    network.meanProbRandomDef += sum(network.popPRD)/network.popSize
    network.meanProbRandomAcc += sum(network.popPRA)/network.popSize
    network.meanProbNeighborCoop += sum(network.popPNC)/network.popSize
    network.meanProbNeighborDef += sum(network.popPND)/network.popSize
    network.meanProbNeighborAcc += sum(network.popPNA)/network.popSize
    
end

function degrees(network::NetworkParameters)
    degTotal = 0
    fitnessTotal = 0
    connDistTotal = 0
    numConnIndivids = 0
    for(i) in 1:network.popSize
        fitnessTotal+=network.popFitness[i]
        degCounter = 0.0
        connDistCounter = 0.0
        connsVec = locSelect(network.popSize, i)
        for(ii) in 1:network.popSize
            if(network.edgeMatrix[i, ii] != 0)
                degCounter += 1.0
                connDistCounter += connsVec[ii] #adds distance between individuals' locations
                #originally connDistCounter += min(abs(network.popLocations[i]-network.popLocations[ii]), abs.(network.popLocations[i]-network.popLocations[ii]-network.popSize),abs.(network.popLocations[i]-network.popLocations[ii]+network.popSize)) 
            end
        end
        degTotal += degCounter
        if(degCounter != 0) #prevents divide by 0
            connDistCounter /= degCounter #divides sum of distances by total number of connections
            numConnIndivids += 1
        end
        connDistTotal += connDistCounter
    end
    degTotal /= network.popSize
    connDistTotal /= numConnIndivids
    fitnessTotal /= network.popSize
    network.meanDegree += degTotal
    network.meanFitness += fitnessTotal
    network.meanDistConnection += connDistTotal
end

#evolution functions

function deathbirth(network::NetworkParameters, findMom::Function, neighborRange::Int64, dbProb::Float64=1.0)
    childID = deathfirst(network)
    parentID = findMom(network, childID, neighborRange)
    birth(network, childID, parentID)
end

function birthdeath(network::NetworkParameters, findMom::Function=anyMom, neighborRange::Int64=1, dbProb::Float64=0.0)
    fitWeights = StatsBase.weights(network.popFitness)
    parentID = sample(1:network.popSize, fitWeights)
    childID = deathsecond(network, parentID)
    birth(network, childID, parentID)
end

function mixeddb(network::NetworkParameters, findMom::Function, neighborRange::Int64, dbProb::Float64)
    if(rand() < dbProb)
        deathbirth(network, findMom, neighborRange)
    else
        birthdeath(network)
    end
end

function deathfirst(network::NetworkParameters)
    deadID = sample(1:network.popSize)
    network.edgeMatrix[deadID, :] .= 0
    network.edgeMatrix[:, deadID] .= 0
    network.popPayoff[deadID] = 0
    network.popFitness[deadID] = 0
    deadID
end

function deathsecond(network::NetworkParameters, parentID::Int64)
    deadID  = (parentID + sample([1, -1]) + 99)%100 + 1 #picks a individual to die randomly from parent's neighbors
    network.edgeMatrix[deadID, :] .= 0
    network.edgeMatrix[:, deadID] .= 0
    network.popPayoff[deadID] = 0
    network.popFitness[deadID] = 0
    deadID
end

function anyMom(network::NetworkParameters, kID::Int64, neighborRange::Int64=50)
    fitWeights = StatsBase.weights(network.popFitness)
    momIndex = sample(1:network.popSize, fitWeights)
    momIndex
end

function neighborMom(network::NetworkParameters, kID::Int64, neighborRange::Int64=1) #100 is hard-coded for popSize
    indexes = zeros(Int64, 2*neighborRange)
    fitnesses = zeros(2*neighborRange)
    for(i) in 1:neighborRange
        indexes[i] = (kID-i+99)%100+1 #index minus 1 through neighborRange, loops around 100
        indexes[i+neighborRange] = (kID+i+99)%100+1 #index plus 1 through neighborRange, loops around 100
        fitnesses[i] = network.popFitness[indexes[i]] #same for fitnesses
        if(indexes[i]==indexes[i+neighborRange]) #necessary for when neighborRange = popSize/2
            break
        end
        fitnesses[i+neighborRange] = network.popFitness[indexes[i+neighborRange]]
    end
    fitWeights = StatsBase.weights(fitnesses)
    momIndex = sample(indexes, fitWeights) 
    momIndex
end

function birth(network::NetworkParameters, child::Int64, parent::Int64)
    
    network.popFitness[child] = 1
    network.popStrategies[child] = network.popStrategies[parent]
    if(rand()<network.muS)
        network.popStrategies[child] -= 1
        network.popStrategies[child] *= -1
    end
    
    network.popPNC[child] = network.popPNC[parent]
    if(rand()<network.muP)
        network.popPNC[child] += randn()*network.sigmapn
        network.popPNC[child] = clamp(network.popPNC[child], 0, 1)
    end
    
    if(network.pnd) #allows separate evolution of PND if pnd==true
        network.popPND[child] = network.popPND[parent]
        if(rand()<network.muP)
            network.popPND[child] += randn()*network.sigmapn
            network.popPND[child] = clamp(network.popPND[child], 0, 1)
        end
    else
        network.popPND[child] = network.popPNC[child]
    end

    network.popPRC[child] = network.popPRC[parent]
    if(rand()<network.muP)
        network.popPRC[child] += randn()*network.sigmapr
        network.popPRC[child] = clamp(network.popPRC[child], 0, 1)
    end

    if(network.prd) #allows separate evolution of PRD if prd==true
        network.popPRD[child] = network.popPRD[parent]
        if(rand()<network.muP)
            network.popPRD[child] += randn()*network.sigmapr
            network.popPRD[child] = clamp(network.popPRD[child], 0, 1)
        end
    else
        network.popPRD[child] = network.popPRC[child]
    end

    if(network.allowReject) #allows evolution of PNA and PRA if true
        network.popPNA[child] = network.popPNA[parent]
        if(rand()<network.muP)
            network.popPNA[child] += randn()*network.sigmapn
            network.popPNA[child] = clamp(network.popPNA[child], 0, 1)
        end
        network.popPRA[child] = network.popPRA[parent]
        if(rand()<network.muP)
            network.popPRA[child] += randn()*network.sigmapr
            network.popPRA[child] = clamp(network.popPRA[child], 0, 1)
        end
    else
        network.popPNA[child] = 1.0
        network.popPRA[child] = 1.0
    end
         
    if(network.distInherit)
        locInherit(network, child, parent)
    else
        nonDistInherit(network, child, parent)
    end

    network.edgeMatrix[parent, child] = 1
    network.edgeMatrix[child, parent] = 1

end

#inheritance functions
function locInherit(network::NetworkParameters, child::Int64, parent::Int64)
    
    dists = network.distFactor.^locSelect(network.popSize, child)
    #originally dists = network.distFactor.^(min.(abs.(network.popLocations.-network.popLocations[child]), abs.(network.popLocations.-network.popLocations[child].-network.popSize),abs.(network.popLocations.-network.popLocations[child].+network.popSize))) #creates weighted list of locations based on distFactor    for(i) in 1:network.popSize

    for(i) in 1:network.popSize
        if(i != child && network.edgeMatrix[i, child] == 0)
            if(network.edgeMatrix[i, parent] != 0)
                if(network.popStrategies[i] == 1) #PNC MODE
                    if(rand() < network.popPNC[child] * network.popPNA[i] * dists[i]) #prob to inherit multiplied by value for distance stored in dists
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPND[child] * network.popPNA[i] * dists[i]) #PND MODE
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
            else
                if(network.popStrategies[i] == 1) #PRC MODE
                    if(rand() < network.popPRC[child] * network.popPRA[i] * dists[i])
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPRD[child] * network.popPRA[i] * dists[i]) #PRD MODE
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
 
            end
        end
    end
end

function nonDistInherit(network::NetworkParameters, child::Int64, parent::Int64)
    for(i) in 1:network.popSize
        if(i != child && network.edgeMatrix[i, child] == 0)
            if(network.edgeMatrix[i, parent] != 0)
                if(network.popStrategies[i] == 1) #PNC MODE
                    if(rand() < network.popPNC[child] * network.popPNA[i]) #PNC of child * PNA of individual newborn is connecting with 
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPND[child] * network.popPNA[i]) #PND of child * PNA of individual newborn is connecting with 
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
            else
                if(network.popStrategies[i] == 1) #PRC of child * PRA of individual newborn is connecting with 
                    if(rand() < network.popPRC[child] * network.popPRA[i])
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                else
                    if(rand() < network.popPRD[child] * network.popPRA[i]) #PRD of child * PRA of individual newborn is connecting with 
                        network.edgeMatrix[i, child] = 1
                        network.edgeMatrix[child, i] = 1
                    end
                end
                #if(rand() < network.popPR[child]) #original code
                #network.edgeMatrix[i, child] = 1
                #network.edgeMatrix[child, i] = 1
                #end
            end
        end
    end
end

function resolveFitnesses(network::NetworkParameters)
    for(i) in 1:network.popSize
        network.popFitness[i] = (1.0 + network.delta) ^ network.popPayoff[i]
    end
    network.popPayoff[:] .= 0.0
end

function cooperate(network::NetworkParameters)
    degs = getDegree(network)
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            B = network.benefit/degs[i]
            D = network.synergism/degs[i]

            for(ii) in 1:network.popSize
                if(network.edgeMatrix[i, ii] != 0)
                    network.popPayoff[ii] += B
                    if(network.popStrategies[ii] == 1)
                        network.popPayoff[ii] += D/degs[ii]
                    end
                end
            end

            network.popPayoff[i] -= network.cost
        end
        network.popPayoff[i] -= network.linkCost * degs[i]
    end
    resolveFitnesses(network)
end 

function getDegree(network::NetworkParameters) #made less efficient by 2 in edgeMatrix; can use sum for 1/0 vals
    degGetter = zeros(100)
    for(i) in 1:100
        for(ii) in 1:100
            if(network.edgeMatrix[i,ii] != 0)
                degGetter[i] += 1
            end
        end
        #degGetter[i] = sum(edgeMatrix[i,:])
    end

    degGetter
end

function graphCalc(network::NetworkParameters) #computes all calculations involving Graphs.jl, including shortest paths & connected components
    edgeGraph = Graph(network.edgeMatrix)

    allPath = 0.0
    disconn = 0 #keeps track of disconnected components
    for(i) in 1:network.popSize
        dsp = dijkstra_shortest_paths(edgeGraph, i)
        paths = dsp.dists
        undefs = zeros(Int, 0) #keeps track of infinite and zero distances between nodes
        for(j) in 1:(size(paths)[1])
            if(paths[j] > (network.popSize+1) || paths[j] == 0)
                push!(undefs, j)
            end
        end
        if(size(undefs)==size(paths)) #if all shortest paths are zero or infinite, disconn increments
            disconn += 1
        else
            for(k) in 1:size(undefs)[1]
                deleteat!(paths, undefs[k]) #removes all infinite and zero distances from paths
                undefs .-= 1 #avoids index out of bounds
            end
            allPath += mean(paths)
        end
    end
    network.meanShortestPaths += (allPath / (network.popSize-disconn))

    cc = connected_components(edgeGraph) #connected_components(edgeGraph) returns a vector of vectors of each connected component
    network.meanConnComponents += size(cc)[1] 
    network.meanConnCompSize += mean(getfield.(size.(cc), 1))
    network.meanLargestConnComp += maximum(getfield.(size.(cc), 1))

end

function runSimsReturn(;B::Float64=2.0, C::Float64=0.5, D::Float64=0.0, CL::Float64=0.0, gen::Int=500, dbOrder::Function=mixeddb, dbProb::Float64=1.0, findMom::Function=anyMom, neighborRange::Int64=1, distInherit::Bool=false, distFactor::Float64=0.975, pn::Float64=0.5, pnd::Bool=false, pr::Float64=0.01, prd::Bool=false, allowReject::Bool=false, pa::Float64=0.75, muP::Float64=0.001, delta::Float64=0.1, sigmapn::Float64=0.05, sigmapr::Float64=0.01, reps::Int64=50)
    dataArray = zeros(15) 
    repSims = reps
    for(x) in 1:repSims

        #initializes globalstuff structure with generic constructor
        network = NetworkParameters(B, C, D, CL, gen, distInherit, distFactor, pn, pnd, pr, prd, allowReject, pa, muP, delta, sigmapn, sigmapr) #if pnd/prd = true, then defector & cooperator probabilities will evolve separately for pn/pr

        #checks efficiency of simulation while running it
        for(g) in 1:(network.numGens * network.popSize)

            if(g > (20))
                cooperate(network)
            end

            dbOrder(network, findMom, neighborRange, dbProb)

            if(g > (network.numGens * network.popSize / 5) && (g % network.popSize) == 0) #metrics
                coopRatio(network)
                probInherit(network)
                degrees(network)
                graphCalc(network)
            end

        end

        #divides meanCooperationRatio by last 400 generations to get a true mean, then outputs
        network.meanProbNeighborCoop /= (network.numGens*0.8)
        network.meanProbNeighborDef /= (network.numGens*0.8)
        network.meanProbNeighborAcc /= (network.numGens*0.8)
        network.meanProbRandomCoop /= (network.numGens*0.8)
        network.meanProbRandomDef /= (network.numGens*0.8)
        network.meanProbRandomAcc /= (network.numGens*0.8)
        network.meanDegree /= (network.numGens*0.8)
        #network.meanAssortment /= (network.numGens*0.8)
        network.meanCoopFreq /= (network.numGens*0.8)
        #network.meanCoopDefDistance /= (network.popSize*network.numGens*0.8)
        network.meanDistInclusion /= (network.popSize*network.numGens*0.8)
        network.meanFitness /= (network.numGens*0.8)
        network.meanShortestPaths /= (network.numGens*0.8)
        network.meanConnComponents /= (network.numGens*0.8)
        network.meanConnCompSize /= (network.numGens*0.8)
        network.meanLargestConnComp /= (network.numGens*0.8)
        network.meanDistConnection /= (network.numGens*0.8)

        dataArray[1] += network.meanProbNeighborCoop
        dataArray[2] += network.meanProbNeighborDef
        dataArray[3] += network.meanProbRandomCoop
        dataArray[4] += network.meanProbRandomDef
        dataArray[5] += network.meanDegree
        dataArray[6] += network.meanProbNeighborAcc
        dataArray[7] += network.meanProbRandomAcc
        dataArray[8] += network.meanDistInclusion
        dataArray[9] += network.meanCoopFreq
        dataArray[10] += network.meanFitness
        dataArray[11] += network.meanShortestPaths
        dataArray[12] += network.meanConnComponents
        dataArray[13] += network.meanConnCompSize
        dataArray[14] += network.meanLargestConnComp
        dataArray[15] += network.meanDistConnection
    end
    dataArray[:] ./= Float64(repSims)
    return dataArray
    
    #save("sim_PNCD$(pnc)_$(pnd)_PR$(pr)_CL$(CL)_B$(BEN)_G$(gen).jld2", "parameters", [CL, BEN], "meanPNI", dataArray[1], "meanPNR", dataArray[2], "meanPR", dataArray[3], "meanDegree", dataArray[4], "meanAssortment", dataArray[5], "meanDistanceFromDefToCoop", dataArray[6], "meanDistanceInclusion", dataArray[7], "meanCooperationRatio", dataArray[8])
end

#@time begin
#    a = runSimsReturn(; B=1.0, C=0.5, D=0.0, CL=0.05, gen=10000, pn=0.5, dbOrder=mixeddb, dbProb=0.5, findMom=anyMom, neighborRange=5, distInherit=true, distFactor=0.975, pnd=false, pr=0.0001, prd=false, allowReject=true, pa=0.75, muP=0.001, delta=0.1, sigmapn=0.01, sigmapr=0.01, reps=1)
#    println(a)
#end
